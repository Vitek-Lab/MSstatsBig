#' Data processing for matter-backed bundles produced by the big converters
#'
#' A streaming, parallel analogue of \code{MSstats::dataProcess} that reads
#' from the bundle written by \code{bigSpectronauttoMSstatsFormat(convert_matter
#' = TRUE)} and emits a single ProteinLevelData file. Per-protein work is
#' dispatched via \code{BiocParallel::bplapply} with \code{matter::SnowfastParam};
#' workers stream their assigned protein slot directly from the
#' \code{intensities.bin} backing file rather than receiving data through
#' sockets.
#'
#' Pipeline per protein (inside each worker):
#' \enumerate{
#'   \item Load \code{intensities[[k]]} → \eqn{(\text{runs} \times
#'     \text{features})} matrix.
#'   \item Equalize-medians normalization: subtract pre-computed per-run
#'     median, add the global median (a single vectorized sweep). Stats come
#'     from \code{run_stats.parquet} produced at conversion time.
#'   \item Top-N feature selection: rank features within the protein by mean
#'     post-normalization abundance, keep top \code{n_top_feature}.
#'   \item Build the V3 packed double vector + meta list (FL × R layout) and
#'     call \code{MSstatsSummarizeSingleTMPV2}, which performs imputation
#'     (AFT survival) and Tukey median-polish summarization in matrix form.
#'   \item Augment per-run summary with \code{NumMeasuredFeature},
#'     \code{MissingPercentage}, \code{more50missing}, \code{NumImputedFeature}.
#' }
#'
#' Assumes (a) one (Run, Feature) value per cell, (b) no fractions, and
#' (c) unlabeled experiment (single \code{IsotopeLabelType}).
#'
#' @param matter_dir path to the bundle directory written by the converter
#'   (must contain intensities.bin, manifest.rds, proteins.parquet,
#'   features.parquet, runs.parquet, run_stats.parquet).
#' @param output_file path to write the assembled ProteinLevelData. Extension
#'   determines format: \code{.rds} (default), \code{.parquet}, or \code{.csv}.
#'   Pass \code{NULL} to skip writing and only return.
#' @param n_top_feature integer; max features per protein after normalization.
#' @param impute logical; AFT-survival imputation of NAs treated as censored.
#' @param censored_symbol character; \code{"NA"} (default) treats NAs as
#'   censored at the detection limit; \code{"0"} treats zeros (already mapped
#'   to NA at conversion if log2_transform=TRUE); \code{NULL} disables.
#' @param remove50missing logical; drop proteins where every run is >50\%
#'   missing.
#' @param aft_iterations integer; max AFT iterations (passed through).
#' @param numberOfCores integer; ignored when \code{BPPARAM} is supplied.
#' @param BPPARAM optional \code{BiocParallel::BiocParallelParam}; defaults to
#'   \code{SnowfastParam(numberOfCores)} or \code{SerialParam()} when
#'   \code{numberOfCores <= 1}.
#' @param verbose logical; print progress.
#'
#' @return invisibly, the assembled ProteinLevelData as a \code{data.table}.
#'
#' @importFrom matter matter_list SnowfastParam
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom data.table data.table as.data.table rbindlist setnames
#' @importFrom arrow read_parquet write_parquet
#' @importFrom stats median
#'
#' @export
dataProcessBig <- function(matter_dir,
                           output_file     = "ProteinLevelData.rds",
                           n_top_feature   = 100L,
                           impute          = TRUE,
                           censored_symbol = "NA",
                           remove50missing = FALSE,
                           aft_iterations  = 90L,
                           numberOfCores   = 1L,
                           BPPARAM         = NULL,
                           verbose         = TRUE) {
    if (!requireNamespace("matter",       quietly = TRUE))
        stop("Package 'matter' is required.")
    if (!requireNamespace("BiocParallel", quietly = TRUE))
        stop("Package 'BiocParallel' is required.")
    if (!requireNamespace("arrow",        quietly = TRUE))
        stop("Package 'arrow' is required.")

    matter_dir <- normalizePath(matter_dir, mustWork = TRUE)

    # ── 1. Load bundle metadata ───────────────────────────────────────────────
    if (verbose) message("[dataProcessBig] reading bundle metadata ...")

    manifest      <- readRDS(file.path(matter_dir, "manifest.rds"))
    proteins_meta <- arrow::read_parquet(
        file.path(matter_dir, "proteins.parquet"))
    features_meta <- arrow::read_parquet(
        file.path(matter_dir, "features.parquet"))
    runs_meta     <- arrow::read_parquet(
        file.path(matter_dir, "runs.parquet"))
    run_stats     <- arrow::read_parquet(
        file.path(matter_dir, "run_stats.parquet"))

    n_proteins <- nrow(proteins_meta)
    n_runs     <- as.integer(manifest$n_runs)
    all_runs   <- as.character(manifest$canonical_run_order)

    backing_file <- file.path(matter_dir, "intensities.bin")
    intensities  <- matter::matter_list(
        NULL,
        type     = "double",
        path     = backing_file,
        lengths  = as.integer(proteins_meta$n_features) * n_runs,
        names    = as.character(proteins_meta$ProteinName),
        readonly = TRUE)

    # ── 2. Pre-compute equalize-medians shifts (global view, vectorized) ─────
    run_median_by_id <- run_stats$median_intensity[
        match(all_runs, run_stats$Run)]
    global_median    <- median(run_median_by_id, na.rm = TRUE)
    norm_shift       <- run_median_by_id - global_median
    norm_shift[is.na(norm_shift)] <- 0

    # ── 3. Pre-split feature metadata into a per-protein lookup ──────────────
    features_by_prot <- split(
        as.data.frame(features_meta[, c("Feature", "PeptideSequence",
                                         "IsotopeLabelType")]),
        features_meta$ProteinName)

    # ── 4. Worker closure ────────────────────────────────────────────────────
    .worker <- local({
        intensities_      <- intensities
        proteins_meta_    <- proteins_meta
        features_by_prot_ <- features_by_prot
        all_runs_         <- all_runs
        norm_shift_       <- norm_shift
        n_runs_           <- n_runs
        n_top_feature_    <- as.integer(n_top_feature)
        impute_           <- impute
        censored_symbol_  <- censored_symbol
        remove50missing_  <- remove50missing
        aft_iterations_   <- as.integer(aft_iterations)

        function(k) {
            library(MSstats, quietly = TRUE, warn.conflicts = FALSE)
            library(matter,  quietly = TRUE, warn.conflicts = FALSE)

            prot_id    <- as.character(proteins_meta_$ProteinName[k])
            n_feat_in  <- as.integer(proteins_meta_$n_features[k])
            feat_meta  <- features_by_prot_[[prot_id]]

            # Load (n_runs × n_features) matrix from disk.
            mat <- matrix(intensities_[[k]],
                          nrow = n_runs_, ncol = n_feat_in)

            # ── 4a. Equalize-medians normalization ────────────────────────────
            mat <- mat - norm_shift_   # NAs propagate, recycled across cols

            # ── 4b. Top-N feature selection (post-normalization) ──────────────
            keep <- if (n_feat_in > n_top_feature_) {
                feat_means <- colMeans(mat, na.rm = TRUE)
                feat_means[is.nan(feat_means)] <- -Inf
                sort(order(feat_means, decreasing = TRUE)[
                    seq_len(n_top_feature_)])
            } else {
                seq_len(n_feat_in)
            }
            mat       <- mat[, keep, drop = FALSE]
            feat_meta <- feat_meta[keep, , drop = FALSE]
            FL        <- ncol(mat)
            R         <- n_runs_

            # ── 4c. Build V3 packed format (FL × R, column-major) ─────────────
            mat_FxR <- t(mat)                      # n_features × n_runs
            cens    <- as.numeric(is.na(mat_FxR))  # 1 = censored, 0 = observed
            cen     <- 1 - cens                    # event indicator

            n_obs_per_feat <- rowSums(!is.na(mat_FxR))
            n_obs_per_run  <- colSums(!is.na(mat_FxR))
            prop_features  <- if (FL == 0L) rep(0, R) else n_obs_per_run / FL

            packed <- c(
                as.double(k),
                as.vector(mat_FxR),                 # newABUNDANCE
                rep(NA_real_, FL * R),               # ABUNDANCE (TMP unused)
                as.vector(cens),                     # censored
                as.vector(cen),                      # cen
                rep(NA_real_, FL * R),               # ANOMALYSCORES (unused)
                as.double(n_obs_per_feat),
                as.double(n_obs_per_run),
                as.double(prop_features)
            )

            flp <- data.frame(
                FEATURE = as.character(feat_meta$Feature),
                LABEL   = as.character(feat_meta$IsotopeLabelType),
                PEPTIDE = as.character(feat_meta$PeptideSequence),
                stringsAsFactors = FALSE)

            meta <- list(
                PROTEIN           = prot_id,
                feat_label_pep    = flp,
                runs              = all_runs_,
                FL                = FL,
                R                 = R,
                is_labeled_ref    = FALSE,
                has_ABUNDANCE     = FALSE,
                has_censored      = any(cens > 0.5),
                has_cen           = TRUE,
                has_anom          = FALSE,
                add_ref_covariate = FALSE)

            # ── 4d. Summarize ────────────────────────────────────────────────
            res <- MSstats:::MSstatsSummarizeSingleTMPV2(
                packed, meta,
                impute          = impute_,
                censored_symbol = censored_symbol_,
                remove50missing = remove50missing_,
                aft_iterations  = aft_iterations_)

            result_dt <- res[[1L]]
            if (is.null(result_dt) || nrow(result_dt) == 0L)
                return(NULL)

            result_dt <- data.table::as.data.table(result_dt)

            # ── 4e. Per-run measurement stats (pre-imputation) ────────────────
            run_idx <- match(as.character(result_dt$RUN), all_runs_)
            num_meas       <- n_obs_per_run[run_idx]
            missing_pct    <- 1 - num_meas / FL
            num_imputed    <- if (impute_)
                                  pmax(FL - num_meas, 0L)
                              else
                                  rep(0L, length(num_meas))

            result_dt[, NumMeasuredFeature := num_meas]
            result_dt[, MissingPercentage  := missing_pct]
            result_dt[, more50missing      := missing_pct > 0.5]
            result_dt[, NumImputedFeature  := num_imputed]

            result_dt
        }
    })

    # ── 5. Dispatch ──────────────────────────────────────────────────────────
    if (is.null(BPPARAM)) {
        BPPARAM <- if (numberOfCores > 1L) {
            matter::SnowfastParam(workers       = numberOfCores,
                                  force.GC      = TRUE,
                                  stop.on.error = FALSE)
        } else {
            BiocParallel::SerialParam()
        }
    }

    if (verbose)
        message(sprintf("[dataProcessBig] processing %d proteins ...",
                        n_proteins))

    results <- BiocParallel::bplapply(
        seq_len(n_proteins), .worker, BPPARAM = BPPARAM)

    # ── 6. Assemble ProteinLevelData ──────────────────────────────────────────
    if (verbose) message("[dataProcessBig] assembling ProteinLevelData ...")

    plv <- data.table::rbindlist(results, fill = TRUE)
    if (nrow(plv) == 0L) {
        warning("[dataProcessBig] no protein produced summarization results.")
        return(invisible(plv))
    }

    runs_dt <- data.table::as.data.table(runs_meta)
    plv <- merge(plv, runs_dt,
                 by.x = "RUN", by.y = "Run", all.x = TRUE)

    data.table::setnames(plv,
        old = c("RUN",         "Condition", "BioReplicate"),
        new = c("originalRUN", "GROUP",     "SUBJECT"),
        skip_absent = TRUE)

    plv[, RUN := originalRUN]

    plv[, TotalGroupMeasurements := sum(NumMeasuredFeature, na.rm = TRUE),
        by = c("Protein", "GROUP", "LABEL")]

    out_cols <- c("Protein", "LABEL", "RUN", "originalRUN",
                  "LogIntensities", "GROUP", "SUBJECT",
                  "TotalGroupMeasurements", "NumMeasuredFeature",
                  "MissingPercentage", "more50missing", "NumImputedFeature")
    out_cols <- intersect(out_cols, colnames(plv))
    plv <- plv[, out_cols, with = FALSE]

    # ── 7. Save ───────────────────────────────────────────────────────────────
    if (!is.null(output_file)) {
        ext <- tolower(tools::file_ext(output_file))
        if (ext == "parquet" || ext == "pq") {
            arrow::write_parquet(plv, output_file)
        } else if (ext == "csv") {
            data.table::fwrite(plv, output_file)
        } else {
            saveRDS(plv, output_file)
        }
        if (verbose)
            message(sprintf("[dataProcessBig] wrote %s (%d rows)",
                            output_file, nrow(plv)))
    }

    invisible(plv)
}
