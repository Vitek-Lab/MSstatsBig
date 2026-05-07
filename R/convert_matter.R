#' Convert preprocessed Arrow output to a matter-backed bundle
#'
#' Materializes the post-preprocessing Arrow dataset as a directory bundle
#' that MSstats can stream protein-by-protein from disk:
#'
#' \itemize{
#'   \item \code{intensities.bin} — \code{matter_list} backing file. Element
#'     \eqn{k} is one protein's \eqn{(\text{runs} \times \text{features})}
#'     matrix flattened column-major as a \code{double} vector.
#'   \item \code{run_stats.parquet} — per-run median intensity and observation
#'     count, computed in a single Arrow scan.
#'   \item \code{proteins.parquet} — \eqn{k}, ProteinName, n_features, n_runs.
#'   \item \code{features.parquet} — column ordering within each protein.
#'   \item \code{runs.parquet} — row ordering and run metadata.
#'   \item \code{manifest.rds} — canonical run order, log2 flag, schema version.
#' }
#'
#' Pipeline (no full-dataset collect at any point):
#' \enumerate{
#'   \item Lazy log2 transform applied to the Arrow plan.
#'   \item One streaming source scan → repartition by Run into Parquet files.
#'     All subsequent reads use this Parquet view.
#'   \item Per-run exact median: one small Parquet partition read per run.
#'   \item Parquet scan → distinct features per protein (small aggregate).
#'   \item Parquet scan → repartition by ProteinName into Parquet files.
#'   \item For each protein: read its (small) partition, scatter into matrix,
#'     write that protein's slot of the matter_list.
#' }
#'
#' Assumes (a) one (Run, Feature) pair has a single intensity, and
#' (b) the experiment has no fractions.
#'
#' @param arrow_data arrow Dataset / lazy reference, output of
#'   \code{MSstatsPreprocessBig}.
#' @param output_dir directory to write the bundle. Created if it does not exist.
#' @param log2_transform logical; if TRUE (default), take \eqn{\log_2}
#'   intensity, mapping non-positive values to NA, before any downstream stats
#'   or matrix writes.
#' @param verbose logical; print progress messages.
#' @param progress_every integer; emit a progress message every N proteins
#'   during the per-protein write loop.
#'
#' @return invisibly, a named list of file paths in the bundle.
#'
#' @keywords internal
.convertArrowToMatterList <- function(arrow_data,
                                      output_dir,
                                      log2_transform = TRUE,
                                      verbose        = TRUE,
                                      progress_every = 1000L) {
    if (!requireNamespace("matter", quietly = TRUE))
        stop("Package 'matter' is required for convert_matter = TRUE.")
    if (!requireNamespace("arrow", quietly = TRUE))
        stop("Package 'arrow' is required for convert_matter = TRUE.")

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # ── 1. Lazy log2 transform ────────────────────────────────────────────────
    if (log2_transform) {
        arrow_data <- dplyr::mutate(
            arrow_data,
            Intensity = dplyr::if_else(Intensity > 0,
                                       log2(Intensity),
                                       NA_real_))
    }

    # ── 1b. Materialize source as Run-partitioned Parquet (one streaming scan)
    #
    # All subsequent reads (per-run median loop, distinct features, ProteinName
    # repartition) go through this Parquet view instead of re-scanning the
    # source CSV.  One Run = one small file → predicate-pushdown filters
    # become ~O(file size for that run), not O(full dataset).
    if (verbose)
        message("[convert_matter] partitioning source by Run (one streaming scan) ...")

    run_partitioned <- file.path(output_dir, ".tmp_run_partitions")
    if (dir.exists(run_partitioned))
        unlink(run_partitioned, recursive = TRUE)

    arrow::write_dataset(
        arrow_data,
        path         = run_partitioned,
        format       = "parquet",
        partitioning = "Run"
    )

    # Replace lazy CSV reference with a lazy Parquet reference for everything
    # downstream.  Hive partitioning auto-recovers Run as a virtual column.
    arrow_data <- arrow::open_dataset(run_partitioned)

    # ── 2. Per-run exact median + run metadata ────────────────────────────────
    #
    # Arrow's groupby `median` translates to `hash_approximate_median`
    # (t-digest), so we pull each run's Intensity values one at a time and
    # compute an EXACT median in R.  Each filter+collect now reads only that
    # run's Parquet partition file.
    run_meta <- arrow_data |>
        dplyr::distinct(Run, Condition, BioReplicate) |>
        dplyr::collect() |>
        dplyr::arrange(Run)

    distinct_runs <- run_meta$Run
    n_runs_total  <- length(distinct_runs)

    if (verbose)
        message(sprintf("[convert_matter] computing exact per-run median across %d runs ...",
                        n_runs_total))

    run_medians_list <- vector("list", n_runs_total)
    for (i in seq_len(n_runs_total)) {
        g    <- distinct_runs[i]
        vals <- arrow_data |>
            dplyr::filter(Run == g) |>
            dplyr::select(Intensity) |>
            dplyr::collect()

        run_medians_list[[i]] <- data.table::data.table(
            Run              = g,
            median_intensity = median(vals$Intensity, na.rm = TRUE),
            n_obs            = sum(!is.na(vals$Intensity)))

        rm(vals)

        if (verbose && i %% 10L == 0L)
            message(sprintf("[convert_matter]   median computed for %d / %d runs",
                            i, n_runs_total))
    }
    run_medians <- data.table::rbindlist(run_medians_list)

    run_meta$row_idx <- seq_len(nrow(run_meta))
    all_runs <- run_meta$Run
    n_runs   <- length(all_runs)

    run_stats <- dplyr::left_join(run_meta, run_medians, by = "Run")
    arrow::write_parquet(run_stats, file.path(output_dir, "run_stats.parquet"))
    arrow::write_parquet(run_meta,  file.path(output_dir, "runs.parquet"))

    # ── 3. Feature catalog per protein (streaming distinct) ───────────────────
    if (verbose) message("[convert_matter] enumerating features per protein ...")

    feature_meta <- arrow_data |>
        dplyr::mutate(Feature = paste(PeptideSequence, PrecursorCharge,
                                       FragmentIon, ProductCharge,
                                       IsotopeLabelType, sep = "_")) |>
        dplyr::distinct(ProteinName, Feature, PeptideSequence,
                        PrecursorCharge, FragmentIon, ProductCharge,
                        IsotopeLabelType) |>
        dplyr::collect() |>
        dplyr::arrange(ProteinName, Feature)

    feature_meta <- feature_meta |>
        dplyr::group_by(ProteinName) |>
        dplyr::mutate(col_idx = seq_along(Feature)) |>
        dplyr::ungroup()

    arrow::write_parquet(feature_meta,
                         file.path(output_dir, "features.parquet"))

    proteins_meta <- feature_meta |>
        dplyr::group_by(ProteinName) |>
        dplyr::summarise(n_features = dplyr::n(), .groups = "drop") |>
        dplyr::arrange(ProteinName) |>
        dplyr::mutate(k = seq_along(ProteinName), n_runs = n_runs)

    arrow::write_parquet(proteins_meta,
                         file.path(output_dir, "proteins.parquet"))

    proteins    <- proteins_meta$ProteinName
    n_proteins  <- length(proteins)
    lengths_vec <- as.numeric(proteins_meta$n_features) * n_runs

    if (verbose)
        message(sprintf(
            "[convert_matter] %d proteins, %d runs, %.1fM cells (%.1f MB on disk)",
            n_proteins, n_runs,
            sum(lengths_vec) / 1e6,
            sum(lengths_vec) * 8 / 1024^2))

    # ── 4. Repartition by ProteinName (Arrow streaming write) ─────────────────
    #
    # Project to the 4 columns we need, then let Arrow stream the source CSV
    # in batches and write one Parquet file per protein.  Peak RAM = one
    # Arrow batch (typically a few hundred MB), not the full dataset.
    if (verbose)
        message("[convert_matter] repartitioning by ProteinName (streaming) ...")

    tmp_partitioned <- file.path(output_dir, ".tmp_partitions")
    if (dir.exists(tmp_partitioned))
        unlink(tmp_partitioned, recursive = TRUE)

    arrow::write_dataset(
        arrow_data |>
            dplyr::mutate(Feature = paste(PeptideSequence, PrecursorCharge,
                                           FragmentIon, ProductCharge,
                                           IsotopeLabelType, sep = "_")) |>
            dplyr::select(ProteinName, Feature, Run, Intensity),
        path         = tmp_partitioned,
        format       = "parquet",
        partitioning = "ProteinName"
    )

    # Re-open the partitioned dataset for predicate-pushdown reads per protein.
    partitioned_ds <- arrow::open_dataset(tmp_partitioned)

    # ── 5. Pre-allocate matter_list (empty backing file, sized slots) ────────
    #
    # matter_list with NULL data + per-element `lengths` creates the backing
    # file with one slot per protein at the right offsets.  Subsequent
    # `intensities[[i]] <- vec` writes directly to disk via matter's C API.
    if (verbose) message("[convert_matter] pre-allocating matter_list ...")

    backing_file <- file.path(output_dir, "intensities.bin")
    intensities <- matter::matter_list(
        NULL,
        type    = "double",
        path    = backing_file,
        lengths = as.integer(lengths_vec),
        names   = proteins
    )

    # ── 6. Stream partitions, write each protein's slot directly ─────────────
    #
    # Peak RAM = one protein's matrix.  No accumulating list — every
    # iteration drops the previous protein's data after writing.
    if (verbose) message("[convert_matter] streaming per-protein matrices to disk ...")

    feat_by_prot <- split(feature_meta$Feature, feature_meta$ProteinName)

    for (i in seq_len(n_proteins)) {
        prot       <- proteins[i]
        prot_feats <- feat_by_prot[[prot]]

        # Single-protein scan via predicate pushdown.
        prot_data <- partitioned_ds |>
            dplyr::filter(ProteinName == prot) |>
            dplyr::select(Run, Feature, Intensity) |>
            dplyr::collect() |>
            data.table::as.data.table()

        run_idx  <- match(prot_data$Run,     all_runs)
        feat_idx <- match(prot_data$Feature, prot_feats)
        valid    <- !is.na(run_idx) & !is.na(feat_idx)

        mat <- matrix(NA_real_, nrow = n_runs, ncol = length(prot_feats))
        if (any(valid))
            mat[cbind(run_idx[valid], feat_idx[valid])] <-
                prot_data$Intensity[valid]

        intensities[[i]] <- as.numeric(mat)   # writes directly to backing file

        rm(prot_data, mat)

        if (verbose && i %% progress_every == 0L)
            message(sprintf("[convert_matter]   wrote %d / %d proteins to disk",
                            i, n_proteins))
    }

    unlink(tmp_partitioned, recursive = TRUE)
    unlink(run_partitioned, recursive = TRUE)

    # ── 7. Manifest ───────────────────────────────────────────────────────────
    manifest <- list(
        schema_version      = "1.0",
        created             = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        log2_intensity      = log2_transform,
        n_proteins          = n_proteins,
        n_runs              = n_runs,
        canonical_run_order = all_runs,
        files = list(
            intensities = "intensities.bin",
            run_stats   = "run_stats.parquet",
            proteins    = "proteins.parquet",
            features    = "features.parquet",
            runs        = "runs.parquet"
        )
    )
    saveRDS(manifest, file.path(output_dir, "manifest.rds"))

    if (verbose)
        message(sprintf("[convert_matter] bundle written to %s",
                        normalizePath(output_dir)))

    invisible(list(
        output_dir  = normalizePath(output_dir),
        intensities = backing_file,
        run_stats   = file.path(output_dir, "run_stats.parquet"),
        proteins    = file.path(output_dir, "proteins.parquet"),
        features    = file.path(output_dir, "features.parquet"),
        runs        = file.path(output_dir, "runs.parquet"),
        manifest    = file.path(output_dir, "manifest.rds")
    ))
}
