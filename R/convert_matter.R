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
#'   \item Streaming Arrow scan → per-run median statistics (small aggregate).
#'   \item Streaming Arrow scan → distinct features per protein (small aggregate).
#'   \item Streaming Arrow scan → repartition by ProteinName into Parquet files.
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

    # ── 2. Per-run median + run metadata (streaming aggregates) ───────────────
    if (verbose) message("[convert_matter] computing per-run median ...")

    run_medians <- arrow_data |>
        dplyr::group_by(Run) |>
        dplyr::summarise(
            median_intensity = median(Intensity, na.rm = TRUE),
            n_obs            = sum(!is.na(Intensity), na.rm = TRUE),
            .groups          = "drop") |>
        dplyr::collect()

    run_meta <- arrow_data |>
        dplyr::distinct(Run, Condition, BioReplicate) |>
        dplyr::collect() |>
        dplyr::arrange(Run)

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

    # ── 5. Stream partitions, build packed_list one protein at a time ────────
    #
    # Reads are streamed (one protein's partition per iteration); only the
    # accumulated packed_list is held in RAM, which equals the eventual
    # on-disk size of intensities.bin (numerics only, no string overhead).
    if (verbose) message("[convert_matter] packing per-protein matrices ...")

    feat_by_prot <- split(feature_meta$Feature, feature_meta$ProteinName)

    packed_list <- vector("list", n_proteins)
    for (i in seq_len(n_proteins)) {
        prot       <- proteins[i]
        prot_feats <- feat_by_prot[[prot]]

        # Single-protein scan via predicate pushdown — only this protein's
        # partition file is read.
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

        packed_list[[i]] <- as.numeric(mat)   # column-major flatten

        if (verbose && i %% progress_every == 0L)
            message(sprintf("[convert_matter]   packed %d / %d proteins",
                            i, n_proteins))
    }

    # ── 6. Materialize matter_list and clean up partition staging dir ────────
    if (verbose) message("[convert_matter] writing matter_list backing file ...")

    backing_file <- file.path(output_dir, "intensities.bin")
    intensities <- matter::matter_list(
        packed_list,
        type  = "double",
        path  = backing_file,
        names = proteins
    )

    rm(packed_list); gc(verbose = FALSE)

    unlink(tmp_partitioned, recursive = TRUE)

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
