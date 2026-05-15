#' @keywords internal
reduceBigSpectronaut <- function(input_file, output_path,
                                 intensity="F.NormalizedPeakArea",
                                 filter_by_excluded = FALSE,
                                 filter_by_identified = FALSE,
                                 filter_by_qvalue = TRUE,
                                 qvalue_cutoff = 0.01,
                                 calculateAnomalyScores=FALSE,
                                 anomalyModelFeatures=c(),
                                 block_size = 16L * 1024L * 1024L) {
  block_size <- as.integer(block_size)
  stopifnot(length(block_size) == 1L, !is.na(block_size), block_size > 0L)

  if (grepl("csv", input_file)) {
    delim <- ","
  } else if (grepl("tsv|xls", input_file)) {
    delim <- "\t"
  } else {
    delim <- ";"
  }

  # Columns cleanSpectronautChunk actually consumes; Arrow's
  # convert_options$include_columns drops everything else at parse time so
  # we never materialize the ~35 unused columns Spectronaut exports.
  needed_cols <- c("R.FileName", "R.Condition", "R.Replicate",
                   "PG.ProteinAccessions", "EG.ModifiedSequence",
                   "FG.LabeledSequence", "FG.Charge",
                   "F.FrgIon", "F.Charge",
                   "EG.Identified", "F.ExcludedFromQuantification",
                   "F.FrgLossType", "PG.Qvalue", "EG.Qvalue",
                   intensity)
  if (calculateAnomalyScores) {
    needed_cols <- c(needed_cols, anomalyModelFeatures)
  }

  # Arrow's CSV reader replaces readr::read_delim_chunked.  Arrow releases
  # per-batch state as soon as a batch is consumed, so peak memory is
  # bounded by one record batch instead of growing with the dataset (readr
  # keeps a string-interning pool that accumulates across chunks).  The
  # `delim` switch above already covers comma / tab / semicolon variants;
  # Arrow's CSV reader handles all three the same way through
  # CsvParseOptions$delimiter.
  parse_opts   <- arrow::CsvParseOptions$create(delimiter = delim)
  convert_opts <- arrow::CsvConvertOptions$create()
  read_opts    <- arrow::CsvReadOptions$create(block_size = block_size)

  ds <- arrow::open_dataset(
    input_file,
    format          = "csv",
    parse_options   = parse_opts,
    convert_options = convert_opts,
    read_options    = read_opts
  )

  reader <- arrow::Scanner$create(ds)$ToRecordBatchReader()

  t_start   <- Sys.time()
  pos       <- 1L
  batch_idx <- 0L
  repeat {
    batch <- reader$read_next_batch()
    if (is.null(batch)) break
    chunk_df <- as.data.frame(batch)
    cleanSpectronautChunk(chunk_df,
                          output_path,
                          intensity,
                          filter_by_excluded,
                          filter_by_identified,
                          filter_by_qvalue,
                          qvalue_cutoff,
                          pos,
                          calculateAnomalyScores,
                          anomalyModelFeatures)
    pos       <- pos + nrow(chunk_df)
    batch_idx <- batch_idx + 1L

    if (batch_idx %% 1000L == 0L) {
      elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
      rate    <- (pos - 1L) / elapsed
      message(sprintf(
        "[reduceBigSpectronaut] %d batches | %s rows | %.1fk rows/s | %.0fs elapsed",
        batch_idx,
        format(pos - 1L, big.mark = ","),
        rate / 1000,
        elapsed))
    }

    rm(batch, chunk_df)
  }

  if (batch_idx %% 1000L != 0L) {
    elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
    rate    <- (pos - 1L) / elapsed
    message(sprintf(
      "[reduceBigSpectronaut] done: %d batches | %s rows | %.1fk rows/s | %.0fs elapsed",
      batch_idx,
      format(pos - 1L, big.mark = ","),
      rate / 1000,
      elapsed))
  }
}

#' @keywords internal
cleanSpectronautChunk = function(input, output_path,
                                 intensity="F.NormalizedPeakArea",
                                 filter_by_excluded = FALSE,
                                 filter_by_identified = FALSE,
                                 filter_by_qvalue = TRUE,
                                 qvalue_cutoff = 0.01,
                                 pos = NULL,
                                 calculateAnomalyScores=FALSE, 
                                 anomalyModelFeatures=c()) {
  all_cols <- c("R.FileName", "R.Condition", "R.Replicate",
                "PG.ProteinAccessions", "EG.ModifiedSequence", "FG.LabeledSequence",
                "FG.Charge", "F.FrgIon", "F.Charge",
                "EG.Identified", "F.ExcludedFromQuantification", "F.FrgLossType",
                "PG.Qvalue", "EG.Qvalue", intensity)
  
  if (calculateAnomalyScores){
    all_cols <- c(all_cols, anomalyModelFeatures)
  }
  
  cols <- intersect(all_cols, colnames(input))
  input <- dplyr::select(input, all_of(cols))
  input <- dplyr::rename_with(input, .fn = MSstatsConvert:::.standardizeColnames)
  
  new_names <- c("Run", "Condition", "BioReplicate", "ProteinName",
                 "PeptideSequence", "LabeledSequence", "PrecursorCharge", "FragmentIon",
                 "ProductCharge", "Identified", "Excluded",
                 "FFrgLossType", "PGQvalue", "EGQvalue",
                 "Intensity")
  if (calculateAnomalyScores){
    new_names <- c(new_names, MSstatsConvert:::.standardizeColnames(anomalyModelFeatures))
  }
  
  # non_standardized =
  old_names <- MSstatsConvert:::.standardizeColnames(all_cols)
  names(old_names) <- new_names
  old_names <- old_names[old_names %in% colnames(input)]
  
  input <- dplyr::rename(input, !!old_names)
  input <- dplyr::mutate(input, Intensity = as.numeric(Intensity))
  
  if (is.character(dplyr::pull(dplyr::collect(head(dplyr::select(input, Excluded))), Excluded))) {
    input <- dplyr::mutate(input, Excluded = Excluded == "True")
  }
  if (is.element("Identified", colnames(input))) {
    if (is.character(dplyr::pull(dplyr::collect(head(dplyr::select(input, Identified))), Identified))) {
      input <- dplyr::mutate(input, Identified = Identified == "True")
    }
  }
  
  if (filter_by_excluded) {
    input <- dplyr::mutate(
      input, Intensity = dplyr::if_else(Excluded, NA_real_, Intensity))
    
  }
  
  if (filter_by_identified) {
    input <- dplyr::mutate(
      input, Intensity = dplyr::if_else(Identified, Intensity, NA_real_))
  }
  
  if (filter_by_qvalue) {
    input <- dplyr::mutate(
      input,
      Intensity = dplyr::if_else(EGQvalue < qvalue_cutoff, Intensity, NA_real_))
    input <- dplyr::mutate(
      input, 
      Intensity = dplyr::if_else(PGQvalue < qvalue_cutoff, Intensity, NA_real_))
  }
  
  input <- dplyr::filter(input, FFrgLossType == "noloss")
  if (is.element("LabeledSequence", colnames(input))) {
    input <- dplyr::mutate(input, IsLabeled = grepl("Lys8", LabeledSequence) | grepl("Arg10", LabeledSequence))
    input <- dplyr::mutate(input, IsotopeLabelType := dplyr::if_else(IsLabeled, "H", "L"))
  } else {
    input <- dplyr::mutate(input, IsotopeLabelType = "L")
  }
  
  select_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                  "ProductCharge", "IsotopeLabelType", "Run", "BioReplicate", "Condition",
                  "Intensity")
  if (calculateAnomalyScores){
    select_cols = c(select_cols, 
                    MSstatsConvert:::.standardizeColnames(anomalyModelFeatures))
  }
  
  input <- dplyr::select(input, select_cols)
  .writeChunkToFile(input, output_path, pos)
  NULL
}
