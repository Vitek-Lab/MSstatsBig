#' @keywords internal
reduceBigSpectronaut <- function(input_file, output_path,
                                 intensity="F.NormalizedPeakArea",
                                 filter_by_excluded = FALSE,
                                 filter_by_identified = FALSE,
                                 filter_by_qvalue = TRUE,
                                 qvalue_cutoff = 0.01,
                                 calculateAnomalyScores=FALSE,
                                 anomalyModelFeatures=c()) {
  if (grepl("csv", input_file)) {
    delim = ","
  } else if (grepl("tsv|xls", input_file)) {
    delim = "\t"
  } else {
    delim <- ";"
  }

  # Restrict parsing to the columns cleanSpectronautChunk actually consumes.
  # Spectronaut exports often have 50+ columns; reading only this subset
  # cuts per-chunk peak memory roughly proportionally to the column ratio.
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

  spec_chunk <- function(x, pos) cleanSpectronautChunk(x,
                                                       output_path,
                                                       intensity,
                                                       filter_by_excluded,
                                                       filter_by_identified,
                                                       filter_by_qvalue,
                                                       qvalue_cutoff,
                                                       pos,
                                                       calculateAnomalyScores,
                                                       anomalyModelFeatures)
  readr::read_delim_chunked(input_file,
                            readr::DataFrameCallback$new(spec_chunk),
                            delim = delim,
                            chunk_size = 1e6,
                            col_names = tidyselect::any_of(needed_cols))
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
