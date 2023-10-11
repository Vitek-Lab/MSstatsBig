#' @keywords internal
reduceBigSpectronaut = function(input_file, output_path,
                                filter_by_excluded = FALSE,
                                filter_by_identified = FALSE,
                                filter_by_qvalue = TRUE,
                                qvalue_cutoff = 0.01) {
  if (grepl("csv", input_file)) {
    delim = ","
  } else if (grepl("tsv", input_file)) {
    delim = "\t"
  } else {
    delim = ";"
  }
  spec_chunk = function(x, pos) cleanSpectronautChunk(x,
                                                      output_path,
                                                      filter_by_excluded,
                                                      filter_by_identified,
                                                      filter_by_qvalue,
                                                      qvalue_cutoff)
  readr::read_delim_chunked(input_file,
                            readr::DataFrameCallback$new(spec_chunk),
                            delim = delim,
                            chunk_size = 1e6)
}

#' @keywords internal
cleanSpectronautChunk = function(input, output_path,
                                 filter_by_excluded = FALSE,
                                 filter_by_identified = FALSE,
                                 filter_by_qvalue = TRUE,
                                 qvalue_cutoff = 0.01) {
  all_cols = c("R.FileName", "R.Condition", "R.Replicate",
               "PG.ProteinAccessions", "EG.ModifiedSequence", "FG.LabeledSequence",
               "FG.Charge", "F.FrgIon", "F.Charge",
               "EG.Identified", "F.ExcludedFromQuantification", "F.FrgLossType",
               "PG.Qvalue", "EG.Qvalue", "F.NormalizedPeakArea")
  cols = intersect(all_cols, colnames(input))
  input = dplyr::select(input, all_of(cols))
  input = dplyr::rename_with(input, .fn = MSstatsConvert:::.standardizeColnames)

  new_names = c("Run", "Condition", "BioReplicate", "ProteinName",
                "PeptideSequence", "LabeledSequence", "PrecursorCharge", "FragmentIon",
                "ProductCharge", "Identified", "Excluded",
                "FFrgLossType", "PGQvalue", "EGQvalue",
                "Intensity")
  # non_standardized =
  old_names = MSstatsConvert:::.standardizeColnames(all_cols)
  names(old_names) = new_names
  old_names = old_names[old_names %in% colnames(input)]

  input = dplyr::rename(input, !!old_names)
  input = dplyr::mutate(input, Intensity = as.numeric(Intensity))

  if (is.character(dplyr::pull(dplyr::collect(head(dplyr::select(input, Excluded))), Excluded))) {
    input = dplyr::mutate(input, Excluded = Excluded == "True")
  }
  if (is.element("Identified", colnames(input))) {
    if (is.character(dplyr::pull(dplyr::collect(head(dplyr::select(input, Identified))), Identified))) {
      input = dplyr::mutate(input, Identified = Identified == "True")
    }
  }

  if (filter_by_excluded) {
    input = dplyr::filter(input, !Excluded)
    input = dplyr::select(input, -Excluded)
  }

  if (filter_by_identified) {
    input = dplyr::filter(input, Identified)
    input = dplyr::select(input, -Identified)
  }

  if (filter_by_qvalue) {
    input = dplyr::mutate(input, Intensity = dplyr::if_else(EGQvalue < qvalue_cutoff, Intensity, 0))
    input = dplyr::mutate(input, Intensity = dplyr::if_else(PGQvalue < qvalue_cutoff, Intensity, NA_real_))
  }

  input = dplyr::filter(input, FFrgLossType == "noloss")
  if (is.element("LabeledSequence", colnames(input))) {
    input = dplyr::mutate(input, IsLabeled = grepl("Lys8", LabeledSequence) | grepl("Arg10", LabeledSequence))
    input = dplyr::mutate(input, IsotopeLabelType := if_else(IsLabeled, "H", "L"))
  } else {
    input = dplyr::mutate(input, IsotopeLabelType = "L")
  }
  input = dplyr::select(input, ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                        ProductCharge, IsotopeLabelType, Run, BioReplicate, Condition,
                        Intensity)
  readr::write_csv(input, file = output_path)
  NULL
}
