#' @keywords internal
reduceBigSpectronaut = function(input_file, output_path,
                               intensity = "F.PeakArea",
                               filter_by_excluded = TRUE,
                               filter_by_identified = TRUE,
                               filter_by_qvalue = TRUE,
                               qvalue_cutoff = 0.01) {
  input = arrow::open_dataset(input_file, format = "tsv")
  cols = c("R.FileName", "R.Condition", "R.Replicate",
           "PG.ProteinAccessions", "EG.ModifiedSequence", "FG.LabeledSequence",
           "FG.Charge", "F.FrgIon", "F.Charge",
           "EG.Identified", "F.ExcludedFromQuantification", "F.FrgLossType",
           "PG.Qvalue", "EG.Qvalue", "F.NormalizedPeakArea", "F.MeasuredRelativeIntensity",
           "F.PeakArea",
           "F.MassAccuracyPPM", "FG.FWHM", "EG.ApexRT", "FG.ShapeQualityScore")
  input = dplyr::select(input, all_of(cols))
  input = dplyr::rename_with(input, .fn = MSstatsConvert:::.standardizeColnames)

  new_names = c("Run", "Condition", "BioReplicate", "ProteinName",
                "PeptideSequence", "LabeledSequence", "PrecursorCharge", "FragmentIon",
                "ProductCharge", "Identified", "Excluded",
                "FFrgLossType", "PGQvalue", "EGQvalue",
                "Intensity", "MeasuredRelativeIntensity", "PeakArea",
                "MassAccuracyPPM", "FWHM", "ApexRt", "ShapeQualityScore")
  old_names = MSstatsConvert:::.standardizeColnames(cols)
  names(old_names) = new_names

  input = dplyr::rename(input, !!old_names)
  input = dplyr::mutate(input, Intensity = as.numeric(Intensity))

  if (is.character(dplyr::pull(dplyr::collect(head(dplyr::select(input, Excluded))), Excluded))) {
    input = dplyr::mutate(input, Excluded = Excluded == "True")
  }
  if (is.character(dplyr::pull(dplyr::collect(head(dplyr::select(input, Identified))), Identified))) {
    input = dplyr::mutate(input, Identified = Identified == "True")
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
    input = dplyr::mutate(input, Intensity = if_else(EGQvalue < qvalue_cutoff, Intensity, 0))
    input = dplyr::mutate(input, Intensity = if_else(PGQvalue < qvalue_cutoff, Intensity, NA_real_))
  }

  input = dplyr::filter(input, FFrgLossType == "noloss")
  input = dplyr::mutate(input, IsLabeled = grepl("Lys8", LabeledSequence) | grepl("Arg10", LabeledSequence))
  input = dplyr::mutate(input, IsotopeLabelType := if_else(IsLabeled, "H", "L"))
  input = dplyr::select(input, ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                        ProductCharge, IsotopeLabelType, Run, BioReplicate, Condition,
                        Intensity, EGQvalue, PGQvalue)
  arrow::write_csv_arrow(input, file = output_path)
  TRUE
}
