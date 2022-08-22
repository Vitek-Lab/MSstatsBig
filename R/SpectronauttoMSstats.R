#' @export
cleanBigSpectronautArrow = function(input_file, output_path,
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


#' @export
#'
BigSpectronauttoMSstatsFormat = function(input_file, output_path) {
  input = arrow::open_dataset(input_file, format = "csv")
  input = dplyr::select(input, ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                        ProductCharge, IsotopeLabelType, Run, BioReplicate, Condition, Intensity)
  input = dplyr::select(input, -BioReplicate, -Condition)
  input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                          ProductCharge, IsotopeLabelType)
  observation_counts = dplyr::summarize(input, NumObs = sum(!is.na(Intensity) & Intensity > 0))
  observation_counts = dplyr::filter(observation_counts, NumObs > 1)
  observation_counts = dplyr::select(observation_counts, -NumObs)
  input = dplyr::inner_join(input, observation_counts,
                            by = c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                                   "ProductCharge", "IsotopeLabelType"))

  pp_df = dplyr::select(input, ProteinName, PeptideSequence)
  pp_df = dplyr::group_by(pp_df, PeptideSequence)
  pp_df = dplyr::summarize(pp_df, NumProteins = n_distinct(ProteinName))
  pp_df = dplyr::filter(pp_df, NumProteins == 1)
  pp_df = dplyr::select(pp_df, -NumProteins)
  input = dplyr::inner_join(input, pp_df, by = "PeptideSequence")

  input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge,
                          FragmentIon, ProductCharge, IsotopeLabelType, Run)
  input = dplyr::summarize(input, Intensity = max(Intensity))

  input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                          ProductCharge, IsotopeLabelType)
  observation_counts = dplyr::summarize(input, NumObs = sum(!is.na(Intensity) & Intensity > 0))
  observation_counts = dplyr::filter(observation_counts, NumObs > 0)  # should be 1
  observation_counts = dplyr::select(observation_counts, -NumObs)
  input = dplyr::inner_join(input, observation_counts,
                            by = c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                                   "ProductCharge", "IsotopeLabelType"))
  arrow::write_csv_arrow(input, output_path)
  TRUE
}
