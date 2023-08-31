#' @export
#'
BigFragPipetoMSstatsFormat = function(input_file,
                                      max_feature_count=20) {
  input = arrow::open_dataset(input_file, format = "csv")

  input = dplyr::mutate(input, Feature=paste(PeptideSequence, PrecursorCharge,
                                             FragmentIon, ProductCharge, sep="_"))

  feature_counts = dplyr::group_by(input, ProteinName, Feature)
  feature_counts = dplyr::summarize(feature_counts,
                                    MeanAbundance = mean(Intensity, na.rm=TRUE))

  feature_counts = dplyr::collect(feature_counts)

  feature_counts = dplyr::mutate(feature_counts,
                                 feature_rank = dplyr::min_rank(dplyr::desc(MeanAbundance)))

  feature_counts = dplyr::filter(feature_counts,
                                 feature_rank <= max_feature_count)

  feature_counts = dplyr::select(feature_counts, -MeanAbundance, -feature_rank)
  input = dplyr::inner_join(input, feature_counts,
                            by = c("ProteinName", "Feature"))
  input = dplyr::select(input, -Feature)
  input = dplyr::collect(input)


  # input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge,
  #                         FragmentIon, ProductCharge, IsotopeLabelType)
  # observation_counts = dplyr::summarize(input, NumObs = sum(!is.na(Intensity) &
  #                                                             Intensity > 0))
  # observation_counts = dplyr::filter(observation_counts, NumObs > 0)
  # observation_counts = dplyr::select(observation_counts, -NumObs)
  # input = dplyr::inner_join(input, observation_counts,
  #                           by = c("ProteinName", "PeptideSequence",
  #                                  "PrecursorCharge", "FragmentIon",
  #                                  "ProductCharge", "IsotopeLabelType"))
  #
  # pp_df = dplyr::select(input, ProteinName, PeptideSequence)
  # pp_df = dplyr::group_by(pp_df, PeptideSequence)
  # pp_df = dplyr::summarize(pp_df, NumProteins = dplyr::n_distinct(ProteinName))
  # pp_df = dplyr::filter(pp_df, NumProteins == 1)
  # pp_df = dplyr::select(pp_df, -NumProteins)
  # input = dplyr::inner_join(input, pp_df, by = "PeptideSequence")
  #
  # input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge,
  #                         FragmentIon, ProductCharge, IsotopeLabelType, Run,
  #                         Condition, BioReplicate)
  # input = dplyr::summarize(input, Intensity = max(Intensity))
  #
  # input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge,
  #                         FragmentIon,
  #                         ProductCharge, IsotopeLabelType)
  # observation_counts = dplyr::summarize(input, NumObs = sum(!is.na(Intensity) &
  #                                                             Intensity > 1))
  # observation_counts = dplyr::filter(observation_counts, NumObs > 2)
  # observation_counts = dplyr::select(observation_counts, -NumObs)
  # input = dplyr::inner_join(input, observation_counts,
  #                           by = c("ProteinName", "PeptideSequence",
  #                                  "PrecursorCharge", "FragmentIon",
  #                                  "ProductCharge", "IsotopeLabelType"))

  # arrow::write_csv_arrow(input, output_path)
  return(input)
}
