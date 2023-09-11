#' @inheritParams MSstatsPreprocessBig
#' @keywords internal
MSstatsPreprocessBigSparklyr = function(connection, input_file, output_file_name,
                                        max_feature_count = 20,
                                        filter_unique_peptides = FALSE,
                                        aggregate_psms = FALSE,
                                        filter_few_obs = FALSE,
                                        remove_annotation = FALSE) {
  tbl_name = getTblName(remove_annotation)

  sparklyr::spark_read_csv(connection, "mstinput", path = input_file,
                           header = TRUE, delimiter = ",",
                           memory = FALSE, overwrite = TRUE,
                           repartition = 100)

  if (remove_annotation) {
    dplyr::tbl(sc, "mstinput") %>%
      dplyr::select(-Condition, -BioReplicate) %>%
      sparklyr::spark_write_table("mstinputnoannot", mode = "overwrite")
    DBI::dbRemoveTable(connection, "mstinput")
  }

  dplyr::tbl(sc, "mstinputnoannot") %>%
    dplyr::group_by(PeptideSequence) %>%
    dplyr::summarize(NumProteins = n_distinct(ProteinName)) %>%
    dplyr::filter(NumProteins > 1) %>%
    sparklyr::spark_write_table("sharedpeptides", mode = "overwrite")
  sparklyr::sdf_repartition(tbl(connection, "mstinputnoannot"), partition_by = "ProteinName")

  dplyr::tbl(sc, "mstinputnoannot") %>%
    dplyr::group_by(ProteinName, PeptideSequence,
                    PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType,
                    Run) %>%
    dplyr::summarize(Intensity = max(Intensity)) %>%
    sparklyr::spark_write_table("mstinputagg", mode = "overwrite")
  DBI::dbRemoveTable(connection, "mstinputnoannot")

  dplyr::tbl(sc, "mstinputagg") %>%
    dplyr::group_by(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType) %>%
    dplyr::summarize(NumObs = sum(as.numeric(!is.na(Intensity) & Intensity > 1))) %>%
    dplyr::filter(NumObs <= 2) %>%
    dplyr::select(-NumObs) %>%
    sparklyr::spark_write_table(name = "less2obs", mode = "overwrite")

  dplyr::anti_join(
    dplyr::tbl(sc, "mstinputagg"),
    tbl(sc, "less2obs"),
    by = c("ProteinName",
           "PeptideSequence", "PrecursorCharge", "FragmentIon",
           "ProductCharge", "IsotopeLabelType")
  ) %>%
    sparklyr::spark_write_table(name = "mstinputfilt", mode = "overwrite")
  DBI::dbRemoveTable(connection, "mstinputagg")
  DBI::dbRemoveTable(connection, "less2obs")

  dplyr::tbl(connection, "mstinputfilt") %>%
    dplyr::group_by(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType) %>%
    dplyr::summarize(Intensity = mean(Intensity)) %>%
    dplyr::ungroup() %>%
    sparklyr::spark_write_table("fortop50", mode = "overwrite", partition_by = "ProteinName")


  dplyr::tbl(connection, "fortop50") %>%
    dplyr::group_by(ProteinName, IsotopeLabelType) %>%
    dplyr::mutate(Rank = min_rank(-Intensity)) %>%
    dplyr::filter(Rank <= max_feature_count) %>%
    dplyr::select(-Rank, -Intensity) %>%
    sparklyr::spark_write_table("top50", mode = "overwrite", partition_by = "ProteinName")

  dplyr::inner_join(
    dplyr::tbl(sc, "mstinputfilt"),
    dplyr::tbl(sc, "top50"),
    by = c("ProteinName", "IsotopeLabelType", "PeptideSequence", "PrecursorCharge",
           "FragmentIon", "ProductCharge")) %>%
    sparklyr::spark_write_table("fgprocessed", mode = "overwrite")
  DBI::dbRemoveTable(connection, "top50")
  DBI::dbRemoveTable(connection, "fortop50")
  DBI::dbRemoveTable(connection, "mstinputfilt")

  sparklyr::spark_write_csv(sparklyr::sdf_repartition(tbl(sc, "fgprocessed"), 1),
                            output_file_name, mode = "overwrite")
  return(TRUE)

}


#' @inheritParams MSstatsPreprocessBig
#' @keywords internal
MSstatsPreprocessBigArrow = function(input_file,
                                     output_file_name,
                                     max_feature_count = 20,
                                     filter_unique_peptides = FALSE,
                                     aggregate_psms = FALSE,
                                     filter_few_obs = FALSE) {
  input = arrow::open_dataset(input_file, format = "csv")

  input = dplyr::mutate(input,
                        Feature = paste(PeptideSequence, PrecursorCharge,
                                        FragmentIon, ProductCharge, sep = "_"))
  feature_counts = dplyr::group_by(input, ProteinName, Feature)
  feature_counts = dplyr::summarize(feature_counts,
                                    MeanAbundance = mean(Intensity,
                                                         na.rm = TRUE))
  feature_counts = dplyr::collect(feature_counts)

  feature_counts = dplyr::mutate(
    feature_counts,
    feature_rank = dplyr::min_rank(dplyr::desc(MeanAbundance)))

  feature_counts = dplyr::filter(feature_counts,
                                 feature_rank <= max_feature_count)

  feature_counts = dplyr::select(feature_counts, -MeanAbundance, -feature_rank)
  input = dplyr::inner_join(input, feature_counts,
                            by = c("ProteinName", "Feature"))
  input = dplyr::select(input, -Feature)
  input = dplyr::collect(input)

  arrow::write_csv_arrow(input, file = paste0("topN_", output_file_name))

  if (filter_unique_peptides) {
    pp_df = dplyr::select(input, ProteinName, PeptideSequence)
    pp_df = dplyr::group_by(pp_df, PeptideSequence)
    pp_df = dplyr::summarize(pp_df, NumProteins = dplyr::n_distinct(ProteinName))
    pp_df = dplyr::filter(pp_df, NumProteins > 1)
    pp_df = dplyr::select(pp_df, -NumProteins)
    input = dplyr::anti_join(input, pp_df, by = "PeptideSequence")
  }

  if (aggregate_psms) {
    input = dplyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge,
                            FragmentIon, ProductCharge, IsotopeLabelType, Run,
                            Condition, BioReplicate)
    input = dplyr::summarize(input, Intensity = max(Intensity, na.rm = TRUE))
  }

  if (filter_few_obs) {
    input = dplyr::group_by(input, ProteinName, Feature)
    observation_counts = dplyr::summarize(input,
                                          NumObs = sum(!is.na(Intensity) &
                                                         Intensity > 0))
    observation_counts = dplyr::filter(observation_counts, NumObs <= 2)
    observation_counts = dplyr::select(observation_counts, -NumObs)
    input = dplyr::anti_join(input, observation_counts,
                             by = c("ProteinName", "Feature"))
  }

  arrow::write_csv_arrow(input, file = output_file_name)
  input
}

getTblName = function(remove_annotation) {
  if (remove_annotation) {
    "mstinput"
  } else {
    "mstinputnoannot"
  }
}
