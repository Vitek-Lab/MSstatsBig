#
# library(sparklyr)
# spark_config()
#
# readDFtoSpark = function(spark_connection, spark_master, spark_home,
#                          table_name, input_path, ...) {
#   sc = sparklyr::spark_connect(master = spark_master,
#                                spark_home = spark_home, ...)
#   spec_tbl = spark_read_csv(sc, name = table_name,
#                             path = input_path, header = TRUE, delimiter = ",")
#
# }
#
#
#
#
# sdf_nrow(spec_tbl)
# spec_tbl %>%
#   filter(!is.na(Intensity)) %>%
#   spark_write_table("cleaned2_part2", mode = "overwrite")
#
#
# .filterFewMeasurementsBig = function(input, min_intensity, remove_few) {
#   input = input %>%
#     mutate(is_bigger = as.numeric(Intensity > min_intensity)) %>%
#     group_by(ProteinName, PeptideSequence, PrecursorCharge,
#              FragmentIon, ProductCharge) %>%
#     mutate(n_obs = sum(is_bigger, na.rm = TRUE)) %>%
#     filter(n_obs > 1)
#   if (remove_few) {
#     cutoff = 2
#   } else {
#     cutoff = 0
#   }
#   input = input %>%
#     filter(n_obs > cutoff) %>%
#     select(-n_obs)
#   input
# }
# .handleSingleFeaturePerProteinBig = function(input, remove_single_feature) {
#   if (remove_single_feature) {
#     proteins = input %>%
#       distinct(ProteinName, PeptideSequence, PrecursorCharge,
#                FragmentIon, ProductCharge) %>%
#       group_by(ProteinName) %>%
#       summarize(NumFeatures = n()) %>%
#       filter(NumFeatures > 1) %>%
#       select(ProteinName)
#     input = inner_join(input, proteins, by = "ProteinName")
#   }
#   input
# }
# .cleanByFeatureBig = function(input, remove_few) {
#   input = input %>%
#     group_by(ProteinName, IsotopeLabelType, Run, PeptideSequence,
#              PrecursorCharge, FragmentIon, ProductCharge, Condition, BioReplicate) %>%
#     summarize(Intensity = max(Intensity, na.rm = TRUE))
#   input = .filterFewMeasurementsBig(input, 0, remove_few)
#   input
# }
# .adjustIntensitiesBig = function(input) {
#   # input = mutate(input,
#   #                Intensity = ifelse(is.finite(Intensity), Intensity, NA))
#   input = mutate(input,
#                  Intensity = ifelse(Intensity > 0 & Intensity <= 1, 0, Intensity))
#   input
# }
#
#
# MSstatsPreprocessBig = function(input, remove_few_measurements,
#                                 remove_single_feature,
#                                 output_name) {
#   input = .filterFewMeasurementsBig(input, 1, FALSE)
#   # input = .handleSharedPeptides(input, remove_shared_peptides)
#   input = .cleanByFeatureBig(input, TRUE)
#   input = .handleSingleFeaturePerProteinBig(input, T)
#   input = .adjustIntensitiesBig(input)
#   sparklyr::spark_write_table(input, output_name)
#   input
# }
#
# input = .filterFewMeasurementsBig(cleaned2, 1, FALSE)
# sdf_nrow(input)
# # input = .handleSharedPeptides(input, remove_shared_peptides)
# input = .cleanByFeatureBig(input, TRUE)
# head(input)
# sdf_nrow(input)
# input = .handleSingleFeaturePerProteinBig(input, T)
# sdf_nrow(input2)
# input = .adjustIntensitiesBig(input)
# spark_write_table(input, "cleaned_msstats_spec")
# sdf_nrow(input)
#
# cleaned2 = spark_read_table(sc, "cleaned2")
# cleaned2 %>%
#   filter(IsotopeLabelType == "H") %>%
#   head()
# cleaned2 %>%
#   filter(IsotopeLabelType == "L") %>%
#   head()
#
#
# counts = cleaned2 %>%
#   group_by(ProteinName, PeptideSequence, PrecursorCharge,
#            FragmentIon, ProductCharge, Run) %>%
#   summarize(num_obs = n_distinct(Intensity))
# spark_write_table(counts, "counts_part")
#
# head(counts)
# counts = spark_read_table(sc, "counts_part", mode = "overwrite")
#
# counts %>%
#   filter(num_obs > 1) %>%
#   head(10)
#
#
# cleaned2 %>%
#   filter(PeptideSequence == "LVPFDHAESTYGLYR",
#          PrecursorCharge == 3,
#          FragmentIon == "y7",
#          ProductCharge == 1) %>%
#   filter(Run == "L_D201015_MDIA_P2531_SExp01-BID035_R01") %>%
#   distinct()
# pull(Run)
#
#
# MSstatsPreprocessBig(cleaned2, T, T, "cleaned_part_new")
#
# msstats_input = spark_read_table(sc, "cleaned_part_new")
# sdf_nrow(msstats_input)
# msstats_input = sdf_repartition(msstats_input, partitions = 1)
# spark_write_csv(msstats_input, "/home/mtst/Projekty/BigMSstats/part_new")
#
# MSstatsGetQC = function(input, output_path) {
#   input = sparklyr::sdf_repartition(input, partitions = 1)
#   heavy_peptides = sparklyr::filter(input, IsotopeLabelType == "H")
#   sparklyr::spark_write_csv(heavy_peptides, output_path)
# }
#
# MSstatsGetQC(spec_tbl, "heavy_part")
#
#
# th = fread("./heavy_part/part-00000-c0782809-7300-4423-99dc-676b8c3c8014-c000.csv")
# dim(th)
# head(th)
#
# head(th)
# table(is.na(th$Intensity))
#
#
#
#
# library(sparklyr)
# spark_config()
# sc <- spark_connect(master = "spark://mtst-acer:7077", spark_home = "/opt/spark", version = "3.2")
# spec_tbl <- spark_read_csv(sc, name = "spec_data2", path = "spectronaut_test_heavy.csv", header = TRUE, delimiter = ",")
# sdf_nrow(spec_tbl)
# spec_tbl %>%
#   filter(!is.na(Intensity)) %>%
#   spark_write_table("cleaned2_part2", mode = "overwrite")
#
# spec_tbl %>%
#   filter(!is.na(Intensity)) %>%
#   spark_write_table("cleaned2_test_h", mode = "overwrite")
#
#
# .filterFewMeasurementsBig = function(input, min_intensity, remove_few) {
#   input = sparklyr::mutate(input, is_bigger = as.numeric(Intensity > min_intensity))
#   input = sparklyr::group_by(input, ProteinName, PeptideSequence, PrecursorCharge,
#                              FragmentIon, ProductCharge)
#   input = sparklyr::mutate(input, n_obs = sum(is_bigger, na.rm = TRUE))
#   input = sparklyr::filter(n_obs > 1)
#   if (remove_few) {
#     cutoff = 2
#   } else {
#     cutoff = 0
#   }
#   input = sparklyr::filter(input, n_obs > cutoff)
#   input = sparklyr::select(input, -n_obs)
#   input
# }
# .handleSingleFeaturePerProteinBig = function(input, remove_single_feature) {
#   if (remove_single_feature) {
#     proteins = distinct(input, ProteinName, PeptideSequence, PrecursorCharge,
#                FragmentIon, ProductCharge)
#     proteins = sparklyr::group_by(proteins, ProteinName)
#     proteins = sparklyr::summarize(proteins, NumFeatures = n())
#     proteins = sparklyr::filter(proteins, NumFeatures > 1)
#     proteins = sparklyr::select(proteins, ProteinName)
#     input = sparklyr::inner_join(input, proteins, by = "ProteinName")
#   }
#   input
# }
# .cleanByFeatureBig = function(input, remove_few) {
#   input = sparklyr::group_by(input,
#                    ProteinName, IsotopeLabelType, Run, PeptideSequence,
#                    PrecursorCharge, FragmentIon, ProductCharge, Condition, BioReplicate)
#   input = summarize(input, Intensity = max(Intensity, na.rm = TRUE))
#   input = .filterFewMeasurementsBig(input, 0, remove_few)
#   input
# }
# .adjustIntensitiesBig = function(input) {
#   # input = mutate(input,
#   #                Intensity = ifelse(is.finite(Intensity), Intensity, NA))
#   input = sparklyr::mutate(input,
#                            Intensity = ifelse(Intensity > 0 & Intensity <= 1, 0, Intensity))
#   input
# }
#
#
# MSstatsPreprocessBig = function(input, remove_few_measurements,
#                                 remove_single_feature,
#                                 output_name) {
#   input = .filterFewMeasurementsBig(input, 1, FALSE)
#   # input = .handleSharedPeptides(input, remove_shared_peptides)
#   input = .cleanByFeatureBig(input, TRUE)
#   input = .handleSingleFeaturePerProteinBig(input, T)
#   input = .adjustIntensitiesBig(input)
#   sparklyr::spark_write_table(input, output_name)
#   input
# }
#
# input = .filterFewMeasurementsBig(cleaned2, 1, FALSE)
# sdf_nrow(input)
# # input = .handleSharedPeptides(input, remove_shared_peptides)
# input = .cleanByFeatureBig(input, TRUE)
# head(input)
# sdf_nrow(input)
# input = .handleSingleFeaturePerProteinBig(input, T)
# sdf_nrow(input2)
# input = .adjustIntensitiesBig(input)
# spark_write_table(input, "cleaned_msstats_spec")
# sdf_nrow(input)
#
# cleaned2 %>%
#   filter(PeptideSequence == "LVPFDHAESTYGLYR",
#          PrecursorCharge == 3,
#          FragmentIon == "y7",
#          ProductCharge == 1) %>%
#   filter(Run == "L_D201015_MDIA_P2531_SExp01-BID035_R01") %>%
#   distinct()
# pull(Run)
#
#
# MSstatsPreprocessBig(cleaned2, T, T, "cleaned_part_new")
#
# msstats_input = spark_read_table(sc, "cleaned_part_new")
# sdf_nrow(msstats_input)
# msstats_input = sdf_repartition(msstats_input, partitions = 1)
# spark_write_csv(msstats_input, "/home/mtst/Projekty/BigMSstats/part_new")
#
# MSstatsGetQC = function(input, output_path) {
#   input = sdf_repartition(input, partitions = 1)
#   heavy_peptides = filter(input, IsotopeLabelType == "H")
#   sparklyr::spark_write_csv(heavy_peptides, output_path, mode = "overwrite")
# }
#
# MSstatsGetQC(spec_tbl, "heavy_part")
#
#
# th = fread("./heavy_part/part-00000-c5ed6180-4b98-4f8d-8ff2-5840c6bd5f8d-c000.csv")
# dim(th)
# head(th)
#
# head(th)
# table(is.na(th$Intensity))
#
#
# th_raw = fread("spectronaut_test_heavy.csv")
# th_h = th_raw[IsotopeLabelType == "H"]
# th_h
