library(testthat)
library(mockery)

context("General converter functions")

test_that("MSstatsAddAnnotationBig adds annotation correctly", {
  input_data <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    Intensity = c(100, 200, 300)
  )

  annotation_data <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    Condition = c("A", "A", "B"),
    BioReplicate = c(1, 2, 1)
  )

  expected_output <- data.frame(
    Run = c("Run1", "Run2", "Run3"),
    Intensity = c(100, 200, 300),
    Condition = c("A", "A", "B"),
    BioReplicate = c(1, 2, 1)
  )

  result <- MSstatsAddAnnotationBig(input_data, annotation_data)

  expect_equal(result, expected_output)
})

test_that("MSstatsPreprocessBig performs feature selection correctly", {
  input_file <- tempfile(fileext = ".csv")
  output_file <- "preprocess_output.csv"

  # P1 has 3 features (frag1, frag2, frag3). frag3 has the highest avg intensity.
  # P2 has 2 features (fragA, fragB). fragB has the highest avg intensity.
  msstats_data <- rbind(
    data.frame(ProteinName = "P1", PeptideSequence = "PEPTIDE", PrecursorCharge = 2, FragmentIon = rep(c("frag1", "frag2", "frag3"), each = 2), ProductCharge = 1, IsotopeLabelType = "L", Condition = "A", BioReplicate = rep(1:2, 3), Run = rep(c("run1", "run2"), 3), Intensity = c(1000, 1100, 500, 550, 2000, 2100)),
    data.frame(ProteinName = "P2", PeptideSequence = "PEPTIDE2", PrecursorCharge = 3, FragmentIon = rep(c("fragA", "fragB"), each = 2), ProductCharge = 1, IsotopeLabelType = "L", Condition = "B", BioReplicate = rep(1:2, 2), Run = rep(c("run1", "run2"), 2), Intensity = c(100, 150, 800, 850))
  )
  readr::write_csv(msstats_data, input_file)

  processed <- MSstatsPreprocessBig(input_file, output_file, backend = "arrow",
                                    max_feature_count = 1)
  result <- dplyr::collect(processed)

  # For P1, frag3 should be selected. For P2, fragB should be selected.
  expect_equal(nrow(result), 4)

  p1_result <- result[result$ProteinName == "P1", ]
  expect_equal(nrow(p1_result), 2)
  expect_true(all(p1_result$FragmentIon == "frag3"))

  p2_result <- result[result$ProteinName == "P2", ]
  expect_equal(nrow(p2_result), 2)
  expect_true(all(p2_result$FragmentIon == "fragB"))

  # Cleanup — output may be a directory when backend = "arrow"
  unlink(input_file, force = TRUE)
  unlink(output_file, recursive = TRUE, force = TRUE)
})

test_that("bigSpectronauttoMSstatsFormat works correctly", {
  # Mock reduceBigSpectronaut as its source is not provided
  mock_reduce <- mock(NULL)

  stub(bigSpectronauttoMSstatsFormat, "reduceBigSpectronaut", function(input_file, output_path, ...) {
    msstats_data <- data.frame(
      ProteinName = "P1", PeptideSequence = "PEPTIDE", PrecursorCharge = 2,
      FragmentIon = rep(c("frag1", "frag2"), each = 2), ProductCharge = 1,
      IsotopeLabelType = "L", Condition = "A", BioReplicate = rep(1:2, 2),
      Run = rep(c("run1", "run2"), 2), Intensity = c(1000, 1100, 2000, 2100) # frag2 is higher
    )
    readr::write_csv(msstats_data, output_path)
  })

  input_file <- "dummy_spectro_input.csv"
  output_file <- "spectro_output.csv"

  processed <- bigSpectronauttoMSstatsFormat(
    input_file = input_file,
    output_file_name = output_file,
    backend = "arrow",
    max_feature_count = 1
  )
  result <- dplyr::collect(processed)

  # The mock reduce function creates a file with 2 features for P1.
  # max_feature_count = 1 should select frag2.
  expect_equal(nrow(result), 2)
  expect_true(all(result$FragmentIon == "frag2"))

  # Cleanup — outputs may be directories when backend = "arrow"
  unlink(output_file, recursive = TRUE, force = TRUE)
  unlink(paste0("reduce_output_", output_file), recursive = TRUE, force = TRUE)
})

test_that("bigSpectronauttoMSstatsFormat overrides Condition/BioReplicate from annotation", {
  stub(bigSpectronauttoMSstatsFormat, "reduceBigSpectronaut", function(input_file, output_path, ...) {
    msstats_data <- data.frame(
      ProteinName = "P1", PeptideSequence = "PEPTIDE", PrecursorCharge = 2,
      FragmentIon = "frag1", ProductCharge = 1,
      IsotopeLabelType = "L",
      Condition = "FROM_SPECTRONAUT", BioReplicate = 999,
      Run = rep(c("run1", "run2"), each = 1),
      Intensity = c(1000, 2000)
    )
    readr::write_csv(msstats_data, output_path)
  })

  input_file <- "dummy_spectro_input.csv"
  output_file <- "spectro_output_annot.csv"

  annotation <- data.frame(
    Run = c("run1", "run2"),
    BioReplicate = c(7L, 8L),
    Condition = c("ctrl", "treat"),
    stringsAsFactors = FALSE
  )

  processed <- bigSpectronauttoMSstatsFormat(
    input_file = input_file,
    annotation = annotation,
    output_file_name = output_file,
    backend = "arrow",
    max_feature_count = 1
  )
  result <- dplyr::collect(processed)
  result <- result[order(result$Run), ]

  expect_equal(result$Condition, c("ctrl", "treat"))
  expect_equal(result$BioReplicate, c(7L, 8L))
  expect_false(any(result$Condition == "FROM_SPECTRONAUT"))
  expect_false(any(result$BioReplicate == 999))

  # Cleanup
  unlink(output_file, recursive = TRUE, force = TRUE)
  unlink(paste0("reduce_output_", output_file), recursive = TRUE, force = TRUE)
})