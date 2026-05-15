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

<<<<<<< HEAD
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
=======
make_spectronaut_input <- function(n = 1L, ...) {
  base <- data.frame(
    R.FileName = rep("run1", n),
    R.Condition = "A",
    R.Replicate = 1L,
    PG.ProteinAccessions = "P1",
    EG.ModifiedSequence = paste0("PEP", seq_len(n)),
    FG.LabeledSequence = paste0("PEP", seq_len(n)),
    FG.Charge = 2L,
    F.FrgIon = "y1",
    F.Charge = 1L,
    EG.Identified = "True",
    F.ExcludedFromQuantification = "False",
    F.FrgLossType = "noloss",
    PG.Qvalue = 0.001,
    EG.Qvalue = 0.001,
    F.NormalizedPeakArea = seq_len(n) * 100,
    stringsAsFactors = FALSE
  )
  overrides <- list(...)
  for (col in names(overrides)) base[[col]] <- overrides[[col]]
  base
}

test_that("cleanSpectronautChunk produces the expected MSstats schema", {
  input <- make_spectronaut_input(n = 1L)
  output_file <- tempfile(fileext = ".csv")
  on.exit(unlink(output_file, force = TRUE), add = TRUE)

  MSstatsBig:::cleanSpectronautChunk(input, output_file, pos = 1L)
  result <- readr::read_csv(output_file, show_col_types = FALSE)

  expected_cols <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
                     "ProductCharge", "IsotopeLabelType", "Run", "BioReplicate", "Condition",
                     "Intensity")
  expect_setequal(colnames(result), expected_cols)
  expect_equal(nrow(result), 1L)
  expect_equal(result$Intensity, 100)
  expect_equal(result$IsotopeLabelType, "L")
})

test_that("cleanSpectronautChunk filter_by_excluded sets Intensity to NA on excluded rows", {
  input <- make_spectronaut_input(
    n = 2L,
    F.ExcludedFromQuantification = c("True", "False"),
    F.NormalizedPeakArea = c(100, 200)
  )
  output_file <- tempfile(fileext = ".csv")
  on.exit(unlink(output_file, force = TRUE), add = TRUE)

  MSstatsBig:::cleanSpectronautChunk(input, output_file, pos = 1L,
                                     filter_by_excluded = TRUE,
                                     filter_by_qvalue = FALSE)
  result <- readr::read_csv(output_file, show_col_types = FALSE)
  result <- result[order(result$PeptideSequence), ]

  expect_true(is.na(result$Intensity[1]))
  expect_equal(result$Intensity[2], 200)
})

test_that("cleanSpectronautChunk filter_by_identified sets Intensity to NA on unidentified rows", {
  input <- make_spectronaut_input(
    n = 2L,
    EG.Identified = c("True", "False"),
    F.NormalizedPeakArea = c(100, 200)
  )
  output_file <- tempfile(fileext = ".csv")
  on.exit(unlink(output_file, force = TRUE), add = TRUE)

  MSstatsBig:::cleanSpectronautChunk(input, output_file, pos = 1L,
                                     filter_by_identified = TRUE,
                                     filter_by_qvalue = FALSE)
  result <- readr::read_csv(output_file, show_col_types = FALSE)
  result <- result[order(result$PeptideSequence), ]

  expect_equal(result$Intensity[1], 100)
  expect_true(is.na(result$Intensity[2]))
})

test_that("cleanSpectronautChunk filter_by_qvalue NA-aware semantics match dplyr::if_else", {
  input <- make_spectronaut_input(
    n = 3L,
    EG.Qvalue = c(0.001, 0.5, NA_real_),
    F.NormalizedPeakArea = c(100, 200, 300)
  )
  output_file <- tempfile(fileext = ".csv")
  on.exit(unlink(output_file, force = TRUE), add = TRUE)

  MSstatsBig:::cleanSpectronautChunk(input, output_file, pos = 1L,
                                     filter_by_qvalue = TRUE,
                                     qvalue_cutoff = 0.01)
  result <- readr::read_csv(output_file, show_col_types = FALSE)
  result <- result[order(result$PeptideSequence), ]

  # PEP1: EGQvalue below cutoff -> kept
  expect_equal(result$Intensity[1], 100)
  # PEP2: EGQvalue above cutoff -> NA
  expect_true(is.na(result$Intensity[2]))
  # PEP3: EGQvalue NA -> NA (preserves dplyr::if_else semantics)
  expect_true(is.na(result$Intensity[3]))
})

test_that("cleanSpectronautChunk drops rows where FFrgLossType != noloss", {
  input <- make_spectronaut_input(
    n = 3L,
    F.FrgLossType = c("noloss", "H2O", "noloss")
  )
  output_file <- tempfile(fileext = ".csv")
  on.exit(unlink(output_file, force = TRUE), add = TRUE)

  MSstatsBig:::cleanSpectronautChunk(input, output_file, pos = 1L,
                                     filter_by_qvalue = FALSE)
  result <- readr::read_csv(output_file, show_col_types = FALSE)

  expect_equal(nrow(result), 2L)
  expect_setequal(result$PeptideSequence, c("PEP1", "PEP3"))
>>>>>>> fdd7476 (add more tests)
})

test_that("reduceBigSpectronaut rejects invalid block_size values", {
  input_file <- tempfile(fileext = ".csv")
  writeLines("a,b\n1,2", input_file)
  output_file <- tempfile()
  on.exit({
    unlink(input_file, force = TRUE)
    unlink(output_file, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  expect_error(reduceBigSpectronaut(input_file, output_file, block_size = -1L))
  expect_error(reduceBigSpectronaut(input_file, output_file, block_size = 0L))
  expect_error(reduceBigSpectronaut(input_file, output_file, block_size = NA_integer_))
  expect_error(reduceBigSpectronaut(input_file, output_file, block_size = c(1L, 2L)))
  expect_error(suppressWarnings(
    reduceBigSpectronaut(input_file, output_file, block_size = "16MB")
  ))
})

test_that("bigSpectronauttoMSstatsFormat plumbs block_size through to reduceBigSpectronaut", {
  captured <- new.env(parent = emptyenv())
  captured$block_size <- NULL

  spy_reduce <- function(input_file, output_path, intensity, filter_by_excluded,
                        filter_by_identified, filter_by_qvalue, qvalue_cutoff,
                        calculateAnomalyScores, anomalyModelFeatures,
                        block_size = 16L * 1024L * 1024L) {
    captured$block_size <- block_size
    msstats_data <- data.frame(
      ProteinName = "P1", PeptideSequence = "PEPTIDE", PrecursorCharge = 2,
      FragmentIon = "frag1", ProductCharge = 1,
      IsotopeLabelType = "L", Condition = "A", BioReplicate = 1,
      Run = "run1", Intensity = 100
    )
    readr::write_csv(msstats_data, output_path)
  }

  input_file <- "dummy_spectro_input.csv"

  # Default forwards 16 MiB.
  stub(bigSpectronauttoMSstatsFormat, "reduceBigSpectronaut", spy_reduce)
  output_file_default <- tempfile(fileext = ".csv")
  on.exit({
    unlink(output_file_default, recursive = TRUE, force = TRUE)
    unlink(paste0("reduce_output_", basename(output_file_default)),
           recursive = TRUE, force = TRUE)
  }, add = TRUE)
  bigSpectronauttoMSstatsFormat(
    input_file = input_file,
    output_file_name = output_file_default,
    backend = "arrow"
  )
  expect_identical(captured$block_size, 16L * 1024L * 1024L)

  # Override forwards the user's value.
  output_file_override <- tempfile(fileext = ".csv")
  on.exit({
    unlink(output_file_override, recursive = TRUE, force = TRUE)
    unlink(paste0("reduce_output_", basename(output_file_override)),
           recursive = TRUE, force = TRUE)
  }, add = TRUE)
  bigSpectronauttoMSstatsFormat(
    input_file = input_file,
    output_file_name = output_file_override,
    backend = "arrow",
    block_size = 8L * 1024L * 1024L
  )
  expect_identical(captured$block_size, 8L * 1024L * 1024L)
})

# test_that("bigDIANNtoMSstatsFormat works with real MSstatsConvert tinytest data", {
#   input_file <- "/Users/rudhikshah/NorthEasternContractWork/MSstatsConvert/inst/tinytest/raw_data/DIANN/diann_input.tsv"
#   annotation_file <- "/Users/rudhikshah/NorthEasternContractWork/MSstatsConvert/inst/tinytest/raw_data/DIANN/annotation.csv"

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