library(testthat)

context("DIANN converter functions")

# Test for the internal cleanDIANNChunk function
test_that("cleanDIANNChunk processes data correctly", {
  output_file <- tempfile(fileext = ".csv")

  diann_chunk_data <- data.frame(
    Run = "run1",
    Protein.Names = "ProteinA",
    Stripped.Sequence = "PEPTIDE",
    Modified.Sequence = "PEPTIDE(mod)",
    Precursor.Charge = 2,
    Fragment.Quant.Corrected = "100;200",
    Q.Value = 0.005,
    Precursor.Mz = 400.5,
    Fragment.Info = "y7^1/1;b3-H2O^1/1", # One valid, one to be filtered
    Lib.Q.Value = 0.01,
    Lib.PG.Q.Value = 0.001,
    stringsAsFactors = FALSE
  )

  # The function is not exported, so we use :::
  MSstatsBig:::cleanDIANNChunk(diann_chunk_data, output_file, MBR = TRUE,
                               quantificationColumn = "Fragment.Quant.Corrected", pos = 1)

  result <- read.csv(output_file)

  expect_equal(nrow(result), 1)
  expect_equal(result$ProteinName, "ProteinA")
  expect_equal(result$PeptideSequence, "PEPTIDE")
  expect_equal(result$Intensity, 100)
  expect_equal(result$FragmentIon, "y7^1/1")
  expect_equal(result$ProductCharge, 1)
  expect_equal(result$IsotopeLabelType, "L")
  expect_true("PeptideModifiedSequence" %in% colnames(result))

  file.remove(output_file)
})

# Test for the internal reduceBigDIANN function
test_that("reduceBigDIANN processes a file correctly", {
  input_file <- tempfile(fileext = ".csv")
  output_file <- tempfile(fileext = ".csv")

  diann_data <- data.frame(
    Run = c("run1", "run1"),
    Protein.Names = c("ProteinA", "ProteinB"),
    Stripped.Sequence = c("PEPTIDE_A", "PEPTIDE_B"),
    Modified.Sequence = c("PEPTIDE_A(mod)", "PEPTIDE_B"),
    Precursor.Charge = c(2, 3),
    Fragment.Quant.Corrected = c("100;200", "300"),
    Q.Value = c(0.005, 0.006),
    Precursor.Mz = c(400.5, 500.5),
    Fragment.Info = c("y7^1/1;b3-H2O^1/1", "y5^1/2"),
    Lib.Q.Value = c(0.01, 0.02),
    Lib.PG.Q.Value = c(0.001, 0.002),
    stringsAsFactors = FALSE
  )
  write.csv(diann_data, input_file, row.names = FALSE)

  MSstatsBig:::reduceBigDIANN(input_file, output_file, MBR = TRUE,
                              quantificationColumn = "Fragment.Quant.Corrected")

  result <- read.csv(output_file)
  expect_equal(nrow(result), 2)
  expect_equal(result$Intensity, c(100, 300))
  expect_equal(result$ProteinName, c("ProteinA", "ProteinB"))
  expect_equal(result$ProductCharge, c(1, 2))
  expect_equal(result$FragmentIon, c("y7^1/1", "y5^1/2"))

  file.remove(input_file)
  file.remove(output_file)
})

# End-to-end test for bigDIANNtoMSstatsFormat
test_that("bigDIANNtoMSstatsFormat works with arrow backend", {
  input_file <- tempfile(fileext = ".csv")
  output_file <- "test_diann_output.csv"

  # 4 features for one protein. Feature selection should pick the top 2.
  diann_data <- rbind(
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(1000, 1100), Q.Value = 0.01, Precursor.Mz = 500, Fragment.Info = "y1", Lib.Q.Value = 0.01, Lib.PG.Q.Value = 0.01),
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(500, 600), Q.Value = 0.01, Precursor.Mz = 500, Fragment.Info = "y2", Lib.Q.Value = 0.01, Lib.PG.Q.Value = 0.01),
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(100, 100), Q.Value = 0.01, Precursor.Mz = 500, Fragment.Info = "y3", Lib.Q.Value = 0.01, Lib.PG.Q.Value = 0.01),
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(2000, 2100), Q.Value = 0.01, Precursor.Mz = 500, Fragment.Info = "y4", Lib.Q.Value = 0.01, Lib.PG.Q.Value = 0.01)
  )
  write.csv(diann_data, input_file, row.names = FALSE)

  converted <- bigDIANNtoMSstatsFormat(
    input_file = input_file,
    output_file_name = output_file,
    backend = "arrow",
    max_feature_count = 2
  )
  result <- dplyr::collect(converted)

  # Avg intensities: y1=1050, y2=550, y3=100, y4=2050.
  # Top 2 features are y4 and y1.
  expect_equal(nrow(result), 4) # 2 features * 2 runs
  expect_true(all(c("y1", "y4") %in% unique(result$FragmentIon)))
  expect_false(any(c("y2", "y3") %in% unique(result$FragmentIon)))

  # Cleanup
  file.remove(input_file)
  if (file.exists(output_file)) file.remove(output_file)
  if (file.exists(paste0("reduce_output_", output_file))) file.remove(paste0("reduce_output_", output_file))
})