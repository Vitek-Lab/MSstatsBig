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
    Lib.Q.Value = 0.001,
    Lib.PG.Q.Value = 0.001,
    stringsAsFactors = FALSE
  )

  # The function is not exported, so we use :::
  MSstatsBig:::cleanDIANNChunk(diann_chunk_data, output_file, MBR = TRUE,
                               quantificationColumn = "FragmentQuantCorrected", pos = 1)

  result <- read.csv(output_file)

  expect_equal(nrow(result), 1)
  expect_equal(result$ProteinName, "ProteinA")
  expect_equal(result$PeptideSequence, "PEPTIDE(mod)")
  expect_equal(result$Intensity, 100)
  expect_equal(result$FragmentIon, "y7^1/1")
  expect_true("PeptideSequence" %in% colnames(result))

  file.remove(output_file)
})

test_that("cleanDIANNChunk handles 'auto' quantification column correctly", {
  output_file <- tempfile(fileext = ".csv")

  # Data with wide format fragment quantification
  diann_chunk_wide <- data.frame(
    Run = "run1",
    Protein.Names = "ProteinA",
    Stripped.Sequence = "PEPTIDE",
    Modified.Sequence = "PEPTIDE",
    Precursor.Charge = 2,
    Fr1Quantity = 100,
    Fr2Quantity = 200,
    Q.Value = 0.005,
    Precursor.Mz = 400.5,
    Fragment.Info = "y1^1/1;y2^1/1",
    Lib.Q.Value = 0.001,
    Lib.PG.Q.Value = 0.001,
    stringsAsFactors = FALSE
  )

  MSstatsBig:::cleanDIANNChunk(diann_chunk_wide, output_file, MBR = TRUE,
                               quantificationColumn = "auto", pos = 1)

  result <- read.csv(output_file)

  expect_equal(nrow(result), 2)
  expect_equal(sort(result$Intensity), c(100, 200))
  expect_equal(sort(result$FragmentIon), c("y1^1/1", "y2^1/1"))

  file.remove(output_file)

  # Test error when columns are missing
  diann_chunk_missing <- diann_chunk_wide[, !grepl("Quantity", names(diann_chunk_wide))]
  expect_error(MSstatsBig:::cleanDIANNChunk(diann_chunk_missing, output_file, MBR = TRUE,
                                            quantificationColumn = "auto", pos = 1),
               "No fragment quantification columns found")
})

test_that("cleanDIANNChunk handles missing Fragment.Info by defaulting ProductCharge to 1", {
  output_file <- tempfile(fileext = ".csv")

  # Data with missing Fragment.Info (simulating it not being present)
  diann_chunk_missing <- data.frame(
    Run = "run1",
    Protein.Names = "ProteinA",
    Stripped.Sequence = "PEPTIDE",
    Modified.Sequence = "PEPTIDE",
    Precursor.Charge = 2,
    Fragment.Quant.Corrected = 100,
    Q.Value = 0.005,
    Precursor.Mz = 400.5,
    # Fragment.Info is missing
    Lib.Q.Value = 0.001,
    Lib.PG.Q.Value = 0.001,
    stringsAsFactors = FALSE
  )

  MSstatsBig:::cleanDIANNChunk(diann_chunk_missing, output_file, MBR = TRUE,
                               quantificationColumn = "FragmentQuantCorrected", pos = 1)

  result <- read.csv(output_file)

  expect_equal(nrow(result), 1)
  expect_equal(result$ProductCharge, 1)
  expect_equal(result$FragmentIon, "Frag1")

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
    Lib.Q.Value = c(0.001, 0.002),
    Lib.PG.Q.Value = c(0.001, 0.002),
    stringsAsFactors = FALSE
  )
  write.csv(diann_data, input_file, row.names = FALSE)

  MSstatsBig:::reduceBigDIANN(input_file, output_file, MBR = TRUE,
                              quantificationColumn = "FragmentQuantCorrected")

  result <- read.csv(output_file)
  expect_equal(nrow(result), 2)
  expect_equal(result$Intensity, c(100, 300))
  expect_equal(result$ProteinName, c("ProteinA", "ProteinB"))
  expect_equal(result$FragmentIon, c("y7^1/1", "y5^1/2"))

  file.remove(input_file)
  file.remove(output_file)
})

# End-to-end test for bigDIANNtoMSstatsFormat
test_that("bigDIANNtoMSstatsFormat works with arrow backend", {
  input_file <- tempfile(fileext = ".csv")
  output_file <- basename(tempfile(fileext = ".csv"))

  # 4 features for one protein. Feature selection should pick the top 2.
  diann_data <- rbind(
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(1000, 1100), Q.Value = 0.001, Precursor.Mz = 500, Fragment.Info = "y1", Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001),
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(500, 600), Q.Value = 0.001, Precursor.Mz = 500, Fragment.Info = "y2", Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001),
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(100, 100), Q.Value = 0.001, Precursor.Mz = 500, Fragment.Info = "y3", Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001),
    data.frame(Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, Fragment.Quant.Corrected = c(2000, 2100), Q.Value = 0.001, Precursor.Mz = 500, Fragment.Info = "y4", Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001)
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

  # Cleanup — outputs may be directories when backend = "arrow"
  unlink(input_file, force = TRUE)
  unlink(output_file, recursive = TRUE, force = TRUE)
  unlink(paste0("reduce_output_", output_file), recursive = TRUE, force = TRUE)
  unlink(paste0("topN_", output_file), recursive = TRUE, force = TRUE)
})

test_that("bigDIANNtoMSstatsFormat works with annotation", {
  input_file <- tempfile(fileext = ".csv")
  output_file <- basename(tempfile(fileext = ".csv"))

  # Minimal DIANN data
  diann_data <- data.frame(
    Run = c("r1", "r2"), Protein.Names = "P1", Stripped.Sequence = "PEPTIDE", 
    Modified.Sequence = "PEPTIDE", Precursor.Charge = 2, 
    Fragment.Quant.Corrected = c(1000, 1100), Q.Value = 0.001, Precursor.Mz = 500, 
    Fragment.Info = "y1", Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001
  )
  write.csv(diann_data, input_file, row.names = FALSE)
  
  # Annotation data
  annot <- data.frame(
    Run = c("r1", "r2"),
    Condition = c("Disease", "Healthy"),
    BioReplicate = c(1, 2)
  )
  
  converted <- bigDIANNtoMSstatsFormat(
    input_file = input_file,
    annotation = annot,
    output_file_name = output_file,
    backend = "arrow"
  )
  result <- dplyr::collect(converted)

  expect_true(all(c("Condition", "BioReplicate") %in% colnames(result)))
  expect_equal(result$Condition[result$Run == "r1"], "Disease")
  expect_equal(result$Condition[result$Run == "r2"], "Healthy")

  # Cleanup — outputs may be directories when backend = "arrow"
  unlink(input_file, force = TRUE)
  unlink(output_file, recursive = TRUE, force = TRUE)
  unlink(paste0("reduce_output_", output_file), recursive = TRUE, force = TRUE)
  unlink(paste0("topN_", output_file), recursive = TRUE, force = TRUE)
  unlink(paste0("cleaned_", output_file), recursive = TRUE, force = TRUE)
})