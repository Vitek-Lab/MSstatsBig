library(testthat)
library(mockery)

context("DIANN cleaning")

test_that("cleanDIANNChunk passes annotation to MSstatsMakeAnnotation", {
  # Prepare data
  input_chunk <- data.frame(Run = "Run1", Intensity = 100)
  annotation <- data.frame(Run = "Run1", Condition = "A", BioReplicate = 1)
  
  # Mocks
  m_import <- mock(input_chunk)
  m_clean <- mock(input_chunk)
  m_annotate <- mock(merge(input_chunk, annotation, by = "Run"))
  m_write <- mock(NULL)
  
  stub(cleanDIANNChunk, "MSstatsImport", m_import)
  stub(cleanDIANNChunk, "MSstatsClean", m_clean)
  stub(cleanDIANNChunk, "MSstatsMakeAnnotation", m_annotate)
  stub(cleanDIANNChunk, ".writeChunkToFile", m_write)
  
  # Execute
  cleanDIANNChunk(input_chunk, "output.csv", MBR = TRUE, 
                  quantificationColumn = "Intensity", pos = 1, 
                  annotation = annotation)
  
  # Verify
  expect_called(m_annotate, 1)
  
  # Check arguments passed to MSstatsMakeAnnotation
  args <- mock_args(m_annotate)[[1]]
  expect_equal(args[[1]], input_chunk) # Input from Clean
  expect_equal(args[[2]], annotation)  # Annotation passed through
  
  # Check that the result of annotation is passed to write
  expect_called(m_write, 1)
  write_args <- mock_args(m_write)[[1]]
  expect_equal(write_args[[1]], merge(input_chunk, annotation, by = "Run"))
})