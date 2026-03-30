#' Read and clean a large DIANN file in chunks
#' 
#' @param input_file Path to the input DIANN file
#' @param output_path Path to the output CSV file
#' @param MBR Boolean, whether MBR was used
#' @param quantificationColumn Name of the column containing intensity values
#' @param global_qvalue_cutoff Global Q-value cutoff
#' @param qvalue_cutoff Q-value cutoff
#' @param pg_qvalue_cutoff Protein group Q-value cutoff
#' @param calculateAnomalyScores Boolean for MSstats+ Model
#' @param anomalyModelFeatures Character vector of features to use for MSstats+ Model
#' @param annotation Annotation file or data frame
#' @return NULL. Writes to file.
#' @keywords internal
reduceBigDIANN <- function(input_file, output_path, MBR = TRUE,
                           quantificationColumn = "FragmentQuantCorrected",
                           global_qvalue_cutoff = 0.01,
                           qvalue_cutoff = 0.01,
                           pg_qvalue_cutoff = 0.01,
                           calculateAnomalyScores=FALSE, 
                           anomalyModelFeatures=c(),
                           annotation = NULL) {
  first_line <- readLines(input_file, n = 1)
  if (grepl("\t", first_line)) {
    delim <- "\t"
  } else if (grepl(",", first_line)) {
    delim <- ","
  } else if (grepl(";", first_line)) {
    delim <- ";"
  } else {
    delim <- "\t"
  }
  
  diann_chunk <- function(x, pos) cleanDIANNChunk(x, output_path, MBR, 
                                                  quantificationColumn, pos,
                                                  global_qvalue_cutoff, 
                                                  qvalue_cutoff, 
                                                  pg_qvalue_cutoff, 
                                                  calculateAnomalyScores,
                                                  anomalyModelFeatures,
                                                  annotation)

  readr::read_delim_chunked(input_file,
                            readr::DataFrameCallback$new(diann_chunk),
                            delim = delim,
                            chunk_size = 1e6)
}

#' Clean a single chunk of DIANN data
#' 
#' @param input Data frame chunk
#' @param output_path Path to output file
#' @param MBR Boolean, whether MBR was used
#' @param quantificationColumn Name of intensity column
#' @param pos Chunk position (1 for first chunk, >1 for subsequent)
#' @param global_qvalue_cutoff Global Q-value cutoff
#' @param qvalue_cutoff Q-value cutoff
#' @param pg_qvalue_cutoff Protein group Q-value cutoff
#' @param calculateAnomalyScores Boolean for MSstats+ Model
#' @param anomalyModelFeatures Character vector of features to use for MSstats+ Model
#' @param annotation Annotation file or data frame
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean MSstatsMakeAnnotation
#' @return NULL
#' @keywords internal
cleanDIANNChunk = function(input, output_path, MBR, quantificationColumn, pos,
                           global_qvalue_cutoff = 0.01,
                           qvalue_cutoff = 0.01,
                           pg_qvalue_cutoff = 0.01,
                           calculateAnomalyScores=FALSE,
                           anomalyModelFeatures = c(),
                           annotation = NULL) {
    input = MSstatsImport(list(input = input),
                          "MSstats", "DIANN")
    input = MSstatsClean(
        input, 
        MBR = MBR, 
        quantificationColumn = quantificationColumn,
        global_qvalue_cutoff = global_qvalue_cutoff, 
        qvalue_cutoff = qvalue_cutoff,
        pg_qvalue_cutoff = pg_qvalue_cutoff,
        calculateAnomalyScores = calculateAnomalyScores,
        anomalyModelFeatures = anomalyModelFeatures
    )
    #input = MSstatsMakeAnnotation(input, annotation)
    .writeChunkToFile(input, output_path, pos)
    NULL
}
