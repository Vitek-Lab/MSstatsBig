#' Read and clean a large DIANN file in chunks
#' 
#' @param input_file Path to the input DIANN file
#' @param output_path Path to the output CSV file
#' @param MBR Boolean, whether MBR was used
#' @param quantificationColumn Name of the column containing intensity values
#' @param global_qvalue_cutoff Global Q-value cutoff
#' @param qvalue_cutoff Q-value cutoff
#' @param pg_qvalue_cutoff Protein group Q-value cutoff
#' @return NULL. Writes to file.
#' @keywords internal
reduceBigDIANN <- function(input_file, output_path, MBR = TRUE,
                           quantificationColumn = "FragmentQuantCorrected",
                           global_qvalue_cutoff = 0.01,
                           qvalue_cutoff = 0.01,
                           pg_qvalue_cutoff = 0.01) {
  if (grepl("csv", input_file)) {
    delim = ","
  } else if (grepl("tsv|xls", input_file)) {
    delim = "\t"
  } else {
    delim <- ";"
  }
  
  diann_chunk <- function(x, pos) cleanDIANNChunk(x, output_path, MBR, quantificationColumn, pos,
                     global_qvalue_cutoff, qvalue_cutoff, pg_qvalue_cutoff)

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
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean
#' @return NULL
#' @keywords internal
cleanDIANNChunk = function(input, output_path, MBR, quantificationColumn, pos,
                           global_qvalue_cutoff = 0.01,
                           qvalue_cutoff = 0.01,
                           pg_qvalue_cutoff = 0.01) {
    input = MSstatsImport(list(input = input),
                          "MSstats", "DIANN")
    input = MSstatsClean(
        input, 
        MBR = MBR, 
        quantificationColumn = quantificationColumn,
        global_qvalue_cutoff = global_qvalue_cutoff, 
        qvalue_cutoff = qvalue_cutoff,
        pg_qvalue_cutoff = pg_qvalue_cutoff
    )
    .writeChunkToFile(input, output_path, pos)
    NULL
}
