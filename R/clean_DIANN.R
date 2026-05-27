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
  # Per-chunk callback shared by both the parquet and delimited-text paths.
  # `pos` drives .writeChunkToFile: pos == 1 overwrites, pos > 1 appends.
  diann_chunk <- function(x, pos) cleanDIANNChunk(x, output_path, MBR,
                                                  quantificationColumn, pos,
                                                  global_qvalue_cutoff,
                                                  qvalue_cutoff,
                                                  pg_qvalue_cutoff,
                                                  calculateAnomalyScores,
                                                  anomalyModelFeatures,
                                                  annotation)

  # Parquet branch (DIANN 2.0+): stream record batches via arrow so the file
  # is never fully materialised. read_delim_chunked can't read parquet bytes.
  if (tolower(tools::file_ext(input_file)) == "parquet") {
    # Lazy handle to the parquet file — no data loaded yet.
    ds <- arrow::open_dataset(input_file, format = "parquet")
    # Scanner + RecordBatchReader yields one batch at a time on demand.
    # batch_size matches the delimited-text path's 1M-row chunks; row-group
    # boundaries in the parquet may cap individual batches below this.
    scanner <- arrow::Scanner$create(ds, batch_size = 1e6)
    reader <- scanner$ToRecordBatchReader()
    pos <- 1
    repeat {
      batch <- reader$read_next_batch()
      if (is.null(batch)) break  # exhausted
      # Materialise just this batch, run it through the shared cleaner.
      diann_chunk(as.data.frame(batch), pos)
      pos <- pos + 1
    }
    return(invisible(NULL))
  }

  # Delimited-text branch (DIANN 1.x TSV/CSV): sniff the delimiter from the
  # first line, defaulting to tab when nothing matches.
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

  # Stream the file in 1M-row chunks, invoking diann_chunk for each.
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
    if (!is.element("IsotopeLabelType", colnames(input))) {
        input <- dplyr::mutate(input, IsotopeLabelType = "L")
    }
    #input = MSstatsMakeAnnotation(input, annotation)
    .writeChunkToFile(input, output_path, pos)
    NULL
}
