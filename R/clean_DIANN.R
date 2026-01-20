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
  
  diann_chunk <- function(x, pos) cleanDIANNChunk(x,
                                                  output_path,
                                                  MBR,
                                                  quantificationColumn,
                                                  pos,
                                                  global_qvalue_cutoff,
                                                  qvalue_cutoff,
                                                  pg_qvalue_cutoff)
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
#' @return NULL
#' @keywords internal
cleanDIANNChunk = function(input, output_path, MBR, quantificationColumn, pos,
                           global_qvalue_cutoff = 0.01,
                           qvalue_cutoff = 0.01,
                           pg_qvalue_cutoff = 0.01) {
  # 1. Handle "auto" quantification column
  processed <- .handleAutoQuantification(input, quantificationColumn)
  input <- processed$input
  quantificationColumn <- processed$quantificationColumn
  
  # 2. Select required columns
  input <- .selectDIANNColumns(input, MBR, quantificationColumn)
  input <- .cleanDIANNAddMissingColumns(input)
  
  # 3. Expand concatenated rows
  input <- .expandDIANNRows(input, quantificationColumn)
  
  # 4. Process fragment info (extract intensity, charge, ion)
  input <- .processDIANNFragmentInfo(input, quantificationColumn)
  
  # 5. Filter invalid fragments
  input <- .filterDIANNFragments(input, quantificationColumn)
  
  # 6. Standardize column names
  input <- .standardizeDIANNColumns(input, quantificationColumn)
  
  # 7. Filter by Q-values
  input <- .filterDIANNByQValues(input, MBR, global_qvalue_cutoff, qvalue_cutoff, pg_qvalue_cutoff)
  
  # 8. Finalize columns (select final set, add IsotopeLabelType)
  input <- .finalizeDIANNColumns(input)
  
  # 9. Write to output
  .writeDIANNChunk(input, output_path, pos)
  
  NULL
}

#' Handle automatic detection of quantification columns
#' 
#' @param input Data frame
#' @param quantificationColumn Name of column or "auto"
#' @return List with input data frame and updated quantification column name
#' @keywords internal
.handleAutoQuantification <- function(input, quantificationColumn) {
  if (quantificationColumn == "auto") {
    fragment_columns <- grep("^Fr[0-9]+Quantity$", colnames(input), value = TRUE)
    if (length(fragment_columns) == 0) {
      stop("No fragment quantification columns found. Please check your input.")
    }
    input <- tidyr::unite(input, "FragmentQuantCorrected", all_of(fragment_columns), sep = ";")
    quantificationColumn <- "FragmentQuantCorrected"
  }
  list(input = input, quantificationColumn = quantificationColumn)
}

#' Select required columns from DIANN output
#' 
#' @param input Data frame
#' @param MBR Boolean
#' @param quantificationColumn Name of intensity column
#' @return Data frame with selected columns
#' @keywords internal
.selectDIANNColumns <- function(input, MBR, quantificationColumn) {
  base_cols <- c('Protein.Names', 'Stripped.Sequence', 'Modified.Sequence', 
                 'Precursor.Charge', quantificationColumn, 'Q.Value', 
                 'Precursor.Mz', 'Fragment.Info', 'Run')
  
  mbr_cols <- if (MBR) {
    c('Lib.Q.Value', 'Lib.PG.Q.Value')
  } else {
    c('Global.Q.Value', 'Global.PG.Q.Value')
  }
  
  req_cols <- intersect(c(base_cols, mbr_cols), colnames(input))
  dplyr::select(input, all_of(req_cols))
}

#' Add missing required columns
#' 
#' @param input Data frame
#' @return Data frame with missing columns added
#' @keywords internal
.cleanDIANNAddMissingColumns <- function(input) {
  if (!"Precursor.Mz" %in% colnames(input)) {
    input <- dplyr::mutate(input, Precursor.Mz = NA)
  }
  if (!"Fragment.Info" %in% colnames(input)) {
    input <- dplyr::mutate(input, Fragment.Info = NA)
  }
  input
}

#' Expand rows with multiple fragments
#' 
#' @param input Data frame
#' @param quantificationColumn Name of intensity column
#' @return Data frame with expanded rows
#' @keywords internal
.expandDIANNRows <- function(input, quantificationColumn) {
  split_cols <- intersect(c(quantificationColumn, "Fragment.Info"), colnames(input))
  if (length(split_cols) > 0) {
    tidyr::separate_rows(input, all_of(split_cols), sep = ";")
  } else {
    input
  }
}

#' Process fragment information strings
#' 
#' @param input Data frame
#' @param quantificationColumn Name of intensity column
#' @return Data frame with FragmentIon and ProductCharge columns
#' @keywords internal
.processDIANNFragmentInfo <- function(input, quantificationColumn) {
  # Convert Intensity to Numeric from Char strings
  input[[quantificationColumn]] <- as.numeric(input[[quantificationColumn]])
  
  dplyr::mutate(
    input,
    FragmentIon = sub('\\^\\.\\*', '', .data$Fragment.Info),
    
    # Extract product charge
    ProductCharge = dplyr::if_else(
      grepl("/", .data$Fragment.Info),
      # Extract charge (number right after "/" in string), default to 1 if parsing fails
      as.integer(stringr::str_extract(.data$Fragment.Info, "(?<=/)[0-9]+")),
      1L,
      missing = 1L
    )
  )
}

#' Filter invalid fragments
#' 
#' @param input Data frame
#' @param quantificationColumn Name of intensity column
#' @return Filtered data frame
#' @keywords internal
.filterDIANNFragments <- function(input, quantificationColumn) {
  dplyr::filter(
    input,
    (!grepl("NH3|H2O", .data$FragmentIon) | is.na(.data$FragmentIon)) & !is.na(.data[[quantificationColumn]])
  )
}

#' Standardize column names to MSstats format
#' 
#' @param input Data frame
#' @param quantificationColumn Name of intensity column
#' @return Data frame with renamed columns
#' @keywords internal
.standardizeDIANNColumns <- function(input, quantificationColumn) {
  input <- dplyr::rename_with(input, .fn = function(x) gsub("\\.", "", x))
  
  # Standardize column names
  clean_quant_col <- gsub("\\.", "", quantificationColumn)
  old_names <- c('ProteinNames', 'StrippedSequence', 'ModifiedSequence',
                 'PrecursorCharge', clean_quant_col, 'QValue',
                 'PrecursorMz', 'FragmentIon', 'Run', 'ProductCharge')
  new_names <- c('ProteinName', 'PeptideSequence', 'PeptideModifiedSequence',
                 'PrecursorCharge', 'Intensity', 'DetectionQValue', 
                 'PrecursorMz', 'FragmentIon', 'Run', 'ProductCharge')
  
  current_names <- colnames(input)
  names_to_rename <- intersect(current_names, old_names)
  
  # Create a named vector for renaming in the format c(new_name = old_name)
  new_names_subset <- new_names[match(names_to_rename, old_names)]
  rename_map <- setNames(names_to_rename, new_names_subset)
  rename_map <- rename_map[!is.na(names(rename_map))]
  
  dplyr::rename(input, any_of(rename_map))
}

#' Filter data by Q-values
#' 
#' @param input Data frame
#' @param MBR Boolean
#' @param global_qvalue_cutoff Numeric
#' @param qvalue_cutoff Numeric
#' @param pg_qvalue_cutoff Numeric
#' @return Filtered data frame
#' @keywords internal
.filterDIANNByQValues <- function(input, MBR, global_qvalue_cutoff, qvalue_cutoff, pg_qvalue_cutoff) {
  input <- dplyr::filter(input, DetectionQValue < global_qvalue_cutoff)
  
  if (MBR) {
    dplyr::filter(input, LibPGQValue < pg_qvalue_cutoff & LibQValue < qvalue_cutoff)
  } else {
    dplyr::filter(input, GlobalPGQValue < pg_qvalue_cutoff & GlobalQValue < qvalue_cutoff)
  }
}

#' Finalize columns for output
#' 
#' @param input Data frame
#' @return Data frame with final columns
#' @keywords internal
.finalizeDIANNColumns <- function(input) {
  # Final column selection for MSstats format
  msstats_cols <- c("ProteinName", "PeptideSequence", "PeptideModifiedSequence", "PrecursorCharge", 
                    "FragmentIon", "ProductCharge", "Run", "Intensity")
  
  
  # Add annotation columns if they exist
  if ("Condition" %in% colnames(input)) msstats_cols <- c(msstats_cols, "Condition")
  if ("BioReplicate" %in% colnames(input)) msstats_cols <- c(msstats_cols, "BioReplicate")
  
  # Add IsotopeLabelType, assuming Light for DIANN
  input$IsotopeLabelType <- "L"
  msstats_cols <- c(msstats_cols, "IsotopeLabelType")
  
  final_cols <- intersect(msstats_cols, colnames(input))
  dplyr::select(input, all_of(final_cols))
}

#' Write chunk to file
#' 
#' @param input Data frame
#' @param output_path Path to output file
#' @param pos Chunk position
#' @return NULL
#' @keywords internal
.writeDIANNChunk <- function(input, output_path, pos) {
  # Write to file
  if (!is.null(pos)) {
    if (pos == 1) {
      readr::write_csv(input, file = output_path, append = FALSE)
    } else {
      readr::write_csv(input, file = output_path, append = TRUE)
    }
  }
}