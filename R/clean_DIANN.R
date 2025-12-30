#' @keywords internal
reduceBigDIANN <- function(input_file, output_path, MBR = TRUE,
                           quantificationColumn = "FragmentQuantCorrected") {
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
                                                  pos)
  readr::read_delim_chunked(input_file,
                            readr::DataFrameCallback$new(diann_chunk),
                            delim = delim,
                            chunk_size = 1e6)
}

#' @keywords internal
cleanDIANNChunk = function(input, output_path, MBR, quantificationColumn, pos) {
  
  # 1. Select required columns
  base_cols <- c('Protein.Names', 'Stripped.Sequence', 'Modified.Sequence', 
                 'Precursor.Charge', quantificationColumn, 'Q.Value', 
                 'Precursor.Mz', 'Fragment.Info', 'Run')
  
  mbr_cols <- if (MBR) {
    c('Lib.Q.Value', 'Lib.PG.Q.Value')
  } else {
    c('Global.Q.Value', 'Global.PG.Q.Value')
  }
  
  req_cols <- intersect(c(base_cols, mbr_cols), colnames(input))
  input <- dplyr::select(input, all_of(req_cols))
  
  # 2. Split concatenated values (un-nest)
  split_cols <- intersect(c(quantificationColumn, "Fragment.Info"), colnames(input))
  if (length(split_cols) > 0) {
    input <- tidyr::separate_rows(input, all_of(split_cols), sep = ";")
  }
  
  # 3. Process fragment information
  input[[quantificationColumn]] <- as.numeric(input[[quantificationColumn]])
  
  input <- dplyr::mutate(
    input,
    FragmentIon = sub('\\^\\.\\*', '', .data$Fragment.Info),
    ProductCharge = dplyr::if_else(
      grepl("/", .data$Fragment.Info),
      # Extract charge, default to 1 if parsing fails
      as.integer(stringr::str_extract(.data$Fragment.Info, "(?<=/)[0-9]+")),
      1L
    )
  )
  
  # 4. Clean and filter data
  input <- dplyr::filter(
    input,
    !grepl("NH3|H2O", .data$FragmentIon) & !is.na(.data[[quantificationColumn]])
  )
  
  # 5. Rename columns to MSstats standard
  input <- dplyr::rename_with(input, .fn = function(x) gsub("\\.", "", x))
  
  # Standardize column names
  old_names <- c('ProteinNames', 'StrippedSequence', 'ModifiedSequence',
                 'PrecursorCharge', gsub("\\.", "", quantificationColumn), 'QValue',
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

  input <- dplyr::rename(input, any_of(rename_map))
  
  # Final column selection for MSstats format
  msstats_cols <- c("ProteinName", "PeptideSequence", "PeptideModifiedSequence", "PrecursorCharge", 
                    "FragmentIon", "ProductCharge", "Run", "Intensity")
  
  #TODO: confirm with Tony -- are these three needed?
  
  # Add annotation columns if they exist
  if ("Condition" %in% colnames(input)) msstats_cols <- c(msstats_cols, "Condition")
  if ("BioReplicate" %in% colnames(input)) msstats_cols <- c(msstats_cols, "BioReplicate")
   
  # Add IsotopeLabelType, assuming Light for DIANN
  input$IsotopeLabelType <- "L"
  msstats_cols <- c(msstats_cols, "IsotopeLabelType")
  
  final_cols <- intersect(msstats_cols, colnames(input))
  input <- dplyr::select(input, all_of(final_cols))
  
  # Write to file
  if (!is.null(pos)) {
    if (pos == 1) {
      readr::write_csv(input, file = output_path, append = FALSE)
    } else {
      readr::write_csv(input, file = output_path, append = TRUE)
    }
  }
  NULL
}