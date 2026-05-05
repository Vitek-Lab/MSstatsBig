#' General converter for larger-than-memory csv files in MSstats format 10-column format
#'
#' @param input_file name of the input text file in 10-column MSstats format.
#' @param output_file_name name of an output file which will be saved after pre-processing
#' @param backend "arrow" or "sparklyr". Option "sparklyr" requires a spark installation
#' and connection to spark instance provided in the `connection` parameter.
#' @param max_feature_count maximum number of features per protein. Features will
#' be selected based on highest average intensity.
#' @param filter_unique_peptides If TRUE, shared peptides will be removed.
#' Please refer to the `Details` section for additional information.
#' @param aggregate_psms If TRUE, multiple measurements per PSM in a Run will
#' be aggregated (by taking maximum value). Please refer to the `Details` section for additional information.
#' @param filter_few_obs If TRUE, feature with less than 3 observations across runs will be removed.
#' Please refer to the `Details` section for additional information.
#' @param remove_annotation If TRUE, columns BioReplicate and Condition will be removed
#' to reduce output file size. These will need to be added manually later before
#' using dataProcess function. Only applicable to sparklyr backend.
#' @param calculateAnomalyScores If TRUE, will carry anomaly model features through pipeline
#' @param anomalyModelFeatures Character vector of column names to be carried through the pipeline
#' @param connection Connection to a spark instance created with the
#' `spark_connect` function from `sparklyr` package.
#'
#' @details Filtering and aggregation may be very time consuming and the ability
#' to perform them in a given R session depends on available memory, settings of
#' external packages, etc. Hence, all value of related parameters (`filter_unique_peptides`,
#' `aggregate_psms`, `filter_few_obs`) are set to FALSE by default and only feature
#' selection is performed, which saves both computation time and memory.
#' Appropriately configured spark backend provides the most consistent way to
#' perform these operations.
#'
#' @return either arrow object or sparklyr table that can be optionally collected
#' into memory by using dplyr::collect function.
#'
#' @examples
#' converted_data <- bigFragPipetoMSstatsFormat(
#'   system.file("extdata", "fgexample.csv", package = "MSstatsBig"),
#'   "tencol_format.csv",
#'   backend="arrow")
#' procd <- MSstatsPreprocessBig("tencol_format.csv", "proc_out.csv", backend = "arrow")
#' head(dplyr::collect(procd))
#'
#' @export
#'
MSstatsPreprocessBig <-  function(input_file,
                                 output_file_name,
                                 backend,
                                 max_feature_count =  100,
                                 filter_unique_peptides =  FALSE,
                                 aggregate_psms =  FALSE,
                                 filter_few_obs =  FALSE,
                                 remove_annotation =  FALSE,
                                 calculateAnomalyScores = FALSE, 
                                 anomalyModelFeatures = c(),
                                 connection =  NULL) {
  if (backend == "arrow") {
    MSstatsPreprocessBigArrow(input_file,
                              output_file_name,
                              max_feature_count,
                              filter_unique_peptides,
                              aggregate_psms,
                              filter_few_obs,
                              calculateAnomalyScores, 
                              anomalyModelFeatures)
  } else if (backend == "sparklyr") {
    MSstatsPreprocessBigSparklyr(connection, input_file, output_file_name,
                                 max_feature_count, filter_unique_peptides,
                                 aggregate_psms, filter_few_obs,
                                 remove_annotation)
  } else {
    stop("backend not implemented")
  }
}

#' Convert out-of-memory FragPipe files to MSstats format.
#'
#' @inheritParams MSstatsPreprocessBig
#'
#' @export
#'
#' @examples
#' converted_data <- bigFragPipetoMSstatsFormat(
#'   system.file("extdata", "fgexample.csv", package = "MSstatsBig"),
#'   "output_file.csv",
#'   backend = "arrow")
#' converted_data <- dplyr::collect(converted_data)
#' head(converted_data)
#'
#' @return either arrow object or sparklyr table that can be optionally collected
#' into memory by using dplyr::collect function.
#'
bigFragPipetoMSstatsFormat <-  function(input_file, output_file_name,
                                       backend,
                                       max_feature_count =  100,
                                       filter_unique_peptides =  FALSE,
                                       aggregate_psms =  FALSE,
                                       filter_few_obs =  FALSE,
                                       remove_annotation =  FALSE,
                                       connection =  NULL) {
  MSstatsPreprocessBig(
    input_file = input_file, 
    output_file_name = output_file_name,
    backend = backend, 
    max_feature_count = max_feature_count, 
    filter_unique_peptides = filter_unique_peptides,
    aggregate_psms = aggregate_psms, 
    filter_few_obs = filter_few_obs, 
    remove_annotation = remove_annotation,
    connection = connection)
}


#' Convert out-of-memory Spectronaut files to MSstats format.
#'
#' @inheritParams MSstatsPreprocessBig
#' @param intensity Name of the intensity column to be used in Spectronaut
#' @param filter_by_excluded if TRUE, will filter by the `F.ExcludedFromQuantification` column.
#' @param filter_by_identified if TRUE, will filter by the `EG.Identified` column.
#' @param filter_by_qvalue if TRUE, will filter by EG.Qvalue and PG.Qvalue columns.
#' @param qvalue_cutoff cutoff which will be used for q-value filtering.
#'
#' @export
#'
#' @examples
#' converted_data <- bigSpectronauttoMSstatsFormat(
#'   system.file("extdata", "spectronaut_input.csv", package = "MSstatsBig"),
#'   "output_file.csv",
#'   backend="arrow")
#' converted_data <- dplyr::collect(converted_data)
#' head(converted_data)
#'
#' @return either arrow object or sparklyr table that can be optionally collected
#' into memory by using dplyr::collect function.
#'
bigSpectronauttoMSstatsFormat <-  function(input_file, output_file_name,
                                          backend,
                                          intensity = "F.NormalizedPeakArea",
                                          filter_by_excluded = FALSE,
                                          filter_by_identified = FALSE,
                                          filter_by_qvalue = FALSE,
                                          qvalue_cutoff = 0.01,
                                          max_feature_count = 100,
                                          filter_unique_peptides =  FALSE,
                                          aggregate_psms =  FALSE,
                                          filter_few_obs =  FALSE,
                                          remove_annotation =  FALSE,
                                          calculateAnomalyScores=FALSE, 
                                          anomalyModelFeatures=c(),
                                          connection =  NULL) {
  reduced_file <- .prefixedPath("reduce_output_", output_file_name)
  reduceBigSpectronaut(input_file, reduced_file,
                       intensity, filter_by_excluded, filter_by_identified,
                       filter_by_qvalue, qvalue_cutoff,
                       calculateAnomalyScores, anomalyModelFeatures)
  msstats_data <- MSstatsPreprocessBig(
    input_file = reduced_file,
    output_file_name = output_file_name, 
    backend = backend, 
    max_feature_count = max_feature_count,
    filter_unique_peptides = filter_unique_peptides,
    aggregate_psms = aggregate_psms, 
    filter_few_obs = filter_few_obs, 
    remove_annotation = remove_annotation, 
    calculateAnomalyScores = calculateAnomalyScores, 
    anomalyModelFeatures = anomalyModelFeatures, 
    connection = connection)
  
  return(msstats_data)
  
}


#' Convert out-of-memory DIANN files to MSstats format.
#'
#' @inheritParams MSstatsPreprocessBig
#' @inheritParams MSstatsConvert::DIANNtoMSstatsFormat
#'
#' @export
#'
#' @return either arrow object or sparklyr table that can be optionally collected
#' into memory by using dplyr::collect function.
#'
bigDIANNtoMSstatsFormat <- function(input_file, 
                                    annotation = NULL,
                                    output_file_name,
                                    backend,
                                    MBR = TRUE,
                                    quantificationColumn = "FragmentQuantCorrected",
                                    global_qvalue_cutoff = 0.01,
                                    qvalue_cutoff = 0.01,
                                    pg_qvalue_cutoff = 0.01,
                                    max_feature_count = 100,
                                    filter_unique_peptides =  FALSE,
                                    aggregate_psms =  FALSE,
                                    filter_few_obs =  FALSE,
                                    remove_annotation =  FALSE,
                                    calculateAnomalyScores=FALSE, 
                                    anomalyModelFeatures=c(),
                                    connection =  NULL) {
  
  # Reduce and clean the DIANN report file in chunks
  reduced_file <- .prefixedPath("reduce_output_", output_file_name)
  reduceBigDIANN(input_file,
                 reduced_file,
                 MBR,
                 quantificationColumn,
                 global_qvalue_cutoff, qvalue_cutoff, pg_qvalue_cutoff,
                 calculateAnomalyScores, anomalyModelFeatures,
                 annotation)

  reduced <- arrow::open_dataset(reduced_file, format = "csv")

  # Identify columns where Arrow inferred 'null' type (all values NA)
  null_cols <- names(reduced$schema)[
    vapply(reduced$schema$fields, function(f) f$type$ToString() == "null", logical(1))
  ]

  if (length(null_cols) > 0) {
    # Drop null-typed columns using a lazy select (no data loaded into memory)
    reduced <- dplyr::select(reduced, -dplyr::all_of(null_cols))

    # Write back using Arrow's streaming writer — stays out-of-memory.
    # write_dataset creates a directory, but open_dataset can read
    # directories just as easily as single files.
    cleaned_file <- .prefixedPath("cleaned_", output_file_name)
    arrow::write_dataset(reduced, cleaned_file, format = "csv")
    reduced_file <- cleaned_file
  }

  # Preprocess the cleaned data (feature selection, etc.)
  msstats_data <- MSstatsPreprocessBig(
    input_file = reduced_file,
    output_file_name = output_file_name,
    backend = backend,
    max_feature_count = max_feature_count,
    filter_unique_peptides = filter_unique_peptides,
    aggregate_psms = aggregate_psms,
    filter_few_obs = filter_few_obs,
    remove_annotation = remove_annotation,
    calculateAnomalyScores = calculateAnomalyScores,
    anomalyModelFeatures = anomalyModelFeatures,
    connection = connection)

  # Merge annotation with the preprocessed data and persist the merge so
  # callers reopening output_file_name see Condition/BioReplicate. The arrow
  # rewrite stays lazy — the underlying source is reduced_file, not
  # output_file_name, so we can safely overwrite the directory we just wrote.
  if (!is.null(annotation)) {
    msstats_data <- MSstatsAddAnnotationBig(msstats_data, annotation)
    if (backend == "arrow") {
      unlink(output_file_name, recursive = TRUE, force = TRUE)
      arrow::write_dataset(msstats_data, output_file_name, format = "csv")
    }
  }
  return(msstats_data)
}


#' Merge annotation to output of MSstatsPreprocessBig
#'
#' @param input output of MSstatsPreprocessBig
#' @param annotation run annotation
#'
#' @export
#'
#' @examples
#' converted_data <- bigFragPipetoMSstatsFormat(
#'   system.file("extdata", "fgexample.csv", package = "MSstatsBig"),
#'   "output_file.csv",
#'   backend = "arrow")
#' converted_data <- dplyr::collect(converted_data)
#' head(converted_data)
#' # Change annotation as an example:
#' converted_data$Condition <- NULL
#' converted_data$BioReplicate <- NULL
#' annot <- data.frame(Run = unique(converted_data[["Run"]]))
#' annot$BioReplicate <- rep(1:53, times = 2)
#' annot$Condition <- rep(1:2, each = 53)
#' head(MSstatsAddAnnotationBig(converted_data, annot))
#'
#' @importFrom MSstats dataProcess groupComparison
#' @importFrom utils head sessionInfo
#'
#' @return table of `input` and `annotation` merged by Run column.
#'
MSstatsAddAnnotationBig <- function(input, annotation) {
  join_keys <- "Run"
  
  # Use tbl_vars which works reliably on both Arrow
  # datasets, arrow_dplyr_query objects, and data frames
  input_cols <- dplyr::tbl_vars(input)
  
  overlap_cols <- setdiff(
    intersect(input_cols, colnames(annotation)),
    join_keys
  )
  if (length(overlap_cols) > 0) {
    input <- dplyr::select(input, -dplyr::all_of(overlap_cols))
  }
  
  dplyr::inner_join(input, annotation, by = join_keys)
}
