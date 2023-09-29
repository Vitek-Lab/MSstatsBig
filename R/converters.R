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
#' converted_data = BigFragPipetoMSstatsFormat(
#'   system.file("extdata", "fgexample.csv", package = "MSstatsBig"),
#'   "tencol_format.csv",
#'   backend="arrow")
#' procd = MSstatsPreprocessBig("tencol_format.csv", "proc_out.csv", backend = "arrow")
#' head(dplyr::collect(procd))
#'
#' @export
#'
MSstatsPreprocessBig = function(input_file,
                                output_file_name,
                                backend,
                                max_feature_count = 20,
                                filter_unique_peptides = FALSE,
                                aggregate_psms = FALSE,
                                filter_few_obs = FALSE,
                                remove_annotation = FALSE,
                                connection = NULL) {
  if (backend == "arrow") {
    MSstatsPreprocessBigArrow(input_file,
                              output_file_name,
                              max_feature_count,
                              filter_unique_peptides,
                              aggregate_psms,
                              filter_few_obs)
  } else if (backend == "sparklyr") {
    MSstatsPreprocessBigSparklyr(connection, input, output_file_name,
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
#' converted_data = BigFragPipetoMSstatsFormat(
#'   system.file("extdata", "fgexample.csv", package = "MSstatsBig"),
#'   "output_file.csv",
#'   backend="arrow")
#' converted_data = dplyr::collect(converted_data)
#' head(converted_data)
#'
#' @return either arrow object or sparklyr table that can be optionally collected
#' into memory by using dplyr::collect function.
#'
BigFragPipetoMSstatsFormat = function(input_file, output_file_name,
                                      backend,
                                      max_feature_count = 20,
                                      filter_unique_peptides = FALSE,
                                      aggregate_psms = FALSE,
                                      filter_few_obs = FALSE,
                                      remove_annotation = FALSE,
                                      connection = NULL) {
  MSstatsPreprocessBig(input_file, output_file_name,
                       backend, max_feature_count, filter_unique_peptides,
                       aggregate_psms, filter_few_obs, remove_annotation,
                       connection)
}


#' Convert out-of-memory Spectronaut files to MSstats format.
#'
#' @inheritParams MSstatsPreprocessBig
#' @param intensity name of a column that will be used as Intensity column.
#' @param filter_by_excluded if TRUE, will filter by the `F.ExcludedFromQuantification` column.
#' @param filter_by_identified if TRUE, will filter by the `EG.Identified` column.
#' @param filter_by_qvalue if TRUE, will filter by EG.Qvalue and PG.Qvalue columns.
#' @param qvalue_cutoff cutoff which will be used for q-value filtering.
#'
#' @export
#'
#' @return either arrow object or sparklyr table that can be optionally collected
#' into memory by using dplyr::collect function.
#'
BigSpectronauttoMSstatsFormat = function(input_file, output_file_name,
                                         backend,
                                         intensity = "F.PeakArea",
                                         filter_by_excluded = TRUE,
                                         filter_by_identified = TRUE,
                                         filter_by_qvalue = TRUE,
                                         qvalue_cutoff = 0.01,
                                         max_feature_count = 20,
                                         filter_unique_peptides = FALSE,
                                         aggregate_psms = FALSE,
                                         filter_few_obs = FALSE,
                                         remove_annotation = FALSE,
                                         connection = NULL) {
  reduceBigSpectronaut(input_file, paste0("reduce_output_", output_file_name), intensity,
                       filter_by_excluded, filter_by_identified,
                       filter_by_qvalue, qvalue_cutoff)
  MSstatsPreprocessBig(paste0("reduce_output_", input_file_path),
                       output_file_path, backend, max_feature_count,
                       aggregate_psms, filter_few_obs,
                       remove_annotation, connection)
}


#' Merge annotation to output of MSstatsPreprocessBig
#'
#' @param input output of MSstatsPreprocessBig
#' @param annotation run annotation
#'
#' @export
#'
#' @examples
#' converted_data = BigFragPipetoMSstatsFormat(
#'   system.file("extdata", "fgexample.csv", package = "MSstatsBig"),
#'   "output_file.csv",
#'   backend="arrow")
#' converted_data = dplyr::collect(converted_data)
#' head(converted_data)
#' # Change annotation as an example:
#' converted_data$Condition = NULL
#' converted_data$BioReplicate = NULL
#' annot = data.frame(Run = unique(converted_data[["Run"]]))
#' annot$BioReplicate = rep(1:53, times = 2)
#' annot$Condition = rep(1:2, each = 53)
#' head(MSstatsAddAnnotationBig(converted_data, annot))
#'
#' @return table of `input` and `annotation` merged by Run column.
#'
MSstatsAddAnnotationBig = function(input, annotation) {
  dplyr::inner_join(input, annotation, by = "Run")
}
