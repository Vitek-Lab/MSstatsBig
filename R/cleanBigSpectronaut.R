#' @export
DiskFrameSetup = function() {
  disk.frame::setup_disk.frame()
  options(future.globals.maxSize = Inf)
}


#' Clean larger-than-memory Spectronaut output
#'
#' @param file_path path to the Spectronaut file
#' @param output_directory path to disk.frame output
#' @param intensity name of the intensity column
#' @param filter_by_excluded logical, if TRUE, only rows with values of column
#' F.ExcludedFromQuantification equal to FALSE will be retained
#' @param filter_by_identified logical, if TRUE, only rows with values of column
#' EG.Identified will be retained
#' @param filter_by_qvalue logical, if TRUE, rows with qvalues in EG.Qvalue
#' and PG.Qvalue columns larger than `qvalue_cutoff` will be replaced
#' by 0 and NA, respectively
#' @param qvalue_cutoff numeric, value that will be used to filter by q-value
#' @param heavy_only logical, if TRUE, only data for heavy-labeled peptides
#' will be saved
#' @param col_Names optional character vector of column names, required if
#' columns of the output are not named
#'
#' @import data.table
#' @export
#'
cleanBigSpectronaut = function(file_path, output_directory,
                               intensity, filter_by_excluded,
                               filter_by_identified,
                               filter_by_qvalue = TRUE,
                               qvalue_cutoff = 0.01,
                               heavy_only = FALSE,
                               col_names = NULL) {
  FFrgLossType = FExcludedFromQuantification = NULL

  disk.frame::csv_to_disk.frame(
    file_path, outdir = output_directory,
    inmapfn = function(x) {
      if (is.null(col_names)) {
        current_cols = colnames(x)
      } else {
        current_cols = col_names
      }
      current_cols = MSstatsConvert:::.standardizeColnames(current_cols)
      colnames(x) = current_cols
      cols = c("R.FileName", "R.Condition", "R.Replicate",
               "PG.ProteinAccessions", "EG.ModifiedSequence", "FG.LabeledSequence",
               "FG.Charge", "F.FrgIon", "F.Charge",
               "EG.Identified", "F.ExcludedFromQuantification", "F.FrgLossType",
               "PG.Qvalue", "EG.Qvalue", "F.NormalizedPeakArea", "F.MeasuredRelativeIntensity",
               "F.PeakArea",
               "F.MassAccuracyPPM", "FG.FWHM", "EG.ApexRT", "FG.ShapeQualityScore")
      cols = MSstatsConvert:::.standardizeColnames(cols)

      # is_present = cols %in% current_cols
      # cols = intersect(cols, current_cols)
      x = x[, cols]
      x = data.table::as.data.table(x)
      # Handle PGQvalue or PGqvalue !!
      # Handle FCharge that could be FFrgZ !! (or not anymore???)
      data.table::setnames(x,
                           c("Run", "Condition", "BioReplicate", "ProteinName",
                             "PeptideSequence", "LabeledSequence", "PrecursorCharge", "FragmentIon",
                             "ProductCharge", "Identified", "Excluded",
                             "FFrgLossType", "PGQvalue", "EGQvalue",
                             "Intensity", "MeasuredRelativeIntensity", "PeakArea",
                             "MassAccuracyPPM", "FWHM", "ApexRt", "ShapeQualityScore"))#[is_present])
      x[, Intensity := as.numeric(Intensity)]
      x[, PeptideSequence := stringi::stri_replace_all(PeptideSequence,
                                                       "", fixed = "_")]
      x[, FragmentIon := stringi::stri_replace_all(FragmentIon,
                                                   "", fixed = "_")]
      if (is.character(x[["Excluded"]])) {
        x[, Excluded := !(Excluded == "False")]
      }
      if (filter_by_excluded) {
        x = x[!(Excluded)]
        x[, Excluded := NULL]
      }
      if (filter_by_identified) {
        x = x[(Identified)]
        x[, Identified := NULL]
      }
      if (filter_by_qvalue) {
        x[, Intensity := ifelse(EGQvalue < qvalue_cutoff, Intensity, 0)]
        x[, Intensity := ifelse(PGQvalue < qvalue_cutoff, Intensity, NA)]
      }
      x = x[FFrgLossType == "noloss"]
      x[, IsLabeled := grepl("Lys8", LabeledSequence) | grepl("Arg10", LabeledSequence)]
      x[, IsotopeLabelType := ifelse(IsLabeled, "H", "L")]
      psm_data = unique(x[, .(ProteinName, PeptideSequence, PrecursorCharge, FragmentIon,
                              ProductCharge, IsotopeLabelType, Run, BioReplicate, Condition,
                              Intensity, MeasuredRelativeIntensity, PeakArea, MassAccuracyPPM, FWHM, ApexRt,
                              ShapeQualityScore, EGQvalue, PGQvalue)])
      if (heavy_only) {
        psm_data = psm_data[IsotopeLabelType == "H"]
      }
      psm_data
    }, shardby = c("ProteinName"), overwrite = TRUE, backend = "readr",
    chunk_reader = "readr", delim = "\t",
    in_chunk_size = 1e6)
  TRUE
}


#' Save cleanBigSpectronaut output to a file
#'
#' @param folder_path path to disk.frame output of cleanBigSpectronaut
#' @param output_file path to output file
#' @param heavy_only logical, if TRUE, only heavy-labeled peptides will be retained
#'
#' @export
#'
saveSpectronautCSV = function(folder_path, output_file, heavy_only = FALSE) {
  disk.frame::setup_disk.frame(workers = 1)
  disk_frame = disk.frame::disk.frame(folder_path)
  disk.frame::cmap(disk_frame, function(chunk) {
    if (heavy_only) {
      chunk = chunk[chunk$IsotopeLabelType == "H", ]
    }
    data.table::fwrite(chunk, output_file, append = TRUE)
    NULL
  }, lazy=FALSE)
  disk.frame::setup_disk.frame() # turn multi worker back on
  TRUE
}
