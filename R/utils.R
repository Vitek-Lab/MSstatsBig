#' Build an intermediate output path by prefixing only the basename.
#'
#' Naive `paste0(prefix, output_file_name)` corrupts paths that contain a
#' directory (`subdir/out.csv` → `topN_subdir/out.csv`,
#' `/tmp/out.csv` → `topN_/tmp/out.csv`). Splitting via dirname/basename keeps
#' the directory component intact so intermediate files land beside the final
#' output.
#'
#' @param prefix Character scalar prepended to the basename.
#' @param path  Output file path supplied by the caller.
#' @return Character scalar.
#' @keywords internal
.prefixedPath <- function(prefix, path) {
    file.path(dirname(path), paste0(prefix, basename(path)))
}

#' Write chunk to file
#'
#' @param input Data frame
#' @param output_path Path to output file
#' @param pos Chunk position
#' @return NULL
#' @keywords internal
.writeChunkToFile <- function(input, output_path, pos) {
    # Write to file
    if (!is.null(pos)) {
        if (pos == 1) {
            readr::write_csv(input, file = output_path, append = FALSE)
        } else {
            readr::write_csv(input, file = output_path, append = TRUE)
        }
    }
}