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