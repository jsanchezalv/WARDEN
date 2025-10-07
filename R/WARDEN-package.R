## usethis namespace: start
#' @useDynLib WARDEN, .registration = TRUE
#' @importFrom Rcpp evalCpp
## usethis namespace: end
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("WARDEN", libpath)
}
