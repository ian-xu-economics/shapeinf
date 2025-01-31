#' @title Package Load Hook
#' @description This function is called when the package is loaded. It registers the S3 method for printing objects of class \code{et}.
#' @param libname The library name where the package is installed.
#' @param pkgname The name of the package.
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  registerS3method("print", "shapeinf", print.shapeinf)
  registerS3method("summary", "shapeinf", summary.shapeinf)
}