#' A subset of WorldClim bioclimatic variables
#'
#' A subset of WorldClim bioclimatic variables cropped on the Central and Western Europe.
#'
#' @docType data
#' @keywords datasets
#' @name Worldclim_tmp
#' @usage data(Worldclim_tmp)
#' @format A data frame obtained from a SpatRaster with 1080 rows, 2160 columns, and 6 layers, namely:"bio1"  "bio3"  "bio9"  "bio12" "bio13" "bio15"
#' @source \code{geodata::worldclim_global(var='bio', res=10, path=getwd())[[c(1, 3,9, 12, 13, 15)]]} 
'Worldclim_tmp'