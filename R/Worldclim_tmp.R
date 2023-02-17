#' A subset of WorldClim bioclimatic variables
#'
#' A subset of WorldClim bioclimatic variables cropped on the Central and Western Europe.
#'
#' @docType data
#' @keywords datasets
#' @name Worldclim_tmp
#' @usage data(Worldclim_tmp)
#' @format A RasterLayer with 1080 rows, 2160 columns, and 5 layers, namely: wc2.1_10m_bio_4, wc2.1_10m_bio_3, wc2.1_10m_bio_14, wc2.1_10m_bio_9, wc2.1_10m_bio_15.
#' @source \code{geodata::worldclim_global(var='bio', res=10, path=getwd())[[c(4,  3, 14,  9, 15)]]} 
'Worldclim_tmp'