#' Principal Component Analysis for Rasters
#' 
#' Calculates R-mode PCA for SpatRaster, RasterBrick or RasterStack and returns a SpatRaster with multiple layers of PCA scores. 
#'  
#' Internally rastPCA relies on the use of \link[stats]{princomp} (R-mode PCA). The covariance matrix is computed using all the observations and will then be used to calculate princomp and predict the full raster.
#' 
#' Pixels with missing values in one or more bands will be set to NA. The built-in check for such pixels can lead to a slow-down of rastPCA.
#' However, if you make sure or know beforehand that all pixels have either only valid values or only NAs throughout all layers you can disable this check
#' by setting \code{naMask=FALSE} which speeds up the computation.
#' 
#' Standardized PCA (\code{stand=TRUE}) can be useful if imagery or bands of different dynamic ranges are combined. In this case, the correlation matrix is computed instead of the covariance matrix, which
#' has the same effect as using normalised bands of unit variance. 
#' 
#' @param env.rast  A RasterStack, RasterBrick or a SpatRaster object comprising the variables describing the environmental space.
#' @param nPC Integer. Number of PCA components to return.
#' @param stand Logical. If \code{TRUE}, perform standardized PCA. Corresponds to centered and scaled input image. This is usually beneficial for equal weighting of all layers. (\code{FALSE} by default)
#' @param naMask Logical. Masks all pixels which have at least one NA (default \code{TRUE} is recommended but introduces a slow-down. 
#' @seealso The \code{rastPCA} function has been conceptualized starting from \code{RStoolbox::rasterPCA} (\url{https://github.com/bleutner/RStoolbox}).
#' @importFrom stats princomp
#' @return Returns a named list containing the PCA model object ($pca) and the SpatRaster with the principal component layers ($PCs).
#' @export 

rastPCA <- function (env.rast, nPC=NULL, naMask=TRUE, stand=FALSE){
  
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  } 
  if (terra::nlyr(env.rast) <= 1) {
    stop("At least two layers are needed to calculate PCA")
  }
  if(is.null(nPC)){
    nPC <- terra::nlyr(env.rast)
    }
  if (nPC > terra::nlyr(env.rast)) {
    nPC <- terra::nlyr(env.rast)
    message(paste0( "\nThe maximum number of PCs that can be estimated is ", terra::nlyr(env.rast),'\n'))
  }
  
  if (sum(terra::global(env.rast, fun="isNA"))==sum(terra::global(env.rast,fun=function(x) terra::ncell(x)))) {
    stop("The layers are empty or contain only NAs")}
  
  if (naMask==TRUE) {
    maskNA <- is.na(env.rast[[1]])
    env.rast <- terra::mask(env.rast, maskNA, maskvalue = NA)
  }
  
  covMatrix <- terra::layerCor(env.rast, fun = "cov", na.rm = TRUE)
  eigenDecomp <- stats::princomp(covmat = covMatrix[[1]], cor = stand)
  eigenDecomp$center <- covMatrix$mean
  eigenDecomp$n.obs <- nrow(as.data.frame(env.rast[[1]]))
  if (stand==TRUE) {
    S <- diag(covMatrix$covariance)
    eigenDecomp$scale <- sqrt(S * (eigenDecomp$n.obs - 1)/eigenDecomp$n.obs)
  }
  pci <- terra::predict(env.rast, eigenDecomp, nPC=nPC, fun=pca_predict)
  names(pci) <- paste0("PC", 1:nPC)
  return(list(call = match.call(), pca = eigenDecomp, PCs = pci))
}