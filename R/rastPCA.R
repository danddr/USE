#' Principal Component Analysis for Rasters
#' 
#' The \code{rastPCA} function calculates the principal component analysis  (PCA) for SpatRaster, RasterBrick, or RasterStack objects and returns a SpatRaster with multiple layers representing the PCA components. Internally, \code{rastPCA} utilizes the \link[stats]{princomp} function for R-mode PCA analysis. The covariance matrix is computed using all the observations within the provided SpatRaster object, which describes the environmental conditions. 
#' The covariance matrix obtained is subsequently utilized as input for the \code{princomp} function, which conducts the PCA. The resulting PCA components are then used to generate the final SpatRaster, consisting of multiple layers that represent the PCA components.
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
    maskNA <-!sum(terra::app(env.rast, is.na))
    env.rast <- terra::mask(env.rast, maskNA, maskvalue = NA)
  }
  
  covMatrix <- terra::layerCor(env.rast, fun = "cov", na.rm = TRUE)
  eigenDecomp <- princompCustom(covmat = covMatrix[[1]], cor = stand)
  eigenDecomp$center <- covMatrix$mean[, 1] # modified after probable update of terra
  eigenDecomp$n.obs  <- terra::global(!any(is.na(env.rast)), sum)$sum
  if (stand==TRUE) {
    S <- diag(covMatrix$covariance)
    eigenDecomp$scale <- sqrt(S * (eigenDecomp$n.obs - 1)/eigenDecomp$n.obs)
  }
  pci <- terra::predict(env.rast, eigenDecomp, nPC=nPC, fun=pca_predict)
  names(pci) <- paste0("PC", 1:nPC)
  return(list(call = match.call(), pca = eigenDecomp, PCs = pci))
}