#' Predict pca
#' @param data  A RasterStack, RasterBrick or a SpatRaster object comprising the variables describing the environmental space.
#' @param model \code{princomp} object.
#' @param nPC Integer. Number of PCA components to return. 
#' @importFrom stats predict
#' @keywords internal
#' @NoRd
pca_predict <- function(data, model, nPC) {
  predict(data, model)[,1:nPC]
}