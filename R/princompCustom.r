#' Custom version of princomp
#' The warning() at L53 substitutes the stop() in the original version of "princomp".
#' @param subset an optional vector used to select rows (observations) of the data matrix x.
#' @param x a numeric matrix or data frame which provides the data for the principal components analysis.
#' @param cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.)
#' @param scores a logical value indicating whether the score on each principal component should be calculated.
#' @param covmat a covariance matrix, or a covariance list as returned by cov.wt (and cov.mve or cov.mcd from package MASS). If supplied, this is used rather than the covariance matrix of x.
#' @param fix_sign Should the signs of the loadings and scores be chosen so that the first element of each loading is non-negative?
#' @importFrom stats cov.wt setNames
#' @keywords internal
#' @export
princompCustom <- function (x, cor = FALSE, scores = TRUE, covmat = NULL, subset = rep_len(TRUE, nrow(as.matrix(x))), fix_sign = TRUE, ...) {
  chkDots(...)
  cl <- match.call()
  cl[[1L]] <- as.name("princomp")
  z <- if (!missing(x)) 
    as.matrix(x)[subset, , drop = FALSE]
  if (is.list(covmat)) {
    if (anyNA(match(c("cov", "n.obs"), names(covmat)))) 
      stop("'covmat' is not a valid covariance list")
    cv <- covmat$cov
    n.obs <- covmat$n.obs
    cen <- covmat$center
  }
  else if (is.matrix(covmat)) {
    if (!missing(x)) 
      warning("both 'x' and 'covmat' were supplied: 'x' will be ignored")
    cv <- covmat
    n.obs <- NA
    cen <- NULL
  }
  else if (is.null(covmat)) {
    dn <- dim(z)
    if (dn[1L] < dn[2L]) 
      stop("'princomp' can only be used with more units than variables")
    covmat <- stats::cov.wt(z)
    n.obs <- covmat$n.obs
    cv <- covmat$cov * (1 - 1/n.obs)
    cen <- covmat$center
  }
  else stop("'covmat' is of unknown type")
  if (!is.numeric(cv)) 
    stop("PCA applies only to numerical variables")
  if (cor) {
    sds <- sqrt(diag(cv))
    if (any(sds == 0)) 
      stop("cannot use 'cor = TRUE' with a constant variable")
    cv <- cv/(sds %o% sds)
  }
  edc <- eigen(cv, symmetric = TRUE)
  ev <- edc$values
  if (any(neg <- ev < 0)) {
    if (any(ev[neg] < -9 * .Machine$double.eps * ev[1L])) 
      warning("covariance matrix is not non-negative definite")
    else ev[neg] <- 0
  }
  cn <- paste0("Comp.", 1L:ncol(cv))
  names(ev) <- cn
  dimnames(edc$vectors) <- if (missing(x)) 
    list(dimnames(cv)[[2L]], cn)
  else list(dimnames(x)[[2L]], cn)
  sdev <- sqrt(ev)
  sc <- stats::setNames(if (cor) 
    sds
    else rep.int(1, ncol(cv)), colnames(cv))
  fix <- if (fix_sign) 
    function(A) {
      mysign <- function(x) ifelse(x < 0, -1, 1)
      A[] <- apply(A, 2L, function(x) x * mysign(x[1L]))
      A
    }
  else identity
  ev <- fix(edc$vectors)
  scr <- if (scores && !missing(x) && !is.null(cen)) 
    scale(z, center = cen, scale = sc) %*% ev
  if (is.null(cen)) 
    cen <- rep(NA_real_, nrow(cv))
  edc <- list(sdev = sdev, loadings = structure(ev, class = "loadings"), 
              center = cen, scale = sc, n.obs = n.obs, scores = scr, 
              call = cl)
  class(edc) <- "princomp"
  edc
}

pca_predict <- function (data, model, nPC) 
{
    predict(data, model)[, 1:nPC]
}

rastPCA <- function (env.rast, nPC = NULL, naMask = TRUE, stand = FALSE) 
{
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  }
  if (terra::nlyr(env.rast) <= 1) {
    stop("At least two layers are needed to calculate PCA")
  }
  if (is.null(nPC)) {
    nPC <- terra::nlyr(env.rast)
  }
  if (nPC > terra::nlyr(env.rast)) {
    nPC <- terra::nlyr(env.rast)
    message(paste0("\nThe maximum number of PCs that can be estimated is ", 
                   terra::nlyr(env.rast), "\n"))
  }
  if (sum(terra::global(env.rast, fun = "isNA")) == sum(terra::global(env.rast, 
                                                                      fun = function(x) terra::ncell(x)))) {
    stop("The layers are empty or contain only NAs")
  }
  if (naMask == TRUE) {
    maskNA <- !sum(terra::app(env.rast, is.na))
    env.rast <- terra::mask(env.rast, maskNA, maskvalue = FALSE)
  }
  covMatrix <- terra::layerCor(env.rast, fun = "cov", na.rm = TRUE)
  eigenDecomp <- princompCustom(covmat = covMatrix$covariance, cor = stand)
  eigenDecomp$center <- covMatrix$mean
  eigenDecomp$n.obs <- nrow(as.data.frame(env.rast[[1]]))
  if (stand == TRUE) {
    S <- diag(covMatrix$covariance)
    eigenDecomp$scale <- sqrt(S * (eigenDecomp$n.obs - 1)/eigenDecomp$n.obs)
  }
  pci <- terra::predict(env.rast, eigenDecomp, nPC = nPC, fun = pca_predict)
  names(pci) <- paste0("PC", 1:nPC)
  return(list(call = match.call(), pca = eigenDecomp, PCs = pci))
}