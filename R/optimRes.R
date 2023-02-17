#' Get optimal resolution of the sampling grid 
#' 
#' \code{optimRes} identifies the optimal resolution of the sampling grid to be used to perform the uniform environmental sampling.
#' The function requires an sf object of geometry type “POINT”. Points here are the scores of the first two axes of a principal component analysis performed on the correlation matrix of a set of environmental variables. A set of candidate resolutions must be provided to find the optimal resolution of the sampling grid. For each candidate resolution, \code{optimRes} computes a metric summarizing the average squared Euclidean distance between points within each cell and the centroid of the convex hull built around the points. This metric is then compared across different sampling grids with increasing resolution (i.e., number of cells). The best resolution is selected as the one allowing finding the best trade-off between the number of cells and the average distance among points within cells. In other words, this is the resolution allowing for uniformly sampling the environmental space without overfitting it. This is obtained by choosing when the average distance among pseudo-absences cannot be reduced by more than 10% (10% is the default value, but another threshold can be set by the user).
#' \code{optimRes} returns a list of length 2 with a matrix reporting the measure computed for each sampling grid at the corresponding resolution and the selected optimal resolution.
#' The function also provides a plot displaying the function values for each resolution, allowing the user to select a resolution on their own.
#' 
#' @param sdf an sf object having point geometry given by the PC-scores values
#' @param grid.res (integer) a vector of resolutions to be tested, i.e seq(1,100, by=1)
#' @param perc.thr rate of change (expressed in percentage) of the function to be minimized for selecting the optimal resolution.  
#' @param cr (integer) number of cores for parallel computing. The default cluster type is PSOCK.
#' @param showOpt (logical) plot the result. 
#' @return It returns a list with: i) a matrix reporting the values of the function to be minimized, along with the corresponding resolution; ii) the optimal resolution.
#' @details In case the function returns NA as the optimal resolution: i) increase the range of \code{grid.res}, ii) increase \code{perc.thr}.
#' @export
optimRes <- function(sdf, grid.res, perc.thr = 10, cr = 1, showOpt = TRUE) {
  stopifnot(exprs = {
    is.numeric(perc.thr)
    is.numeric(grid.res)
    is.logical(showOpt)
    inherits(sdf, "sf")
  })
  if(is.null(cr)) stop("Please, provide a cluster")
  cl <- parallel::makeCluster(spec=cr, type="PSOCK", nnodes=cr, outfile="")
  grid.res <- sort(grid.res, decreasing = FALSE)
  SS_vec <- parallel::parSapply(cl, grid.res, function(res) {
    Grd <- sf::st_make_grid(sdf, n = res)
    SS_mean <- mean(sapply(Grd, function(i) {
      X <- sdf[i, ]
      if (nrow(X) >= 2) {
        X_c <- sf::st_centroid(sf::st_convex_hull(sf::st_union(X)))
        D <- sum(as.numeric(sf::st_distance(X, X_c))^2)/nrow(X)
        return(D)
      }
      else {
        return(NA_real_)
      }
    }), na.rm = TRUE)
    return(SS_mean)
  })
  parallel::stopCluster(cl) #stop clusters used in parSapply
  SS_vec <- sqrt(SS_vec)
  Rat_chg <- SS_vec[1:(length(SS_vec) - 1)]/SS_vec[2:length(SS_vec)]
  Res_opt <- grid.res[(which(Rat_chg <= ((perc.thr/100) + 1))[1])]
  if (showOpt) {
    plot(grid.res, SS_vec, type = "b", xlab = "Grid resolution", 
         ylab = "Function", main = "Optimal grid resolution")
    graphics::abline(v = Res_opt, col = "red")
  }
  return(list(F_val = cbind(Fval = SS_vec, Res = grid.res), 
              Opt_res = Res_opt))
}
