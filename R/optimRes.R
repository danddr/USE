#' Get optimal resolution of the sampling grid 
#' 
#' \code{optimRes} identifies the optimal resolution of the sampling grid to be used to perform the uniform environmental sampling. 
#' To find this optimal resolution, a set of candidate resolutions must be provided. For each candidate resolution, \code{optimRes} calculates a metric that summarizes the average squared Euclidean distance between the observations (PC-scores of the first two principal components) within each cell and the centroid of the convex hull encompassing the points. It's important to note that the centroid is specific to each cell.
#' 
#' This metric is then compared across different sampling grids with increasing resolution, i.e., an increasing number of cells. The best resolution is selected based on the trade-off between the number of cells and the average distance among observations within each cell. Essentially, the goal is to find the finest resolution of the sampling grid that enables uniform sampling of the environmental space without overfitting it.
#' 
#' By default, the optimal resolution is determined as the one where the average distance among observations and the cell-specific centroids cannot be reduced by more than 10%. However, users have the flexibility to adjust this setting according to their needs. The \code{optimRes} function returns a list with two elements. The first element is a matrix that reports the metric calculated for each sampling grid at the corresponding resolution. The second element is the selected optimal resolution.
#' 
#' Additionally, the function provides a plot that displays the metric values for each resolution. This allows users to visually analyze the relationship between resolution and the associated metric, thereby empowering them to make an informed decision when selecting a resolution.

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
