#' Get optimal resolution of the sampling grid 
#'
#' This function finds the optimal resolution to perform the uniform environmental sampling. 
#' It requires an sf object having point geometry given by the PC-scores values (the first two axes of a PCA computed on the environmental variables of interest), 
#' and a set of proposed resolution for the sampling grid to test.
#' For a given resolution, the function iterates across each cell and: i) finds the centroid of the PC-scores within the cell; 
#' ii) computes the average (squared) euclidean distance between each PC-score and their centroid.
#' Finally, the cell-specific distances are averaged by the number of cells.  
#' It returns a list with: i) a matrix reporting the values of the function to be minimized, along with the corresponding resolution; ii) the optimal resolution.
#' The function also provides a plot displaying the function values for each resolution, allowing the user to select a resolution on their own. 
#'
#' @param sdf an sf object having point geometry given by the PC-scores values
#' @param grid.res (integer) a vector of resolutions to be tested, i.e seq(1,100, by=1)
#' @param perc.thr rate of change (expressed in percentage) of the function to be minimized for selecting the optimal resolution.  
#' @param showOpt (logical) plot the result. 
#' @return It returns a list with: i) a matrix reporting the values of the function to be minimized, along with the corresponding resolution; ii) the optimal resolution.
#' @export
optim_res <- function(sdf, grid.res, perc.thr = 10, showOpt = TRUE) {
  if(!require(sf)) install.packages('sf')
  stopifnot(exprs = {
    is.numeric(perc.thr)
    is.numeric(grid.res)
    is.logical(showOpt)
    inherits(sdf, "sf")
  })
  grid.res <- sort(grid.res, decreasing = FALSE)
  SS_vec <- sapply(grid.res, function(res) { 
    Grd <- st_make_grid(sdf, n = res)
    SS_mean <- mean(sapply(Grd, function(i) {
      X <- sdf[i, ]
      if(nrow(X) >= 2) {
        X_c <- st_centroid(st_convex_hull(st_union(X)))
        D <- sum(as.numeric(st_distance(X, X_c))^2)/nrow(X) 
        return(D)
      } else {
        return(NA_real_)
      }
    }), na.rm = TRUE)
    return(SS_mean)
  })
  SS_vec <- sqrt(SS_vec)
  Rat_chg <- SS_vec[1:(length(SS_vec)-1)]/SS_vec[2:length(SS_vec)]
  Res_opt <- grid.res[(which(Rat_chg <= ((perc.thr/100)+1))[1])] 
  if(showOpt) {
    plot(grid.res, SS_vec, type = "b", xlab = "Grid resolution", ylab = "Function", main = "Optimal grid resolution")
    abline(v = Res_opt, col = "red")
  }
  return(list(F_val = cbind(Fval = SS_vec, Res = grid.res), Opt_res = Res_opt))
}