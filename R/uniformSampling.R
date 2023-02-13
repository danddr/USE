#' Uniform sampling of the environmental space
#'
#'\code{uniformSampling} performs the uniform sampling of the environmental space.
#'It requires an sf object of geometry type “POINT”. Points are the scores of the first two axes of a principal component analysis performed on the correlation matrix of a set of environmental variables. The resolution of the sampling grid used to systematically collect observations within the environmental space must be provided. The resolution can be either set to an arbitrary value by the user or selected using the \code{optimRes()} function.
#'\code{uniformSampling} allows collecting observations to populate both training and (optional) a testing dataset.
#'In both cases, the user must provide a number of points to be sampled in each cell of the sampling grid (\code{n.tr:} points for the training dataset; \code{n.ts}: points for the testing dataset).
#'Note that the number of pseudo-absences eventually sampled depends on the spatial configuration of the points within the environmental space. Usually, some cells of the sampling grid are empty (i.e., do not include points), as points are not covering the whole extent of the environmental space. For this reason, the number of pseudo-absences eventually collected is likely to be lower than the product between the resolution of the sampling grid and \code{n.tr}(or \code{n.ts}). 
#'
#' @param sdf an sf object having point geometry given by the PC-scores values
#' @param grid.res (integer) resolution of the sampling grid. The resolution can be arbitrarily selected or defined using the \code{optimRes()} function. 
#' @param n.tr (integer) number of points for the training dataset to sample in each cell of the sampling grid
#' @param n.tr (integer; optional) number of expected points given a certain prevalence threshold for the training dataset.
#' @param n.ts (integer; optional) number of  points for the testing dataset to sample in each cell of the sampling grid. sub.ts argument must be TRUE.
#' @param n.prev (double) sample prevalence
#' @param sub.ts (logical) sample the validation points
#' @param plot_proc (logical) plot progress of the sampling
#' @param verbose (logical) Print verbose
#' @return A spatial point data frame with the coordinates of the sampled points both in the geographical and environmental space
#' @export
uniformSampling <- function(sdf, grid.res, n.tr = 5, n.prev = NULL, sub.ts = FALSE, n.ts = 5, plot_proc = FALSE, verbose = FALSE) {
  if(!(all(sf::st_is(sdf, "POINT")))) {
    stop("sdf object must have a sf POINT GEOMETRY class")
  }
  if(!is.numeric(n.tr)) stop(paste(n.tr, "is not of class 'numeric'.", sep = " "))
  if(!is.logical(plot_proc)) stop("plot_proc is not of class 'logical'; it has class 'numeric'.")
  grid <- st_make_grid(sdf, n = grid.res)
  sdf$ID <- row.names(sdf)
  res <- do.call(rbind, lapply(seq_len(length(grid)), function(i) {
    if(isTRUE(verbose)) message(paste("Processing tile", i, sep = " "))
    if(isTRUE(plot_proc)) {
      if(i == 1) {
        plot(grid, border = "black")
        plot(grid[i], col = "green", add = TRUE)
      } else {
        plot(grid[i], col = "green", add = TRUE)
      }
    }
    subs <- sdf[grid[i], ]
    if(nrow(subs) <= n.tr) {
      return(subs)
    } else {
      subs <- subs[sample(nrow(subs), n.tr, replace = FALSE), ]
      return(subs)
    }
  }))
  if(!is.null(n.prev)) {
    if(nrow(res) < n.prev) {
      n.cell <- length(grid)
      dif.abs <- n.prev - nrow(res) 
      if(dif.abs < length(grid)) {
        Subs <- sdf[!(sdf$ID %in% res$ID), ] 
        while(dif.abs > 0) { 
          ID_cell <- sample(n.cell, size = 1) 
          abs <- Subs[grid[ID_cell], ] 
          if(nrow(abs) == 0) next 
          abs <- abs[sample(nrow(abs), 1), ] 
          absID <- abs$ID #get the obs ID
          dif.abs <- (dif.abs - 1) 
          Subs <- Subs[!(Subs$ID %in% absID), ] 
          res <- rbind(res, abs) 
          if(nrow(Subs) == 0) { 
            message("There are not points left in the environmental space to reach prevalence")
            break 
          }
        }
      } else {
        Subs <- sdf[!(sdf$ID %in% res$ID), ] 
        ratio <- floor(n.prev/length(grid)) 
        while(dif.abs > 0) { 
          if(ratio == 0 ) { 
            message("There are not points left in the environmental space to reach prevalence")
            break
          }
          new_set <- do.call(rbind, lapply(seq_len(n.cell), function(.) { 
            subs <- Subs[grid[.], ]
            if(nrow(subs) == 0) {
              return(subs)
            } else {
              subs <- subs[sample(nrow(subs), 1), ]
              return(subs)
            }
          }))
          Subs <- Subs[!(Subs$ID %in% new_set$ID), ]
          res <- rbind(res, new_set)
          ratio <- (ratio - 1)
          dif.abs <- (dif.abs-nrow(new_set))
        }
      }
    } else {
      message("You should increase prevalence by randomly excluding some absence from the output dataset")
    }
  }
  if(isTRUE(sub.ts)) {
    abs_val <- sdf[!(sdf$ID %in% res$ID), ]
    res_val <- do.call(rbind, lapply(rev(seq_len(length(grid))), function(i) {
      if(isTRUE(verbose)) message(paste("Processing tile", i, sep = " "))
      if(isTRUE(plot_proc)) plot(grid[i], col = "red", add = TRUE)
      subs <- abs_val[grid[i], ]
      if(nrow(subs) <= n.ts) {
        return(subs)
      } else {
        subs <- subs[sample(nrow(subs), n.ts, replace = FALSE), ]
        return(subs)
      }
    }))
    return(list(Bkg.tr = res, Bkg.ts = res_val))
  } else {
    return(res)
  }
}