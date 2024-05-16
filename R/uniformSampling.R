#' Uniform sampling of the environmental space
#'
#'\code{uniformSampling} performs the uniform sampling of observations within the environmental space. Note that \code{uniformSampling} can be more generally used to sample observations (not necessarily associated with species occurrence data) within bi-dimensional spaces (e.g., vegetation plots). Being designed with species distribution models in mind, \code{uniformSampling} allows collectively sampling observations for both the training and testing dataset (optional). 
#'In both cases, the user must provide a number of observations that will be sampled in each cell of the sampling grid (\code{n.tr}: points for the training dataset; \code{n.ts}: points for the testing dataset). Note that the optimal resolution of the sampling grid can be found using the \code{optimRes} function. 

#' @param sdf an sf object having point geometry given by the PC-scores values
#' @param grid.res (integer) resolution of the sampling grid. The resolution can be arbitrarily selected or defined using the \code{optimRes()} function. 
#' @param n.tr (integer) number of points for the training dataset to sample in each cell of the sampling grid
#' @param n.tr (integer; optional) number of expected points given a certain prevalence threshold for the training dataset.
#' @param n.ts (integer; optional) number of  points for the testing dataset to sample in each cell of the sampling grid. sub.ts argument must be TRUE.
#' @param n.prev (double) sample prevalence
#' @param sub.ts (logical) sample the validation points
#' @param plot_proc (logical) plot progress of the sampling
#' @param verbose (logical) Print verbose
#' @return An sf object with the coordinates of the sampled points both in the geographical and environmental space
#' @export
uniformSampling <- function(sdf, grid.res, n.tr = 5, n.prev = NULL, sub.ts = FALSE, n.ts = 5, plot_proc = FALSE, verbose = FALSE) {
  if(!(all(sf::st_is(sdf, "POINT")))) {
    stop("sdf object must have a sf POINT GEOMETRY class")
  }
  if(!is.numeric(n.tr)) stop(paste(n.tr, "is not of class 'numeric'.", sep = " "))
  if(!is.logical(plot_proc)) stop("plot_proc is not of class 'logical'; it has class 'numeric'.")
  grid <- sf::st_make_grid(sdf, n = grid.res)
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
  if(any(duplicated(res$ID))){
    res <- res[!duplicated(res$ID),]
    warning("Warning: pseudo-replicate sampled and removed from the final training dataset ")
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
    if(any(duplicated(res_val$ID))){
      res_val <- res_val[!duplicated(res_val$ID),]
      warning("Warning: pseudo-replicate sampled and removed from the final testing dataset ")
    }
    return(list(obs.tr = res, obs.ts = res_val))
  } else {
    return(res)
  }
}