#' Uniform sampling of the environmental space
#'
#' This function allows to perform a uniform sampling of the environmental space. 
#' It requires an sf object having point geometry given by the PC-scores values (the first two axes of a PCA computed on the environmental variables of interest), and the 
#' bounding box of the environmental space (defined by the extent of the PC-scores values).
#' A resolution of the sampling grid must be provided by the user. The resolution can be arbitrarily selected or 
#' or defined using the Optim_res() function. 
#' The functions allows to retrieve both a training and a testing dataset (optional) of uniformely sampled background points. 
#' The user must provide suitable numbers of background points to be sampled in each cell of the sampling grid 
#' (n.tr: background points for the training dataset; n.ts: background points for the testing dataset). 
#' Notice that the final number of sampled background points depends on the PC-scores configuration in the environmental space. 
#'
#' @param sdf an sf object having point geometry given by the PC-scores values
#' @param grid.res (integer) resolution of the sampling grid. The resolution can be arbitrarily selected or defined using the Optim_res() function. 
#' @param n.tr (integer) number of background points for the training dataset to sample in each cell of the sampling grid
#' @param n.tr (integer; optional) number of expected background given a certain prevalence thresholdpoints for the training dataset
#' @param n.ts (integer; optional) number of background points for the testing dataset to sample in each cell of the sampling grid. sub.ts argument must be TRUE.
#' @param sub.ts (logical) sample the validation background points
#' @param plot_proc (logical) plot progress of the sampling
#' @param verbose (logical) Print verbose
#' @return A spatial point data frame with the coordinates of the sampled background points both in the geographical and environmental space
#' @export
uesampling <- function(sdf, grid.res, n.tr = 5, n.prev = NULL, sub.ts = FALSE, n.ts = 5, plot_proc = FALSE, verbose = FALSE) {
  if(!require(sf)) install.packages('sf')
  if(!require(assertive)) install.packages('assertive')
  if(!(all(st_is(sdf, "POINT")))) {
    stop("sdf object must have a sf POINT GEOMETRY class")
  }
  if(!is_numeric(n.tr)) stop(paste(n.tr, "is not of class 'numeric'.", sep = " "))
  if(!is_logical(plot_proc)) stop("plot_proc is not of class 'logical'; it has class 'numeric'.")
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