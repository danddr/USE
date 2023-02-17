#' Sampling pseudo-absences for the training and testing datasets.  
#'
#' \code{paSampling} internally calls the \code{uniformSampling} function, which performs the uniform sampling of the environmental space. 
#' Before performing the uniform sampling of the environmental space, the environmental space is partitioned as follows: a kernel density estimate of the species presence observations identifies and excludes the portion of the environmental space associated with environmental conditions likely suitable for the species to establish using a threshold specified by the user. 
#' The bandwidth of the kernel can be estimated automatically or defined by the user, allowing to account for meta-population dynamics (e.g. sink-source populations).
#' 
#' @param env.rast A RasterStack, RasterBrick or a SpatRaster object comprising the variables describing the environmental space. 
#' @param pres A SpatialPointsDataframe, a SpatVector or an sf object including the presence-only observations of the species of interest.
#' @param thres (double) This value identifies the quantile value used to specify the boundary of the kernel density estimate (default \code{thres=0.75} ). Thus, probability values higher than the threshold should indicate portions of the multivariate space likely associated with presence points.
#' @param H The kernel bandwidth (i.e., the width of the kernel density function that defines its shape) excluding the portion of the environmental space associated with environmental conditions likely suitable for the species. It can be either defined by the user or automatically estimated by \code{paSampling} via \code{ks::Hpi}. 
#' @param grid.res (integer) resolution of the sampling grid. The resolution can be arbitrarily selected or defined using the \code{optimRes} function. 
#' @param n.tr (integer) number of pseudo-absences for the training dataset to sample in each cell of the sampling grid
#' @param n.ts (integer; optional) number of pseudo-absences for the testing dataset to sample in each cell of the sampling grid. sub.ts argument must be TRUE.
#' @param sub.ts (logical) sample the validation pseudo-absences
#' @param prev (double) prevalence value to be specified instead of n.tr and n.ts
#' @param plot_proc (logical) plot progress of the sampling, default FALSE
#' @param verbose (logical) Print verbose
#' @importFrom stats na.omit quantile
#' @return An sf object with the coordinates of the pseudo-absences both in the geographical and environmental space.
#' @export
paSampling <- function(env.rast, pres=NULL, thres=0.75, H=NULL, grid.res=NULL, n.tr = 5, sub.ts=FALSE, n.ts=5, prev=NULL, plot_proc=FALSE, verbose=FALSE) {
  if (!inherits(env.rast, "BasicRaster") && !inherits(env.rast, "SpatRaster") ) {
    stop("Environmental data provided in an unconvenient form")
  }
  if(is.null(pres)) {
    stop('Species occurrences must be provided by the user')
  }
  if (!inherits(pres, "SpatialPoints") && !inherits(pres, "SpatialPointsDataFrame") && !inherits(pres, "SpatVector") && !inherits(pres, "sf")){
    stop('Occurrences must be provided as spatial object')
  }
  if (is.null(grid.res)){
    stop('A grid resolution must be provided as length-one integer')
  }
  if (is.null(prev)){
    # n.tr= n.tr 
    # n.ts= n.ts
    estPrev=round(nrow(pres)/(n.tr*(grid.res^2)),2)
    message(paste("Estimated prevalence of", estPrev))
  } else {
    n.tr=(nrow(pres)/prev)/(grid.res^2)
    n.ts=(nrow(pres)/prev)/(grid.res^2)
    message(paste("User-defined prevalence of", prev))
  }
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
    } 
  if (!inherits(pres, "SpatVector")) {
    occ.vec <-terra::vect(pres)
  } else {
    occ.vec <- pres
  }
  message("Computing PCA and presences kernel density estimation in the multivariate space")
  grid.vec <-terra::as.points(env.rast)
  rpc <- RStoolbox::rasterPCA(raster::stack(env.rast), spca = TRUE)
  dt <- raster::as.data.frame(rpc$map[[c("PC1", "PC2")]], xy = TRUE)
  dt$myID <- seq_len(nrow(dt))
  
  id <- dt[,c("x", "y", "myID")]
  id_rast <- terra::rast(id, digits = 10, type="xyz")
  # id_rast<-terra::rast(id, resolution=terra::res(env.rast), digits = 10, type="xyz")
  id_rast <- terra::resample(id_rast, env.rast)
  terra::ext(id_rast) <- terra::ext(env.rast)
  
  PC12 <- c(terra::rast(rpc$map[[c("PC1", "PC2")]]), id_rast)
  abio.st <- c(id_rast, env.rast)
  abio.ex <- na.omit(terra::extract(abio.st, grid.vec, cells=FALSE, df=TRUE))
  PC12ex <- na.omit(terra::extract(PC12, grid.vec, cells=FALSE, df=TRUE))
  PC12occ <- terra::extract(PC12, occ.vec, cells=FALSE, df=TRUE)
 
  PC12ex <- merge(x=PC12ex, y=PC12occ, by='myID', all.x=TRUE )
  PC12ex$PA <- ifelse(is.na(PC12ex$ID.y),0,1)
  PC12ex <- PC12ex[, c("ID.x", "PC1.x", "PC2.x", "myID", "PA" )]
  names(PC12ex) <- c("ID", "PC1", "PC2" ,"myID", "PA")
  PC12ex <- na.omit(PC12ex)
  
  if (is.null(H)) {
    H <- ks::Hpi(x = PC12ex[,c( "PC1", "PC2")])
  } 
  
  estimate <- data.frame(KDE=ks::kde(PC12ex[,c( "PC1", "PC2")], 
                                 eval.points = PC12ex[,c( "PC1", "PC2")],h=H)$estimate, PC12ex[,-1])
  quantP <- quantile(estimate[estimate$PA == '1',"KDE"], thres)
  estimate$percP <- ifelse(estimate$KDE <= unname(quantP[1]),'out','in')
  
  fullDB.sp <- merge(x=estimate, y=dt[,c("x", "y", "myID")], by='myID', all.x=TRUE )
  fullDB.sp <- merge(x=fullDB.sp, y=abio.ex[,2:ncol(abio.ex)], by='myID', all.x=TRUE )
  fullDB.sp <- subset(fullDB.sp, fullDB.sp$PA == 0 & fullDB.sp$percP == "out" )
  fullDB.sp <- na.omit(fullDB.sp)
  fullDB.sp <- sf::st_as_sf(fullDB.sp, coords = c("PC1", "PC2"))
  
  if(is.null(prev)) {
  myPas <- NULL  
  } else {
    myPas <- floor(nrow(pres)/prev) 
  }
  
  message("\nPerforming pseudo-absences sampling in the environmental space\n")
  Res <- uniformSampling(sdf = fullDB.sp, 
                         grid.res=grid.res,  
                         n.tr = n.tr, 
                         n.prev = myPas, 
                         sub.ts = sub.ts, 
                         n.ts = n.ts,
                         plot_proc = plot_proc, 
                         verbose=verbose)
 
  if(sub.ts) {
    message("\n", paste(nrow(Res$obs.tr), "training pseudo-absences sampled in the environmental space, \n and", nrow(Res$obs.ts), "testing pseudo-absences sampled in the environmental space.",  sep = " "), "\n")  
    message("\n",paste("Estimated final prevalence of", round(nrow(pres)/nrow(Res$obs.tr),2), "instead of", prev ), "\n")
    } else {
      message("\n",paste(nrow(Res), "training pseudo-absences sampled in the environmental space", sep = " "), "\n") 
      message("\n",paste("Estimated final prevalence of", round(nrow(pres)/nrow(Res),2), "instead of", prev ),"\n")
    }
   return(Res)
}
