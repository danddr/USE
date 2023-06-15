#' Sampling pseudo-absences for the training and testing datasets.  
#'
#' \code{paSampling} performs a two-step procedure for uniformly sampling pseudo-absences within the environmental space.
#' In the initial step, a kernel-based filter is utilized to determine the subset of the environmental space that will be subsequently sampled. The kernel-based filter calculates the probability function based on the presence observations, enabling the identification of areas within the environmental space that likely exhibit suitable conditions for the species. To achieve this, a probability threshold value is utilized to assign observations to the corresponding portion of the environmental space. These areas, deemed to have suitable environmental conditions, are excluded from the subsequent uniform sampling process conducted in the second step using the \code{uniformSampling} function, which is internally called. 
#' The bandwidth of the kernel can be automatically estimated from the presence observations or directly set by the user, providing flexibility in determining the scope and precision of the filter.
#' 
#' Being designed with species distribution models in mind, \code{paSampling} allows collectively sampling pseudo-absences for both the training and testing dataset (optional). In both cases, the user must provide a number of observations that will be sampled in each cell of the sampling grid (\code{n.tr}: points for the training dataset; \code{n.ts}: points for the testing dataset). 
#' Note that the optimal resolution of the sampling grid can be found using the \code{optimRes} function. Also, note that the number of pseudo-absences eventually sampled in each cell by the internally-called \code{uniformSampling} function depends on the spatial configuration of the observations within the environmental space. Indeed, in most cases some cells of the sampling grid will be empty (i.e., those at the boundary of the environmental space). For this reason, the number of pseudo-absences returned by \code{paSampling} is likely to be lower than the product between the number of cells of the sampling gird and \code{n.tr}(or \code{n.ts}).
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
#' 
paSampling <- function (env.rast=NULL, pres = NULL, thres = 0.75, H = NULL, grid.res = NULL, 
                         n.tr = 5, sub.ts = FALSE, n.ts = 5, prev = NULL, plot_proc = FALSE, 
                         verbose = FALSE) 
{
  if (!inherits(env.rast, "BasicRaster") && !inherits(env.rast, 
                                                      "SpatRaster")) {
    stop("Environmental data provided in an unconvenient form")
  }
  if (is.null(pres)) {
    stop("Species occurrences must be provided by the user")
  }
  if (!inherits(pres, "SpatialPoints") && !inherits(pres, "SpatialPointsDataFrame") && 
      !inherits(pres, "SpatVector") && !inherits(pres, "sf")) {
    stop("Occurrences must be provided as spatial object")
  }
  if (is.null(grid.res)) {
    stop("A grid resolution must be provided as length-one integer")
  }
  if (is.null(prev)) {
    estPrev <- round(nrow(pres)/(n.tr * (grid.res^2)), 2)
    message(paste("Estimated prevalence of", estPrev))
  }
  else {
    n.tr <- (nrow(pres)/prev)/(grid.res^2)
    n.ts <- (nrow(pres)/prev)/(grid.res^2)
    message(paste("User-defined prevalence of", prev))
  }
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  }
  if (!inherits(pres, "SpatVector")) {
    occ.vec <- terra::vect(pres)
  }
  else {
    occ.vec <- pres
  }
  message("Computing PCA and presences kernel density estimation in the multivariate space")
  grid.vec <- terra::as.points(env.rast)
  rpc <- rastPCA(env.rast, stand = TRUE)
  dt <- terra::as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE)
  dt$myID <- seq_len(nrow(dt))
  id <- dt[, c("x", "y", "myID")]
  id_rast <- terra::rast(id, digits = 10, type = "xyz")
  id_rast <- terra::resample(id_rast, env.rast)
  terra::ext(id_rast) <- terra::ext(env.rast)
  PC12 <- c(rpc$PCs[[c("PC1", "PC2")]], id_rast)
  abio.st <- c(id_rast, env.rast)
  abio.ex <- na.omit(terra::extract(abio.st, grid.vec, cells = FALSE, 
                                    df = TRUE))
  PC12ex <- na.omit(terra::extract(PC12, grid.vec, cells = FALSE, 
                                   df = TRUE))
  PC12occ <- terra::extract(PC12, occ.vec, cells = FALSE, df = TRUE)
  PC12ex <- merge(x = PC12ex, y = PC12occ, by = "myID", all.x = TRUE)
  PC12ex$PA <- ifelse(is.na(PC12ex$ID.y), 0, 1)
  PC12ex <- PC12ex[, c("ID.x", "PC1.x", "PC2.x", "myID", "PA")]
  names(PC12ex) <- c("ID", "PC1", "PC2", "myID", "PA")
  PC12ex <- na.omit(PC12ex)
  if (is.null(H)) {
    H <- ks::Hpi(x = PC12ex[, c("PC1", "PC2")])
  }
  estimate <- data.frame(KDE = ks::kde(PC12ex[PC12ex$PA == 
                                                1, c("PC1", "PC2")], eval.points = PC12ex[PC12ex$PA == 
                                                                                            1, c("PC1", "PC2")], h = H)$estimate, PC12ex[PC12ex$PA == 
                                                                                                                                           1, -1])
  quantP <- quantile(estimate[, "KDE"], thres)
  estimate$percP <- ifelse(estimate$KDE <= unname(quantP[1]), "out", "in")
  estimate <- merge(PC12ex, estimate[estimate$PA == 1,c("myID","percP")], by = "myID", all.x = TRUE)
  chull <- sf::st_as_sf(subset(estimate, estimate$percP=="in", select=c( "PC1","PC2" )), coords=c( "PC1","PC2" )) 
  chull <- sf::st_union(chull)
  chull <- sf::st_convex_hull(chull)
  estimate <- sf::st_as_sf(estimate, coords=c( "PC1","PC2" ))
  chull <- sf::st_filter(estimate, chull)
  chull$percP <- "in"
  chull <- cbind.data.frame(sf::st_drop_geometry(chull), 
                           data.frame("PC1"=sf::st_coordinates(chull)[,1], 
                                      "PC2"=sf::st_coordinates(chull)[,2]))
  estimate<- cbind.data.frame(sf::st_drop_geometry(estimate), 
                           data.frame("PC1"=sf::st_coordinates(estimate)[,1], 
                                      "PC2"=sf::st_coordinates(estimate)[,2]))
  estimate<-estimate[ !estimate$myID %in% chull$myID, ]
  estimate$percP <- "out"
  estimate<-rbind.data.frame(estimate, chull)
  fullDB.sp <- merge(x = estimate, y = dt[, c("x", "y", "myID")], 
                     by = "myID", all.x = TRUE)
  fullDB.sp <- merge(x = fullDB.sp, y = abio.ex[, 2:ncol(abio.ex)], 
                     by = "myID", all.x = TRUE)
  fullDB.sp <- na.omit(fullDB.sp)
  fullDB.sp <- sf::st_as_sf(fullDB.sp, coords = c("PC1", "PC2"))
  if (is.null(prev)) {
    myPas <- NULL
  }
  else {
    myPas <- floor(nrow(pres)/prev)
  }
  message("\nPerforming pseudo-absences sampling in the environmental space\n")
  fullDB.sp <- subset(fullDB.sp, fullDB.sp$percP == "out" & fullDB.sp$PA == 0)
  Res <- uniformSampling(sdf = fullDB.sp, grid.res = grid.res, 
                         n.tr = n.tr, n.prev = myPas, sub.ts = sub.ts, n.ts = n.ts, 
                         plot_proc = plot_proc, verbose = verbose)
  if (sub.ts) {
    message("\n", paste(nrow(Res$obs.tr), "training pseudo-absences sampled in the environmental space, \n and", 
                        nrow(Res$obs.ts), "testing pseudo-absences sampled in the environmental space.", 
                        sep = " "), "\n")
    message("\n", paste("Estimated final prevalence of", 
                        round(nrow(pres)/nrow(Res$obs.tr), 2), "instead of", 
                        prev), "\n")
  }
  else {
    message("\n", paste(nrow(Res), "training pseudo-absences sampled in the environmental space", 
                        sep = " "), "\n")
    message("\n", paste("Estimated final prevalence of", 
                        round(nrow(pres)/nrow(Res), 2), "instead of", prev), 
            "\n")
  }
  return(Res)
}
