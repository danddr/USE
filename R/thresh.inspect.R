#' Inspect the effect of the kernel threshold parameter on the environmental space partitioning
#'
#'\code{thresh.inspect} function allows for a pre-inspection of the impact that selecting a specific threshold for the kernel-based filter will have on the exclusion of the environmental space in the subsequent uniform sampling of the pseudo-absences process (see \code{paSampling}). By providing a range of threshold values, the function generates a plot that illustrates the entire environmental space, including the portion delineated by the kernel-based filter and the associated convex-hull. This plot helps visualize the areas that will be excluded from the uniform sampling of the pseudo-absences.
#'This functionality proves particularly valuable in determining a meaningful threshold for the kernel-based filter in specific ecological scenarios. For instance, when dealing with sink populations, selecting the appropriate threshold enables the exclusion of environmental space regions where the species is present, but the conditions are unsuitable. This allows for a more accurate sampling of pseudo-absences, considering the unique requirements of different ecological contexts.
#' @param env.rast A RasterStack, RasterBrick or a SpatRaster object comprising the variables describing the environmental space. 
#' @param pres A SpatialPointsDataframe, a SpatVector or an sf object including the presence-only observations of the species of interest.
#' @param thres (double) This value or vector of values identifies the quantile value used to specify the boundary of the kernel density estimate (default \code{thres=0.75} ). Thus, probability values higher than the threshold should indicate portions of the multivariate space likely associated with presence points.
#' @param H The kernel bandwidth (i.e., the width of the kernel density function that defines its shape) excluding the portion of the environmental space associated with environmental conditions likely suitable for the species. It can be either defined by the user or automatically estimated by \code{paSampling} via \code{ks::Hpi}. 
#' @importFrom stats na.omit quantile
#' @return A ggplot2 object showing how the environmental space is partitioned accordingly to the selected \code{thres} values.
#' @export
#' 
# thresh.inspect(
# env.rast=env_datasets[[2]], pres=myPres, thres=c(0.1,0.5), H=NULL)
thresh.inspect <- function (env.rast, pres = NULL, thres = 0.75, H = NULL) 
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
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  }
  if (!inherits(pres, "SpatVector")) {
    occ.vec <- terra::vect(pres)
  }
  else {
    occ.vec <- pres
  }
  cols <- c("In" = "darkorange", "Out" = "lightblue", "Presences" = "black")
  sizes <- c("In" = 0.5, "Out" = 0.5, "Presences" = 0.2)
  
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
  estimate_list <- vector("list", length(thres))
  for(i in 1:length(thres)){
      estimate.i <- estimate
      quantP <- quantile(estimate.i[, "KDE"], thres[i])
      estimate.i$percP <- ifelse(estimate.i$KDE <= unname(quantP[1]), "out", "in")
      estimate.i <- merge(PC12ex, estimate.i[estimate.i$PA == 1,c("myID","percP")], by = "myID", all.x = TRUE)
      chull <- sf::st_as_sf(subset(estimate.i, estimate.i$percP=="in", select=c( "PC1","PC2" )), coords=c( "PC1","PC2" )) 
      chull <- sf::st_union(chull)
      chull <- sf::st_convex_hull(chull)
      estimate.i <- sf::st_as_sf(estimate.i, coords=c( "PC1","PC2" ))
      chull <- sf::st_filter(estimate.i, chull)
      chull$percP <- "In"
      chull <- cbind.data.frame(sf::st_drop_geometry(chull), 
                                data.frame("PC1"=sf::st_coordinates(chull)[,1], 
                                           "PC2"=sf::st_coordinates(chull)[,2]))
      estimate.i <- cbind.data.frame(sf::st_drop_geometry(estimate.i), 
                                  data.frame("PC1"=sf::st_coordinates(estimate.i)[,1], 
                                             "PC2"=sf::st_coordinates(estimate.i)[,2]))
      estimate.i <- estimate.i[ !estimate.i$myID %in% chull$myID, ]
      estimate.i$percP <- "Out"
      estimate.i<-rbind.data.frame(estimate.i, chull)
      estimate.i <- estimate.i[, c("percP", "PC1", "PC2")]
      estimate.i <- rbind(estimate.i, data.frame(PC1=PC12occ$PC1, 
                                                 PC2=PC12occ$PC2, 
                                                 percP="Presences"))
      estimate.i$thres <- thres[i]
      estimate_list[[i]] <- estimate.i
    }
  estimate_list <- do.call(rbind.data.frame, estimate_list)
  p <- ggplot2::ggplot(estimate_list, ggplot2::aes(x = estimate_list$PC1, y = estimate_list$PC2, 
                                          col = estimate_list$percP, size = estimate_list$percP))+
    ggplot2::geom_point()+
    ggplot2::scale_colour_manual(values = cols)+
    ggplot2::scale_size_manual(values = sizes, guide = 'none')+
    ggplot2::labs(x = "PC1", y = "PC2", col = "Kernel threshold")+
    ggplot2::facet_wrap(~estimate_list$thres)+
    ggplot2::theme_classic()+
    ggplot2::theme(plot.title =  ggplot2::element_text(size=16, face = 'bold'),
          legend.background = ggplot2::element_blank(),
          strip.text.x = ggplot2::element_text(size=12, face = 'bold'),
          legend.position = 'bottom',
          text = ggplot2::element_text(size=16), 
          panel.grid = ggplot2::element_blank(),
          aspect.ratio = 1,
          legend.text = ggplot2::element_text(size=16, angle = 0), legend.title = ggplot2::element_text(size=16),
          legend.key.size = ggplot2::unit(2, 'cm'))
  return(p)
}
