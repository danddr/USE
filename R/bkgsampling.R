#' Background points training and testing dataset sampling 
#'
#' This function internally calls the uesampling() function to perform the uniform sampling of the environmental space. 
#' Before performing the uniform sampling of the environmental space, a set of PC-scores within a core area associated with the species presences 
#' is filtered out using a kernel density estimation approach. 
#' The bandwith of the kernel can be either defined by the user or automatically estimated, allowing the user to account for metapopulation dynamics (e.g. sink-source populations).
#'
#' @param env.rast A raster stack or raster brick object comprising the variables describing the environmental space. 
#' @param pres A SpatialPointDataFrame including species of interest presence-only.
#' @param thres Kernel density threshold SPIEGARE
#' @param H Parameter defining kernel bandwidth, it can be provided by the user or automatically estimated via SPIEGARE
#' @param grid.res (integer) resolution of the sampling grid. The resolution can be arbitrarily selected or defined using the Optim_res() function. 
#' @param n.tr (integer) number of background points for the training dataset to sample in each cell of the sampling grid
#' @param n.ts (integer; optional) number of background points for the testing dataset to sample in each cell of the sampling grid. sub.ts argument must be TRUE.
#' @param sub.ts (logical) sample the validation background points
#' @param prev (double) prevalence value to be specified instead of n.tr and n.ts
#' @param plot_proc (logical) plot progress of the sampling, default FALSE
#' @return A spatial point data frame with the coordinates of the background points both in the geographical and environmental space.
#' @export
bkgsampling <- function(env.rast, pres=NULL, thres=0.75, H=NULL, grid.res=NULL, n.tr = 5, sub.ts=FALSE, n.ts=5, prev=NULL, plot_proc=FALSE) {
  if(!require(raster)) install.packages('raster')
  if(!require(RStoolbox)) install.packages('RStoolbox')
  if(!require(terra)) install.packages('terra')
  if(!require(tidyverse)) install.packages('tidyverse') 
  if(!require(data.table)) install.packages('data.table')
  if(!require(ks)) install.packages('ks')  
  if(!require(sf)) install.packages('sf')  
  if (!inherits(env.rast, "BasicRaster")) {
    stop("Data provided in a non convenient from")
  }
  if(is.null(pres)) {
    stop('Species occurrences must be provided by the user')
  }
  if (!inherits(pres, "SpatialPoints") && !inherits(pres, "SpatialPointsDataFrame")){
    stop('Occurrences must be provided as spatial object')
  }
  
  if (is.null(grid.res)){
    stop('A grid resolution must be provided as length-one integer')
  }
  if (is.null(prev)){
    n.tr= n.tr 
    n.ts= n.ts
    estPrev=round(nrow(pres)/(n.tr*(grid.res^2)),2)
    message(paste("Estimated prevalence of", estPrev))
  } else {
    n.tr=(nrow(pres)/prev)/(grid.res^2)
    n.ts=(nrow(pres)/prev)/(grid.res^2)
    message(paste("User-defined prevalence of", prev))
  }
  message("Computing PCA and presences kernel density estimation in the PC-space")
  grid.vec<-vect(rasterToPoints(env.rast, spatial = TRUE))
  pen.vec<-vect(pres)
  rpc<-rasterPCA(env.rast,spca = TRUE)
  dt <- data.table(raster::as.data.frame(rpc$map[[c("PC1", "PC2")]], xy = TRUE))
  dt$myID<-seq_len(nrow(dt))
  
  id<-dt[,c("x", "y", "myID")]
  id_rast<-rasterFromXYZ(id,res=raster::res(env.rast),digits = 10) 
  id_rast<-raster::resample(id_rast,env.rast)
  raster::extent(id_rast)<-raster::extent(env.rast)
  
  PC12<-stack(rpc$map[[c("PC1", "PC2")]],id_rast)
  abio.st<-rast(stack(id_rast, env.rast))
  abio.ex<-na.omit(terra::extract(abio.st, grid.vec, cells=FALSE, df=TRUE))
  PC12ex<-na.omit(terra::extract(rast(PC12), grid.vec, cells=FALSE,df=TRUE))
  PC12pen<-terra::extract(rast(PC12),pen.vec,cells=FALSE,df=TRUE)
 
  PC12ex=PC12ex %>% 
    left_join(PC12pen,by='myID') %>% 
    mutate(PA=ifelse(is.na(ID.y),0,1)) %>% 
    dplyr::select("ID.x", "PC1.x", "PC2.x", "myID", "PA" ) %>% 
    rename(ID="ID.x", PC1="PC1.x" ,PC2=   "PC2.x") %>% 
    drop_na()

  if (is.null(H)) {
    H <- Hpi(x=PC12ex[,c( "PC1", "PC2")])
  } 
  
  estimate <- data.frame(KDE=kde(PC12ex[,c( "PC1", "PC2")], 
                                 eval.points = PC12ex[,c( "PC1", "PC2")],h=H)$estimate, PC12ex[,-1])
  quantP<-quantile(estimate[estimate$PA=='1',1], thres)
  estimate$percP= ifelse(estimate$KDE <= unname(quantP[1]),'out','in')
  
  fullDB.sp=estimate %>% 
    left_join(dt[,c("x", "y", "myID")],by='myID') %>% 
    left_join(abio.ex[,2:ncol(abio.ex)],by='myID') %>% 
    filter(PA == 0 & percP == "out") %>% 
    drop_na() %>% 
    st_as_sf(coords = c("PC1", "PC2"))

  message("\nPerforming background points sampling in the environmental space\n")
  Res <- uesampling(sdf = fullDB.sp, grid.res=grid.res,  n.tr = n.tr, sub.ts = sub.ts, n.ts = n.ts,
                    plot_proc = plot_proc)
 
  if(sub.ts) {
    
    message("\n", paste(nrow(Res$Bkg.tr), "training background points sampled in the environmental space, \n and", nrow(Res$Bkg.ts), "testing background points sampled in the environmental space.",  sep = " "), "\n")  
    message("\n",paste("Estimated final prevalence of", round(nrow(pres)/nrow(Res$Bkg.tr),2), "instead of", prev ), "\n")
    } else {
    
      message("\n",paste(nrow(Res), "training background points sampled in the environmental space", sep = " "), "\n") 
      message("\n",paste("Estimated final prevalence of", round(nrow(pres)/nrow(Res),2), "instead of", prev ),"\n")
  
      }
  
  return(Res)
}
