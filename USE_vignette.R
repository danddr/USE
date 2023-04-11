## ----setup, message=FALSE, echo=TRUE, warning=FALSE, results = FALSE------------------------------------------------------
# Sys.setlocale("LC_ALL", "English")
library(geodata)
library(USE)
library(terra)
library(raster)
library(virtualspecies)
library(sf)
library(ggplot2)


## ---- eval=FALSE, message=FALSE-------------------------------------------------------------------------------------------
## Worldclim <- geodata::worldclim_global(var='bio', res=10, path=getwd())
## envData <- terra::crop(Worldclim, terra::ext(-12, 25, 36, 60))

## ---- eval=TRUE, echo=FALSE, message=FALSE, results=FALSE-----------------------------------------------------------------
envData <- USE::Worldclim_tmp


## ---- eval=FALSE, message=FALSE-------------------------------------------------------------------------------------------
## #create virtual species
## myRandNum <- sample(1:19,size=5, replace = FALSE)
## envData <- envData[[myRandNum]]

## ---- eval=TRUE, message=FALSE, warning=FALSE, fig.height = 8, fig.width = 8, fig.align='center'--------------------------
set.seed(123)
random.sp <- virtualspecies::generateRandomSp(raster::stack(envData), 
                                              convert.to.PA = FALSE, 
                                              species.type = "additive",
                                              realistic.sp = TRUE, 
                                              plot = FALSE)
#reclassify suitability raster using a probability conversion rule
new.pres <- virtualspecies::convertToPA(x=random.sp, 
                      beta=0.55,
                      alpha = -0.05, plot = FALSE)
#Sample true occurrences
presence.points <- virtualspecies::sampleOccurrences(new.pres,
                                     n = 300, # The number of points to sample
                                     type = "presence-absence",
                                     sample.prevalence = 0.6,
                                     detection.probability = 1,
                                     correct.by.suitability = TRUE,
                                     plot = TRUE)  


## ---- eval=TRUE-----------------------------------------------------------------------------------------------------------
myPres <- presence.points$sample.points[which(presence.points$sample.points$Observed==1), c( "x", "y",  "Observed")]
myPres <- st_as_sf(myPres, coords=c("x", "y"), crs=4326)


## ---- eval=TRUE-----------------------------------------------------------------------------------------------------------
rpc <- rastPCA(envData, stand = TRUE)
dt <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))
dt <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))


## ---- eval=FALSE, echo=TRUE-----------------------------------------------------------------------------------------------
## myRes <- USE::optimRes(sdf=dt,
##                     grid.res=c(1:10),
##                     perc.thr = 20,
##                     showOpt = TRUE,
##                     cr=5)


## ---- eval=TRUE, echo=FALSE, message=FALSE, results="hide"----------------------------------------------------------------
myRes <- list()
myRes$Opt_res <- 5

## ---- eval=TRUE-----------------------------------------------------------------------------------------------------------
myRes$Opt_res


## ---- eval=TRUE, message=FALSE--------------------------------------------------------------------------------------------
myObs <- USE::uniformSampling(sdf=dt, 
                              grid.res=myRes$Opt_res,
                              n.tr = 5,
                              sub.ts = TRUE,
                              n.ts = 2,
                              plot_proc = FALSE)


## ---- eval=TRUE-----------------------------------------------------------------------------------------------------------
head(myObs$obs.tr)


## ---- eval=TRUE, message=FALSE, warning=FALSE, fig.height = 6, fig.width = 6, fig.align='center'--------------------------
env_pca <- c(rpc$PCs$PC1, rpc$PCs$PC2)
env_pca <- na.omit(as.data.frame(env_pca))

ggplot(env_pca, aes(x=PC1))+
  geom_density(aes(color="Environment"), size=1 )+
  geom_density(data=data.frame(st_coordinates(myObs$obs.tr)), 
               aes(x=X,  color="Uniform"), size=1)+
  scale_color_manual(name=NULL, 
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60'))+     
  labs(y="Density of PC-scores")+
  xlim(-5,3)+ ylim(0,1)+
  theme_classic()+
  theme(legend.pos="bottom",  
        text = element_text(size=14),  
        legend.text=element_text(size=12))

ggplot(env_pca, aes(x=PC2))+
  geom_density(aes(color="Environment"), size=1 )+
  geom_density(data=data.frame(st_coordinates(myObs$obs.tr)), 
               aes(x=Y,  color="Uniform"), size=1)+
  scale_color_manual(name=NULL, 
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60'))+     
  labs(y="Density of PC-scores")+
  xlim(-5,3)+ ylim(0,1)+
  theme_classic()+
  theme(legend.pos="bottom",  
        text = element_text(size=14),  
        legend.text=element_text(size=12))


## ---- eval=TRUE, message=FALSE--------------------------------------------------------------------------------------------
myGrid.psAbs <- USE::paSampling(env.rast=envData,
                                pres=myPres,
                                thres=0.75,
                                H=NULL,
                                grid.res=as.numeric(myRes$Opt_res),
                                n.tr = 5,
                                prev=0.3,
                                sub.ts=TRUE,
                                n.ts=5,
                                plot_proc=FALSE,
                                verbose=FALSE)


## ---- eval=TRUE, message=FALSE, warning=FALSE, fig.height = 6, fig.width = 6, fig.align='center'--------------------------
ggplot(env_pca, aes(x=PC1))+
  geom_density(aes(color="Environment"), size=1 )+
  geom_density(data=data.frame(st_coordinates(myGrid.psAbs$obs.tr)), 
               aes(x=X,  color="Uniform"), size=1)+
  geom_density(data=terra::extract(c(rpc$PCs$PC1, rpc$PCs$PC2), myPres, df=TRUE), 
               aes(x=PC1, color="Presence"), size=1 )+
  scale_color_manual(name=NULL, 
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60', "Presence"="black"))+     
  labs(y="Density of PC-scores")+
  xlim(-5,3)+ ylim(0,1)+
  theme_classic()+
  theme(legend.pos="bottom",  
        text = element_text(size=14),  
        legend.text=element_text(size=12))

ggplot(env_pca, aes(x=PC2))+
  geom_density(aes(color="Environment"), size=1 )+
  geom_density(data=data.frame(st_coordinates(myGrid.psAbs$obs.tr)), 
               aes(x=Y,  color="Uniform"), size=1)+
  geom_density(data=terra::extract(c(rpc$PCs$PC1, rpc$PCs$PC2), myPres, df=TRUE), 
               aes(x=PC2, color="Presence"), size=1 )+
  scale_color_manual(name=NULL, 
                     values=c('Environment'='#1E88E5', 'Uniform'='#D81B60', "Presence"="black"))+     
  labs(y="Density of PC-scores")+
  xlim(-5,3)+ ylim(0,1)+
  theme_classic()+
  theme(legend.pos="bottom",  
        text = element_text(size=14),  
        legend.text=element_text(size=12))


## ---- eval=TRUE, message=FALSE, warning=FALSE, fig.height = 6, fig.width = 6, fig.align='center'--------------------------
ggplot()+
  stars::geom_stars(data = stars::st_as_stars(new.pres$pa.raster), alpha = 0.5 )+
  scale_fill_manual(values =viridis::viridis(2),  na.value = "transparent", 
                  labels = c("Absence", "Presence", ""))+
  geom_sf(data=myPres, 
          aes(color= "Presences"), alpha=1, size=2, shape= 19)+
  geom_sf(data=st_as_sf(st_drop_geometry(myGrid.psAbs$obs.tr), 
                        coords = c("x","y"), crs=4326), 
          aes(color="Pseudo-absences"), alpha=0.8, size=2, shape = 19 )+
    scale_colour_manual(name=NULL, 
                      values=c('Presences'='steelblue', 'Pseudo-absences'='#A41616'))+
  labs(x="Longitude",y="Latitude", fill="Virtual species")+
  theme_light()+
  theme(legend.pos="bottom",  
        legend.background=element_blank(),
        legend.box="vertical",
        panel.grid = element_blank(),
        text = element_text(size=14),  
        legend.text=element_text(size=14), 
        aspect.ratio = 1, 
        panel.spacing.y = unit(2, "lines"))

