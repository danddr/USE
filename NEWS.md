# USE 0.1.1
* Adding thresh.inspect function to investigate the effect of the kernel threshold on the partition of the environmental space. 

# USE 0.1.0
* Changing the filter on the presences applied to the kernel density: in USE v0.0.0.9000, the kernel density was computed using all the pixels in the whole environmental space and then a filter was applied to the PC-scores associated with the presences. In the new version, the kernel density is computed only in the environmental space associated with the presence observations. 
* Adding the rastPCA function to compute a PCA on SpatRaster and Raster objects. 
* Dependency from RSToolbox removed

# USE 0.0.0.9000
* First release 
