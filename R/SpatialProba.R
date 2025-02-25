#' Virtual species probability of occurrence
#' 
#' The \code{SpatialProba}  function calculates the simulated probability of occurrence of a virtual species based on an additive model that incorporates environmental variables. The model considers both linear and quadratic relationships between the environmental factors and the species' probability of presence. 
#' This function uses environmental data provided as a SpatRaster object (e.g., temperature, precipitation) to compute the probability of species presence across a defined area of interest. 
#' The resulting probabilities are mapped to a range between 0 and 1, representing the likelihood of species occurrence in the given locations.

#' @param coefs a named vector of regression parameters. Names must match those of the environmental layers (except for intercept, and quadratic terms). Parameters for quadratic terms must have the prefix 'quadr_' (e.g., `quadr_bio1`).
#' @param env.rast a SpatRaster object with environmental layers to generate the spatial layer of probabilities.
#' @param quadr_term a named vector with names of coefs for which a quadratic term is specified (without prefix 'quadr_').
#' @param marginalPlots logical, if TRUE, returns marginal plots.
#' @importFrom grDevices recordPlot
#' @importFrom stats plogis
#' @export
#' @usage
#' SpatialProba(coefs, env.rast, quadr_term, marginalPlots)
SpatialProba <- function(coefs=NULL, env.rast=NULL, quadr_term = NULL, marginalPlots=TRUE) {
  #check if names(coefs) is null
  if(isTRUE(is.null(names(coefs)))) stop("coefs must be a named vector")
  #check if the input env.rast is a SpatRaster
  if (inherits(env.rast, "BasicRaster")) {
    env.rast <- terra::rast(env.rast)
  } 
  #check that if quadr_term is not null a 'quad_' name is in coefs
  if(isTRUE(!is.null(quadr_term))) {
    if(!isTRUE(any(grepl("quad", x = names(coefs))))) stop("quadr term specified, but no term with prefix 'quad' found in coefs")
  } 
  #get names of predictors (excluding intercept and names of quadratic terms, if any)
  if(isTRUE(!is.null(quadr_term))) {
    coefs_names <- names(coefs[!grepl("quad|intercept", x = names(coefs))]) 
  } else {
    coefs_names <- names(coefs[!grepl("intercept", x = names(coefs))]) 
  }
  #check if all coefs_names are in names(env.rast)
  if(isTRUE(!all(coefs_names %in% names(env.rast)))) stop("not all names in coefs are found in env.rast")
  #extract intercept
  mod_intrcpt <- coefs[["intercept"]]
  #subset env.rast if needed
  if(isTRUE(length(coefs_names) < terra::nlyr(env.rast))) env.rast <- env.rast[[coefs_names]]
  #coerce env.rast to a data.frame
  env_df <- terra::as.data.frame(env.rast, na.rm = TRUE, xy = TRUE)
  #get coords
  env_coords <- as.matrix(env_df[c("x", "y")])
  #get rid of coords
  env_df <- env_df[-c(1, 2)]
  #add quadr_term(s) (if any)
  if(!is.null(quadr_term)) {
    message(paste("Adding quadratic term for:", paste(quadr_term, collapse = " ")))
    quad_df <- data.frame(lapply(quadr_term, function(coef_nm) (env_df[[coef_nm]])^2))
    colnames(quad_df) <- paste0("quad_", quadr_term)
    env_df <- cbind(env_df, quad_df)
  }
  #coerce env_df to a matrix
  env_mat <- as.matrix(cbind("intercept" = 1, env_df))
  #re-order cols in env_mat to match coefs names
  env_mat <- env_mat[, names(coefs)]
  #check
  if(isTRUE(!all(colnames(env_mat) == names(coefs)))) stop("Names do not match between environmental matrix and coef vector")
  #compute probabilities (logit scale)
  proba_link <- env_mat%*%unname(coefs)
  #back-transform to proba (response) scale
  proba_resp <- stats::plogis(proba_link)
  #rasterize to get the spatial layer
  spatial_proba <- terra::rasterize(x = env_coords, y = env.rast[[1]], values = proba_resp)
  names(spatial_proba) <- "TrueProba"
  
  #compute marginal effects: assuming only additive models now 
  if(isTRUE(marginalPlots)) {
    plotOut <- list()
    if(isTRUE(!is.null(quadr_term))) {
      #compute marginal effect for the quadratic term
      env_mat.tmp<-env_mat
      #fix the non quadratic predictors taking the mean
      env_mat.tmp[ ,coefs_names[!coefs_names %in% c(quadr_term, paste0("quad_", quadr_term))]] <- mean(env_mat.tmp[ ,coefs_names[!coefs_names %in% c(quadr_term, paste0("quad_", quadr_term))]])
      proba_link <- env_mat.tmp%*%unname(coefs)
      proba_resp <- plogis(proba_link)
      marg.eff<-data.frame(quadr_term=env_mat.tmp[,quadr_term], PA=proba_resp)
      colnames(marg.eff)<-c(quadr_term, "PA")
      marg.eff <- marg.eff[order(marg.eff[, 1]), ] 
      # plot
     plot(x = marg.eff[, 1], y = marg.eff[, 2], 
           xlab=quadr_term, ylab="Probability of Presence", 
           ylim=c(0,1), type="l",  bty = "n")
     p <- grDevices::recordPlot()
      plotOut[[quadr_term]]<- p 
    }
  
    #compute marginal effect for the others term
    env_mat.tmp<-env_mat
    #fix the other predictors taking the mean
    env_mat.tmp <- cbind(env_mat.tmp[ , which(!colnames(env_mat.tmp) %in% c(quadr_term, paste0("quad_", quadr_term)))],  
                     data.frame(t(colMeans(env_mat.tmp[ , c(quadr_term, paste0("quad_", quadr_term))]))))
    #re-order cols in env_mat to match coefs names
    env_mat.tmp <- as.matrix(env_mat.tmp[, names(coefs)])
    proba_link <- env_mat.tmp%*%unname(coefs)
    proba_resp <- plogis(proba_link)
    myname <- names(coefs[which(!colnames(env_mat.tmp) %in% c(quadr_term, paste0("quad_", quadr_term), "intercept"))])
    marg.eff<-data.frame(x=env_mat.tmp[, myname], PA=proba_resp)
    colnames(marg.eff)<-c(myname, "PA")
    marg.eff <- marg.eff[order(marg.eff[, 1]), ] 
    
    # plot
    p <- plot(x = marg.eff[, 1], y = marg.eff[, 2], 
              xlab=myname, ylab="Probability of Presence", 
              ylim=c(0,1), type="l",  bty = "n")
    p <- recordPlot()
    plotOut[[myname]] <- p
    plotOut<-cowplot::plot_grid(plotlist = plotOut, nrow=1, ncol=length(plotOut), labels = "AUTO")
    
    return(list(rast=spatial_proba, margEff = plotOut))

    } else {
    return(spatial_proba)
  }
}


