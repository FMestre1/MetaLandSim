range_raster <-
  function(presences.map, re.out, mask.map=NULL, plot=TRUE)
  {
    
    mask_raster <- terra::rast(mask.map)
    presences_raster <- terra::rast(presences.map)
    terra::NAflag(mask_raster) <- 0
    
    #F1
    fit.sigmoid <- function(y, x, start.params = list(a = 1, b = 0.5, c = 0))
    {
      fitmodel <- minpack.lm::nlsLM(y ~ a / (1 + exp(b * (x - c))), start = start.params, control = nls.lm.control(maxiter = 1000))
      return(coef(fitmodel))
    }
    #F2
    predict.sigmoid <- function(params, x)
    {
      return(params[1] / (1 + exp(params[2] * (x - params[3]))))
    }
    
    params.fit <- fit.sigmoid(re.out$PROPORTION, re.out$DISTANCE/1000)
    
    summed_up <- terra::mosaic(presences_raster, mask_raster, fun="sum")
    cost.dist <- terra::costDist(summed_up, target = 2) #distance to cells with value 2
    cost.dist <- cost.dist/1000 

    out_raster <- params.fit[1]/(1 + exp(params.fit[2] * (cost.dist - params.fit[3])))

    if (plot == TRUE) { 
    par(mfrow=c(1,2))
    plot(re.out$DISTANCE,
         predict.sigmoid(params.fit, re.out$DISTANCE/1000),
         type = "l", 
         col = "blue",
         lwd = 2,
         main="Dispersal probability", 
         xlab="Distance (m)", 
         ylab="Probability"
         )
    terra::image(out_raster, 
                 main="Spatial projection"
                 )
    terra::contour(out_raster, add=TRUE)
    }
    
    model_formula <- "probability ~ a / (1 + exp(b * (distance - c)))"
    out1 <- list(model_formula, params.fit, out_raster)
    names(out1) <- c("MODEL", "MODEL_PARAMETERS", "RASTER_PROJECTIONS")
    message("Writing dispersal model rasters.")
    terra::writeRaster(out_raster, filename = "PROB.tif", overwrite=TRUE)
    
    return(out1)
  }
