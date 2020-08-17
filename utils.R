#' Sigmoid in Earth Engine
ee_sigmoid <- function(z) {
  ee_one <-  ee$Image(1)
  ee_z <- ee$Image(exp(1))$pow(z$multiply(-1))
  ee_one$divide(ee_one$add(ee_z))
}

#' Logistic Regression Predict
ee_predict <- function(coefficients, ee_predictor) {
  # Name of predictors
  bandnames <- ee_predictor$bandNames()$getInfo()
  coefficients_names <- coefficients$bandNames()$getInfo()
  # Coefficients
  slope <- coefficients$select(0:(length(coefficients_names) - 2))
  intercept <- coefficients$select(length(coefficients_names) -1)

  if (!(length(coefficients_names) - 1) == length(bandnames)) {
    stop("Diferent number of predictors!")
  }

  # # From Local to Earth Engine
  # ee_coef <- do.call(ee$Image$cat, lapply(slope, function(x) ee$Image(x)))

  # Linear Regression and apply sigmoid
  range_coef <- seq_len(length(bandnames)) - 1
  ee_img_mul <- function(x) slope$select(x)$multiply(ee_predictor$select(x))
  ee_lr <- do.call(
    what = ee$Image$cat,
    args = lapply(range_coef, ee_img_mul)
  )
  ee_sigmoid(ee_lr$reduce(ee$Reducer$sum())$add(ee$Image(intercept)))
}

# create_points_grid <- function(x, points = 5000) {
#   read_stars(x, proxy = TRUE) %>%
#     st_bbox() %>%
#     st_as_sfc() %>%
#     st_transform(4326) ->
#     roi
#   potential_points <- st_sample(roi,size = points)
#   return(potential_points)
# }


generate_target <- function(bnb_2018) {
  other_classes <- bnb_2018$eq(2)
  forest <- bnb_2018$eq(3)$multiply(1)
  water <- bnb_2018$eq(4)
  target <- bnb_2018$eq(5)$multiply(2)

  dataset <- target$add(forest)
  dataset <- dataset$updateMask(dataset)$add(-1)
  list(dataset = dataset, water = water, other = other_classes)
}


study_area <- function() {
  area <- c(-79.980469, -15.004206, -68.071289,   0.129071 )
  class(area) <- "bbox"
  roi <- st_as_sfc(area)
  st_crs(roi) <- "EPSG:4326"
  roi
}

img_to_vector <- function(img) {
  img$updateMask(img)$reduceToVectors(
    geometry = sf_as_ee(study_area()),
    crs = "EPSG:4326",
    scale = 1000,
    bestEffort = TRUE,
    geometryType = "polygon",
    eightConnected = FALSE,
    reducer = NULL
  )
}

image_area <- function(img) {
  img$multiply(ee$Image$pixelArea())$reduceRegion(
    reducer = ee$Reducer$sum(),
    geometry = sf_as_ee(study_area()),
    scale = 250,
    maxPixels = 10**13
  )
}

split_extract <- function(x, y, by = 1000) {
  y_len <- y$size()$getInfo()
  for (r_index in seq(1, y_len, by)) {
    index <- r_index - 1
    print(sprintf("Extracting information [%s/%s] ...",index, y_len))
    ee_potential_points <- ee$FeatureCollection(y) %>%
      ee$FeatureCollection$toList(by,index) %>%
      ee$FeatureCollection()
    if (r_index == 1) {
      dataset <- ee_extract(
        x = x,
        y = ee_potential_points,
        sf = FALSE
      )
    } else {
      db_local <- ee_extract(x = x, y = ee_potential_points, sf = FALSE)
      dataset <- rbind(dataset, db_local)
    }
  }
  return(dataset)
}

create_predictors2 <- function() {
  dem <- ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
    ee$Image$unitScale(0, 2500) %>%
    ee$Image$rename("dem")

  slope <- ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
    ee$Terrain$slope() %>%
    ee$Image$divide(180) %>%
    ee$Image$multiply(pi) %>%  # in radians
    ee$Image$tan() %>%  # in percentage <0-1>
    ee$Image$multiply(100) %>%
    ee$Image$expression(
      paste0(
        "( b('slope') > 45) ? 3",
        ": (b('slope') > 25) ? 2",
        ": (b('slope') > 15) ? 1",
        ": 0"
      )
    ) %>%
    ee$Image$rename("slp")

  tpi <- ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI") %>%
    ee$Image$unitScale(-100, 100) %>%
    ee$Image$rename("tpi")

  cpo <- ee$Image("users/csaybar/forest/predictors/cp") %>%
    ee$Image$unitScale(0, 1000) %>%
    ee$Image$rename("cpo")

  rio <- ee$Image("users/csaybar/forest/predictors/rios") %>%
    ee$Image$unitScale(0, 1000) %>%
    ee$Image$rename("rio")

  anp <- ee$Image("users/csaybar/forest/predictors/anp") %>%
    ee$Image$unitScale(0, 2000) %>%
    ee$Image$rename("anp")

  acr <- ee$Image("users/csaybar/forest/predictors/acr") %>%
    ee$Image$unitScale(0, 2000) %>%
    ee$Image$rename("acr")

  via <- ee$Image("users/csaybar/forest/predictors/via") %>%
    ee$Image$unitScale(0, 1000) %>%
    ee$Image$rename("via")

  cag <- ee$Image("users/csaybar/forest/predictors/cag") %>%
    ee$Image$unitScale(0, 1000) %>%
    ee$Image$rename("cag")

  cat <- ee$Image("users/csaybar/forest/predictors/ctm") %>%
    ee$Image$unitScale(0, 1000) %>%
    ee$Image$rename("cat")

  # nat <- ee$Image("users/csaybar/forest/predictors/nat") %>%
  #   ee$Image$unitScale(0, 1000) %>%
  #   ee$Image$rename("nat")

  ee$ImageCollection(list(dem, slope, tpi, cpo, rio, anp, acr,
                          via, cag, cat))
}

create_predictors <- function(roi, ic, resolution = 100, monitoring = TRUE) {
  try(ee_manage_delete(ic),silent = TRUE)
  ee_manage_create(ic, asset_type = "ImageCollection")

  ##### Population -----------------------------------------
  ee$ImageCollection("WorldPop/GP/100m/pop") %>%
    ee$ImageCollection$max() %>%
    ee$Image$gt(10) %>%
    ee$Image$clip(roi) %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost population",
      assetId = sprintf("%s/c_pop",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading cost population")
  if (monitoring) {
    ee_monitoring()
  }


  ##### Protected natural areas  ------------------------------------
  ee$Image("users/csaybar/forest/pred/anp") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost ANP",
      assetId = sprintf("%s/c_anp",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading ANP")
  if (monitoring) {
    ee_monitoring()
  }


  #####  Regional Conservation Area ----------------------------------
  ee$Image("users/csaybar/forest/pred/acr") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  500 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost ACR",
      assetId = sprintf("%s/c_acr",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading ACR")
  if (monitoring) {
    ee_monitoring()
  }


  # Agriculture ------------------------------------------------------
  ee$Image("users/csaybar/forest/pred/cag") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost CAG",
      assetId = sprintf("%s/c_cag",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading CAG")
  if (monitoring) {
    ee_monitoring()
  }

  # Water ------------------------------------------------------------
    ee$Image("users/csaybar/forest/pred/prs") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost WAT",
      assetId = sprintf("%s/c_wat",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading WAT")
  if (monitoring) {
    ee_monitoring()
  }


  # Via 1 ------------------------------------------------------------
  ee$Image("users/csaybar/forest/pred/vi1") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost VI1",
      assetId = sprintf("%s/c_vi1",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading VI1")
  if (monitoring) {
    ee_monitoring()
  }


  # Via 2 ------------------------------------------------------------
  ee$Image("users/csaybar/forest/pred/vi2") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost VI2",
      assetId = sprintf("%s/c_vi2",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading VI2")
  if (monitoring) {
    ee_monitoring()
  }


  # Mining -----------------------------------------------------------
  ee$Image("users/csaybar/forest/pred/cat") %>%
    ee$Image$cumulativeCost(
      cost = cost,
      source = .,
      maxDistance =  300 * 1000
    ) %>%
    ee$Image$unitScale(0, 1e4) %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading cost CAT",
      assetId = sprintf("%s/c_cat",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading CAT")
  if (monitoring) {
    ee_monitoring()
  }


  # Elevation -----------------------------------------------------------
  ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
    ee$Image$unitScale(0, 3000) %>%
    ee$Image$rename("ele") %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading ELE",
      assetId = sprintf("%s/ele",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading ELE")
  if (monitoring) {
    ee_monitoring()
  }


  # Slope ----------------------------------------------------------------
  ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
    ee$Terrain$slope() %>%
    ee$Image$divide(180) %>%
    ee$Image$multiply(pi) %>%  # in radians
    ee$Image$tan() %>%  # in percentage <0-1>
    ee$Image$multiply(100) %>%
    ee$Image$unitScale(0, 45) %>%
    ee$Image$rename("slp") %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading SLP",
      assetId = sprintf("%s/slp",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading SLP")
  if (monitoring) {
    ee_monitoring()
  }


  # TPI ----------------------------------------------------------------
  ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI") %>%
    ee$Image$unitScale(-100, 100) %>%
    ee$Image$rename("tpi") %>%
    ee_image_to_asset(
      scale = resolution,
      region = roi,
      description = "uploading TPI",
      assetId = sprintf("%s/tpi",ic),
      maxPixels = 10**13,
      overwrite = TRUE
    ) %>%
    ee$batch$Task$start()

  message("Uploading TPI")
  if (monitoring) {
    ee_monitoring()
  }
}

## gaussian: wgt = exp(-.5*(vdist/bw)^2);
near_select <- function(distance, n, bw = bw) {
  distance_v <- as.numeric(distance)
  names(distance_v) <- 1:length(distance_v)
  list(position = as.numeric(names(sort(distance_v)[1:n])),
       value = exp(-.5*(sort(distance_v)[1:n]/bw)^2))
}


GWR_basic <- function(dataset, grid_data, bw = 15000, balance_ratio = 0.5, npoints = 200) {
  target_one <- npoints*balance_ratio #lazzy TODO random sample
  target_zero <- npoints - target_one

  dataset_one <- filter(dataset, target == 1)
  dataset_zero <- filter(dataset, target == 0)


  save_results <- list()
  #save_results <- foreach(nr = 1:100) %dopar%  {
  for (nr in seq_len(nrow(my_grid))) {
    if (nr%%10 == 0) {
      print(nr)
    }

    # one dataset
    ## gaussian: wgt = exp(-.5*(vdist/bw)^2);
    one_d <- sf::st_distance(dataset_one, grid_data[nr, ]$geometry)
    selected_points_one <- near_select(one_d, target_one, bw = bw)
    zero_d <- st_distance(dataset_zero, grid_data[nr, ]$geometry)
    selected_points_zero <- near_select(zero_d, target_zero, bw = bw)

    fdb <- rbind(
      dataset_one[selected_points_one$position, ],
      dataset_zero[selected_points_zero$position, ]
    )

    # First Statistic
    # spconsistency_index <- max(
    #   mean(selected_points_one$value)/mean(selected_points_zero$value),
    #   mean(selected_points_zero$value)/mean(selected_points_one$value)
    # )

    fdb <- na.omit(fdb)
    unq_val <- unique(fdb$slp)
    if (length(unq_val) == 1) {
      if (unq_val == 0) {
        fdb$slp[sample(200,1)] <- 1
      } else if (unq_val == 1) {
        fdb$slp[sample(200,1)] <- 0
      } else if (unq_val == 2) {
        fdb$slp[sample(200,1)] <- 0
      } else if (unq_val == 3) {
        fdb$slp[sample(200,1)] <- 0
      }
    }
    # Main model

    main_model <- suppressWarnings(glm(
      formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
      data = fdb,
      family = binomial(),
      weights = c(selected_points_one$value, selected_points_zero$value)
    ))

    prediction <- as.numeric(predict(main_model,fdb) > 0.5)
    accuracy <- sum(prediction == fdb$target)/length(fdb$target)

    # Contrasted model
    ref_model <-  glm(target ~ 1, data = fdb, family = "binomial")

    # Pvalue
    pvalue <- anova(ref_model, main_model, test = 'LRT')$'Pr(>Chi)'[2]

    model_params <- main_model$coefficients
    intercept <- model_params[1]
    slope <- model_params[2:length(model_params)]
    db_drop <- suppressWarnings(drop1(main_model, test = "LRT"))
    save_results[[nr]] <- list(
      accuracy = accuracy,
      pvalue = pvalue,
      model_coef = list(intercept = intercept, slope = slope),
      LRT = db_drop$LRT[2:length(db_drop$LRT)]
    )
  }
  save_results
}


sf_results <- function(model_results, geom_sfc) {
  geom_sfc <- st_as_sf(geom_sfc)
  geom_sfc['accuracy'] <- sapply(model_results, function(x) x$accuracy)
  geom_sfc['pvalue'] <- sapply(model_results, function(x) x$pvalue)
  geom_sfc['spconsistency_index'] <- sapply(model_results, function(x) x$spconsistency_index)

  # Model
  geom_sfc['intercept'] <- sapply(model_results, function(x) x$model_coef$intercept)
  geom_sfc['slope_acr'] <- sapply(model_results, function(x) x$model_coef$slope[1])
  geom_sfc['slope_anp'] <- sapply(model_results, function(x) x$model_coef$slope[2])
  geom_sfc['slope_cag'] <- sapply(model_results, function(x) x$model_coef$slope[3])
  geom_sfc['slope_cat'] <- sapply(model_results, function(x) x$model_coef$slope[4])
  geom_sfc['slope_cpo'] <- sapply(model_results, function(x) x$model_coef$slope[5])
  geom_sfc['slope_via'] <- sapply(model_results, function(x) x$model_coef$slope[6])
  geom_sfc['slope_rio'] <- sapply(model_results, function(x) x$model_coef$slope[7])
  geom_sfc['slope_dem'] <- sapply(model_results, function(x) x$model_coef$slope[8])
  geom_sfc['slope_slp'] <- sapply(model_results, function(x) x$model_coef$slope[9])
  geom_sfc['slope_tpi'] <- sapply(model_results, function(x) x$model_coef$slope[10])

  # LRT
  geom_sfc['lrt_acr'] <- sapply(model_results, function(x) x$LRT[1])
  geom_sfc['lrt_anp'] <- sapply(model_results, function(x) x$LRT[2])
  geom_sfc['lrt_cag'] <- sapply(model_results, function(x) x$LRT[3])
  geom_sfc['lrt_cat'] <- sapply(model_results, function(x) x$LRT[4])
  geom_sfc['lrt_cpo'] <- sapply(model_results, function(x) x$LRT[5])
  geom_sfc['lrt_via'] <- sapply(model_results, function(x) x$LRT[6])
  geom_sfc['lrt_rio'] <- sapply(model_results, function(x) x$LRT[7])
  geom_sfc['lrt_dem'] <- sapply(model_results, function(x) x$LRT[8])
  geom_sfc['lrt_slp'] <- sapply(model_results, function(x) x$LRT[9])
  geom_sfc['lrt_tpi'] <- sapply(model_results, function(x) x$LRT[10])
  geom_sfc
}

ee_spatial <- function(sf_var, asset_id, bucket, output = "results/", upload = TRUE) {
  sf_var <- na.omit(sf_var)
  to_sp <- st_sf(var = sf_var[[1]], geometry = st_geometry(sf_var))
  to_sp$var[to_sp$var > 100] = 100
  to_sp$var[to_sp$var < -100] = -100
  tmp_file <- tempfile()
  my_file <- sprintf("%s.shp", tmp_file)
  write_sf(to_sp, my_file)

  sp_var <- as(to_sp,"Spatial")
  slope_variogram <- autofitVariogram(var~1, sp_var)
  png(sprintf("%s/%s.png",output, names(sf_var[1])[1]))
  plot(slope_variogram)
  dev.off()
  # results <- ee$FeatureCollection("users/csaybar/forest/forest_risk_results")
  # results$kriging(
  #   propertyName = "accurcy",
  #   shape = "spherical",
  #   range = slope_variogram$var_model$range[2],
  #   sill = slope_variogram$var_model$psill[2],
  #   nugget = slope_variogram$var_model$psill[1],
  #   maxDistance = 100*1000,
  #   reducer = "mean"
  # ) -> var_raster
  # ssss <- ee_image_to_asset(
  #   image = var_raster,
  #   region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  #   scale = 5000,
  #   assetId = "users/csaybar/demo"
  # )
  system(
    sprintf(
      "gdal_grid -l %s -zfield var -a invdist:power=2.0:smothing=1.5:radius1=0.0:radius2=0.0:angle=0.0:max_points=100:min_points=0:nodata=0.0 -ot Float32 -of GTiff %s %s",
      basename(tmp_file),
      my_file,
      sprintf("%s/%s.tif",normalizePath(output), names(sf_var[1])[1])
    )
  )
  if (upload) {
    stars_obj <- stars::read_stars(sprintf("%s/%s.tif",normalizePath(output), names(sf_var[1])[1]))
    raster_as_ee(stars_obj, asset_id, bucket = bucket, overwrite = TRUE)
  }
}


best_predictor <- function(sf_model_results, output = "results/",res = 0.005) {
  ltr_var <- sf_model_results[16:25]
  ltr_var$x <- NULL
  best_predictor <- apply(ltr_var, 1, which.max)
  bp_x <- st_sf(best_p = best_predictor, geometry = sf_model_results$x)
  factor_p <- factor(best_predictor)
  levels(factor_p) <- colnames(ltr_var)
  bp_x$best_pname <- as.character(factor_p)
  # write_sf(bp_x,"data/best_predictor.geojson", delete_dsn = TRUE)
  # #mapview::mapview(bp_x, zcol = "best_p")
  # gdal_grid(
  #   src_datasource = normalizePath("data/best_predictor.geojson"),
  #   dst_filename = sprintf("%s/best_predictor.tif",normalizePath(output)),
  #   a = "nearest:radius1=0.0:radius2=0.0:angle=0.0:nodata=0.0",
  #   txe = c(-79.541512, -68.711505),
  #   tye = c(-14.472808, -0.064539),
  #   outsize = c(1000, 1000),
  #   of="GTiff",
  #   ot="Float64",
  #   l="best_predictor",
  #   zfield = "best_p"
  #)
}


# upload_coeficients<- function(output = "results/") {
#   coef_id <- "users/csaybar/forest/coeficients"
#   try(ee_manage_delete("users/csaybar/forest/coeficients"))
#   ee_manage_create("users/csaybar/forest/coeficients","ImageCollection")
#   raster_files <- list.files(output,"slope_.*\\.tif$",full.names = TRUE)
#   input_file_order <- c(1,2,3,4,6,9,10,11,5,7,8)
#
#   dddd <- raster(raster_files[input_file_order][1])
#   ddddd <- flip(dddd,direction = "y")
#   writeRaster(ddddd, "results/demo.tif",overwrite=TRUE)
#
#   raster_as_ee(
#     x = raster("results/demo.tif"),
#     assetId = sprintf("%s/%s",coef_id,"demo2"),
#     bucket = "rgee_dev"
#   )
# }
