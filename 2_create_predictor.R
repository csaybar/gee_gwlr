#' This script re-scale and save all the predictors (ee$Image) to create
#' the forest model risk map in a ee$ImageCollection

library(rgee)
ee_Initialize("data.colec.fbf")

source("utils.R")

gge <- sf_as_ee(st_as_sf(study_area()))$geometry()
cost <- ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
  ee$Image$unitScale(0, 3000)

# Population
amz_population <- ee$ImageCollection("WorldPop/GP/100m/pop")$max()$gt(10) %>% ee$Image$clip(gge)
cost_pop <- cost$cumulativeCost(
  source = amz_population,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 1e5)

task <- ee_image_to_asset(
  image = cost_pop,
  scale = 1000,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost",
  assetId = "users/datacolecfbf/c_pop",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()

# ANP
anp <- ee$Image("users/csaybar/forest/pred/anp")
cost_anp <- cost$cumulativeCost(
  source = anp,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 1e5)
Map$addLayer(cost_anp)

task <- ee_image_to_asset(
  image = cost_anp,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost",
  assetId = "users/csaybar/forest/predictors/c_anp",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()

# ACR
acr <- ee$Image("users/csaybar/forest/pred/acr")
cost_acr <- cost$cumulativeCost(
  source = acr,
  maxDistance =  500 * 1000
) %>%
  ee$Image$unitScale(0, 1e5)
Map$addLayer(cost_acr)

task <- ee_image_to_asset(
  image = cost_acr,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost",
  assetId = "users/csaybar/forest/predictors/c_acr",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()

#Agriculture
cag <- ee$Image("users/csaybar/forest/pred/cag")
cost_cag <- cost$cumulativeCost(
  source = cag,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 5e4)

task <- ee_image_to_asset(
  image = cost_cag,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost",
  assetId = "users/csaybar/forest/predictors/c_cag",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()

# Water
bnb_2018 <- ee$ImageCollection("users/csaybar/forest/PeruForest_bnb2018")
water <- generate_target(bnb_2018$mosaic())$water
water <- water$updateMask(water)

cost_water <- cost$cumulativeCost(
  source = water,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 1e4)

task <- ee_image_to_asset(
  image = cost_water,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost water",
  assetId = "users/csaybar/forest/predictors/c_water",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()

# Via 1
vi1 <- ee$Image("users/csaybar/forest/pred/vi1")
cost_vi1 <- cost$cumulativeCost(
  source = vi1,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 1e5)

task <- ee_image_to_asset(
  image = cost_vi1,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost via1",
  assetId = "users/csaybar/forest/predictors/c_vi1",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()


# Via 2
vi2 <- ee$Image("users/csaybar/forest/pred/vi2")
# Map$addLayer(vi2,list(min=0,max=0))
cost_vi2 <- cost$cumulativeCost(
  source = vi2,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 5e4)

task <- ee_image_to_asset(
  image = cost_vi2,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost via2",
  assetId = "users/csaybar/forest/predictors/c_vi2",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()


# Mining
cat <- ee$Image("users/csaybar/forest/pred/cat")
cost_cat <- cost$cumulativeCost(
  source = cat,
  maxDistance =  300 * 1000
) %>%
  ee$Image$unitScale(0, 5e4)

task <- ee_image_to_asset(
  image = cost_cat,
  region = sf_as_ee(st_as_sf(study_area()))$geometry(),
  description = "uploading cost cat",
  assetId = "users/csaybar/forest/predictors/cat",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()





# Elevation
ee_elevation <- ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
  ee$Image$unitScale(0, 3000) %>%
  ee$Image$rename("ele")
# Map$addLayer(
#   eeObject = ee_elevation_scale$clip(sf_as_ee(study_area())$geometry()),
#   visParams = list(min=0,max=1)
# )
task <- ee_image_to_asset(
  image = ee_elevation,
  region = sf_as_ee(study_area())$geometry(),
  description = "uploading ele",
  assetId = "users/csaybar/forest/predictors/ele",
  maxPixels = 10**13,
  overwrite = TRUE
)
task$start()
ee_monitoring()

# Slope
ee_slope <- ee$Image("WWF/HydroSHEDS/03CONDEM") %>%
  ee$Terrain$slope() %>%
  ee$Image$divide(180) %>%
  ee$Image$multiply(pi) %>%  # in radians
  ee$Image$tan() %>%  # in percentage <0-1>
  ee$Image$multiply(100) %>%
  ee$Image$unitScale(0, 45) %>%
  ee$Image$rename("slp")

task <- ee_image_to_asset(
  image = ee_slope,
  region = sf_as_ee(study_area())$geometry(),
  description = "uploading slp",
  assetId = "users/csaybar/forest/predictors/slp",
  maxPixels = 10**13,
  overwrite = TRUE
)

# task$start()
# ee_monitoring()


# The mTPI distinguishes ridge from valley forms
ee_mtpi <- ee$Image('CSP/ERGo/1_0/Global/ALOS_mTPI') %>%
  ee$Image$unitScale(-100, 100) %>%
  ee$Image$rename("tpi")

task <- ee_image_to_asset(
  image = ee_mtpi,
  region = sf_as_ee(study_area())$geometry(),
  description = "uploading tpi",
  assetId = "users/csaybar/forest/predictors/tpi",
  maxPixels = 10**13,
  overwrite = TRUE
)
# task$start()


dsadas <- cost_cat$projection()
dsadas$getInfo()
