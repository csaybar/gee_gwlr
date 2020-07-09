library(automap)
library(mapview)
library(dplyr)
library(rgee)
library(stars)
library(gstat)
library(sf)

source("utils.R")

ee_Initialize("csaybar", gcs = TRUE)


# 1. Load target dataset
target_dataset <- read_sf("data/train_dataset.geojson")
# sf_as_ee(
#   x = target_dataset,
#   via = "gcs_to_asset",
#   bucket = "rgee_dev",
#   assetId = "users/csaybar/forest/forest_risk_points",
#   overwrite = TRUE,
#   monitoring = FALSE
# )


# 2. Create Predictors
x <- create_predictors2()$toBands()


# 3. Extract values
y <- ee$FeatureCollection("users/csaybar/forest/forest_risk_points")
dataset <- split_extract(x, y, by = 1000)
names(dataset) <- c("id", "target", "dem", "slp", "tpi", "cpo",
                    "rio", "anp", "acr", "via", "cag", "cat", "geometry")
# write_sf(dataset, "data/dataset.geojson", delete_dsn =TRUE)

# 3.1 create a grid
amazon_roi <- read_sf("data/Selva_reg.geojson") %>%
  st_buffer(0.5) %>%
  '['("geometry")

my_grid <- st_as_stars(study_area(), dx=0.05, dy=0.05) %>%
  st_coordinates() %>%
  st_as_sf(coords = c(1,2)) %>%
  `st_crs<-`(4326) %>%
  st_intersection(amazon_roi)

# 4. Train models
model_results <- GWR_basic(
  dataset = dataset,
  grid_data = my_grid,
  bw = 111111,
  balance_ratio = 0.5,
  npoints = 200
)
# save(model_results,file = "results/model_results.Rdata")

# 5. save results
sf_model_results <- sf_results(model_results, my_grid$geometry)
results_raster <- as(sf_model_results,"Spatial")
save_raster <- function(x) {
  r_upload <- results_raster[x]
  r_name <- sprintf("results/%s.tif",names(r_upload))
  writeRaster(rasterFromXYZ(r_upload), r_name)
}
lapply(1:13,save_raster)


# 6. Upload results (raster)
# ee_manage_create("users/csaybar/forest/coeficients", "ImageCollection")
# ee_manage_delete("users/csaybar/forest/coeficients", "ImageCollection")
upload_results <- function() {
  files_to_upload <- list.files("results/","slope|intercept",full.names = TRUE)
  for (index in seq_along(files_to_upload)) {
    file_to_upload <- files_to_upload[index]
    stars_obj <- raster(file_to_upload)
    crs(stars_obj) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    var_name <- gsub("slope_|\\.tif$","",basename(file_to_upload))
    asset_id <- sprintf("users/csaybar/forest/coeficients/%s", var_name)
    message("Uploading ... ", var_name)
    raster_as_ee(st_as_stars(stars_obj), asset_id, bucket = "rgee_dev", overwrite = TRUE)
  }
}
upload_results()


# 8. Final Map
ee_predictor <- create_predictors2()$toBands()
coefficients <- ee$ImageCollection("users/csaybar/forest/coeficients")$toBands()
coefficients <- coefficients$select(c(5,8,9,4,7,1,0,10,2,3,6))
# ee_predictor$bandNames()$getInfo()
# coefficients$bandNames()$getInfo()
final_map <- ee_predict(coefficients, ee_predictor)
# ee_task <- ee_image_to_asset(
#   image = final_map,
#   assetId = "users/csaybar/forest/finalmap",
#   overwrite = FALSE,
#   scale = 100
# )
# ee_task$start()

mask <- ee$ImageCollection("users/csaybar/forest/PeruForest_bnb2018")$mosaic() %>%
  generate_target() %>%
  "[["("dataset") %>%
  ee$Image$multiply(0) %>%
  ee$Image$add(1) %>%
  ee$Image$unmask(0) %>%
  ee$Image$updateMask(.,.)

final_map_raster <- ee_as_raster(final_map2, via = "drive", scale = 1000,dsn = "/home/aybarpc01/Github/gee_gwlr/results/final_map.tif")

Map$addLayer(final_map2,list(min=0, max=1))
