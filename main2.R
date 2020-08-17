library(GWmodel)
library(stars)
library(rgee)
library(raster)
library(sf)
library(sp)

source("utils.R")
ee_Initialize("csaybar", gcs = TRUE)

# 1. Load geometry dataset (see main.R)
dataset <- read_sf("data/dataset.geojson")
grid <- read_sf("data/new_grid_15km.geojson")
#district_selva <- read_sf("data/Distritos/Distritos.shp")
# study_area <- district_selva[["geometry"]] %>%
#   st_buffer(.1) %>%
#   st_union()
# my_grid <- st_as_stars(study_area, dx=0.15, dy=0.15) %>%
#   st_coordinates() %>%
#   st_as_sf(coords = c(1,2)) %>%
#   `st_crs<-`(4326) %>%
#   st_intersection(study_area)
# write_sf(my_grid, "data/new_grid_15km.geojson")

# 4. GWR
dt_sf_p <- st_cast(st_centroid(dataset[["geometry"]]),"POINT")
DM1 <- gw.dist(dp.locat = st_coordinates(dt_sf_p))
DM2 <- gw.dist(dp.locat = st_coordinates(dt_sf_p), rp.locat = st_coordinates(grid))

# 4.1 Calculate kernel bandwidth (0.225)
bw1 <- bw.gwr(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset, "Spatial"),
  kernel = "gaussian",
  dMat = DM,
  adaptive = TRUE,
  longlat = TRUE
)


# 4.2 Calculate collineallity
coll <- gwr.collin.diagno(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset, "Spatial"),
  bw = bw1,
  kernel = "gaussian",
  adaptive = TRUE,
  dMat = DM
)

# 4.3 Model Selection
# model.sel<-gwr.model.selection(
#   DeVar = "target",
#   InDeVars = c("dem", "slp", "tpi", "cpo", "rio", "anp", "acr", "via", "cag", "cat"),
#   data = as(dataset2, "Spatial"),
#   kernel = "gaussian",
#   dMat=DM
# )


# 4.4 GWR calibrartion
gwr.res1 <- ggwr.basic(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset, "Spatial"),
  bw = bw1,
  family = "binomial",
  adaptive = TRUE,
  kernel = "gaussian",
  longlat = TRUE,
  dMat = DM
)

gwr.res2 <- ggwr.basic(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset, "Spatial"),
  regression.points = as(grid, "Spatial"),
  bw = bw1,
  family = "binomial",
  adaptive = TRUE,
  kernel = "gaussian",
  longlat = TRUE,
  dMat = DM2
)

# 5. Save results
results_raster <- gwr.res2$SDF
save_raster <- function(x) {
  r_upload <- results_raster[x]
  r_name <- sprintf("results/%s.tif",names(r_upload))
  writeRaster(rasterFromXYZ(r_upload), r_name, overwrite = TRUE)
}
lapply(1:11, save_raster)

# 6. Upload results (raster)
# ee_manage_create("users/csaybar/forest/coeficients", "ImageCollection")
# ee_manage_delete("users/csaybar/forest/coeficients", "ImageCollection")
upload_results <- function(path, ee_path) {
  predictors <- c("dem", "slp", "tpi", "cpo", "rio", "anp", "acr", "via", "cag", "cat", "Intercept")
  files_to_upload <- sprintf("%s/%s.tif", normalizePath(path), predictors)
  for (index in seq_along(files_to_upload)) {
    file_to_upload <- files_to_upload[index]
    stars_obj <- read_stars(file_to_upload)
    st_crs(stars_obj) <- 4326
    asset_id <- paste0(ee_path,"/", gsub("\\.tif$", "", basename(file_to_upload)))
    raster_as_ee(stars_obj, asset_id, bucket = "rgee_dev", overwrite = TRUE)
  }
}

path <- "/home/csaybar/Documents/Github/gee_gwlr/results/"
ee_path <- paste0(ee_get_assethome(),"/forest/coeficients")
upload_results(path, ee_path)

# 7. Upload results (raster)
ee_predictor <- create_predictors2()$toBands()
coefficients <- ee$ImageCollection("users/csaybar/forest/coeficients")$toBands()
coefficients <- coefficients$select(c(6,8,9,5,7,2,1,10,3,4,0))

# ee_predictor$bandNames()$getInfo()
# coefficients$bandNames()$getInfo()
final_map <- ee_predict(coefficients, ee_predictor)

# Mask to values
mask <- ee$ImageCollection("users/csaybar/forest/PeruForest_bnb2018")$mosaic() %>%
  generate_target() %>%
  "[["("dataset") %>%
  ee$Image$multiply(0) %>%
  ee$Image$add(1) %>%
  ee$Image$unmask(0) %>%
  ee$Image$updateMask(.,.)

#Map$addLayer(final_map$updateMask(mask),list(min=0,max=0.7,palette=c("blue","green","red")))
final_map_raster <- ee_as_raster(
  image = final_map$updateMask(mask),
  via = "drive",
  scale = 1000,
  dsn = "/home/aybarpc01/Github/gee_gwlr/results/final_map.tif"
)
