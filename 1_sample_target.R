#' This script creates ee_potential_points. It is a ee$FeatureCollection
#' with 10 000 points containing information about forest (0) and
#' deforestation (1).
library(rgee)

source("utils.R")
ee_Initialize("csaybar@gmail.com", gcs = TRUE)


bnb_2018 <- ee$ImageCollection("users/csaybar/forest/PeruForest_bnb2018")
point_number <- 12000

# total area
total_area <- generate_target(bnb_2018$mosaic())$dataset
area_total <- image_area(total_area$multiply(0)$add(1))$getInfo()$b1


# points by tile
ee_seed <- 1000
bnb_2018 <- ee$ImageCollection("users/csaybar/forest/PeruForest_bnb2018")$toList(20)

for (r_index in seq_len(20)) {
  index <- r_index - 1
  print(sprintf("Title [%s/20]: Obtaining points ...", r_index))
  # <0-1> forest vs deforestation raster
  images <- generate_target(ee$Image(bnb_2018$get(index)))$dataset
  # Image area
  img_area <- image_area(images$multiply(0)$add(1))$getInfo()$b1
  # Number of points by tile
  tile_points <- ceiling(img_area/area_total*point_number)
  # Stratified Sample
  points <- images$stratifiedSample(
    numPoints = tile_points,
    scale = 30,
    projection="EPSG:4326",
    seed = ee_seed,
    geometries = TRUE
  )
  # Save results in local
  if (r_index == 1) {
    db_points <- ee_as_sf(points, maxFeatures = 10000)["b1"]
  } else {
    amaz_points <- ee_as_sf(points, maxFeatures = 10000)
    if (nrow(amaz_points) == 0) {
      next
    }
    db_points <- rbind(db_points, amaz_points["b1"])
  }
}

potential_points <- db_points[sample(nrow(db_points),10000),]
#plot(potential_points)
#table(potential_points$b1)

names(potential_points) <- c("target", "geometry")
ee_potential_points <- sf_as_ee(
  x = potential_points,
  via = "gcs_to_asset",
  assetId = "users/csaybar/forest/forest_risk_points",
  overwrite = TRUE,
  bucket = "rgee_dev"
)

# ee_monitoring()
# Map$addLayer(ee_potential_points)
