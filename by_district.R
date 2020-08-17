#' GWLR by district
#' @csaybar

library(automap)
library(mapview)
library(spgwr)
library(dplyr)
library(rgee)
library(stars)
library(gstat)
library(sf)
source("utils.R")

ee_Initialize("csaybar", gcs = TRUE)

# 1. Load geometries
district <- read_sf("data/Distritos/Distritos.shp")
y <- ee$FeatureCollection("users/csaybar/distritos_selva")

# 2. Load predictors
bnb_2018 <- ee$ImageCollection("users/csaybar/forest/PeruForest_bnb2018")
ee_target <- generate_target(bnb_2018$mosaic())$dataset
x <- ee_target$addBands(create_predictors2()$toBands())

# 3. Extract values
dataset <- split_extract(x$unmask(), y, by = 100)
dataset2 <- st_sf(dataset, geometry = district$geometry)
colnames(dataset2) <- c("target", "dem", "slp", "tpi", "cpo", "rio", "anp", "acr", "via", "cag", "cat", "geometry")

# 4. GWR
# 4.1 Calculate kernel bandwidth
dt_sf_p <- st_cast(st_centroid(dataset2[["geometry"]]),"POINT")
DM <- gw.dist(dp.locat=st_coordinates(dt_sf_p))
bw1 <- bw.gwr(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset2, "Spatial"),
  kernel = "gaussian",
  dMat = DM
)

# 4.2 Calculate collineallity
coll <- gwr.collin.diagno(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset2, "Spatial"),
  bw = bw1,
  kernel = "gaussian",
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
gwr.res1 <- gwr.basic(
  formula = target ~ acr + anp + cag + cat + cpo + via + rio + dem + slp + tpi,
  data = as(dataset2, "Spatial"),
  bw = bw1,
  kernel = "gaussian",
  longlat = TRUE,
  dMat = DM
)
