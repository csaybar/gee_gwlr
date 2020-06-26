library(mlbench)     # for PimaIndiansDiabetes2 dataset
library(dplyr)       # for data manipulation (dplyr)
library(broom)       # for making model summary tidy
library(visreg)      # for potting logodds and probability
library(margins)     # to calculate Average Marginal Effects
library(rcompanion)  # to calculate pseudo R2
library(ROCR)        # to compute and plot Reciever Opering Curve
library(rgee)        # Google Earth Engine for R
library(stars)       # datacubes in R
library(sf)          # simple features in R

source("utils.R")

ee_Initialize()

x <-"/home/aybarpc01/Downloads/wetransfer-0be8b3/Bosque_No_Bosque_2018_raster/Bosque_No_Bosque_2018.tif"
points <- create_points_grid(x)


