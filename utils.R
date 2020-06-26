#' Sigmoid in Earth Engine
ee_sigmoid <- function(z) {
  ee_one <-  ee$Image(1)
  ee_z <- ee$Image(exp(1))$pow(z$multiply(-1))
  ee_one$divide(ee_one$add(ee_z))
}

#' Logistic Regression Predict
logistic_regression <- function(rmodel, ee_predictor) {
  # Name of predictors
  bandnames <- ee_predictor$bandNames()$getInfo()

  # Coefficients
  coef <- rmodel$coefficients
  names(coef) <- NULL
  slope <- coef[2:length(coef)]
  intercept <- coef[1]

  if (!(length(rmodel$coefficients) - 1) == length(bandnames)) {
    stop("Diferent number of predictors!")
  }

  # From Local to Earth Engine
  ee_coef <- do.call(ee$Image$cat, lapply(slope, function(x) ee$Image(x)))

  # Linear Regression and apply sigmoid
  range_coef <- seq_len(length(slope)) - 1
  ee_img_mul <- function(x) ee_coef$select(x)$multiply(ee_predictor$select(x))
  ee_lr <- do.call(
    what = ee$Image$cat,
    args = lapply(range_coef, ee_img_mul)
  )
  ee_sigmoid(ee_lr$reduce(ee$Reducer$sum())$add(ee$Image(intercept)))
}


create_points_grid <- function(x, points = 5000) {
  read_stars(x, proxy = TRUE) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_transform(4326) ->
    roi
  potential_points <- st_sample(roi,size = points)
  return(potential_points)
}
