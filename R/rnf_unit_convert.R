#' Runoff unit conversion
#'
#' Runoff units can be compared with precipitation by transforming it scale into
#' comparable scales. Usually rainfall is measured in level terms meaning the
#' rise of water level in a theoretical volume of base equal to one square meter.
#' @name rnf_unit_convert
#'
#' @description The function to transform runoff variables from/to cubic meters
#' per second to water level in mm per unit time on a square meter.
#' @title Runoff conversion Volumen per time from/to level
#' @param x Vector numeric. Runoff variable
#' @param from Character. String to convert 'vol' or 'lev'. Default 'vol'
#' @param to.time.unit Character. String containing the time units of conversion.
#' See more in details
#' @param area numeric, catchment area  in suare meters.
#' @keywords conversion, rainfall-runoff, plot
#' @author Gabriel Gaona
#' @importFrom stringr str_extract
#' @importFrom stringr str_detect
#' @importFrom readr parse_number
#' @details
#'
#' `from` argument require `"vol"` meaning "volume per second" (m3/s) or `"lev"`
#' meaning "level on a time over area" (mm day).
#'
#' `to.time.unit` is a time fraction string i.e. "day", "3 hour", "10 minutes", etc.
#'
#' @export
#' @examples
#'
#' # In a catchment area of 1.4 km2 a measure of water outcomes by runoff is
#' # 0.023 'm3/s'
#' area <- 1.4 * 1000000 #convert km2 to m2
#' r_vol <- 0.023
#' # let's convert to "mm dia"
#' r_lev <- rnf_unit_convert(r_vol, from = "vol",
#'                           to.time.unit = "1 day", area = area)
#' # Reverse process
#' r_lev <- rnf_unit_convert(r_lev, from = "lev",
#'                           to.time.unit = "1 day", area = area)

rnf_unit_convert <- function(x,
                             from = "vol", #cubic meters per second
                             to.time.unit = "day",
                             area # runoff basin area in square meters
                             ) {
  time.unit <- stringr::str_extract(
    to.time.unit,
    "second|minute|hour|day|month|year|decade") %>%
    paste("mm", .)

  if(stringr::str_detect(to.time.unit, "[0-9]+")){
    ntime <- readr::parse_number(to.time.unit)
  } else {
    ntime <- 1
  }
  secs <- c("mm second" = 1,
            "mm minute" = 60,
            "mm hour" = 3600,
            "mm day" = 86400,
            "mm month" = 2592000,
            "mm year" = 31104000,
            "mm decade" = 311040000)
  if(from == "vol") {
    (x * 1000 * secs[time.unit] * ntime) / area
  } else if (from == "lev") {
    (x * area) / (1000 * secs[time.unit] * ntime)
  } else {
    stop("'from' argument: \n\t'vol' (volume/time) or 'lev' (mm over sq meter) are allowed")
  }
}
