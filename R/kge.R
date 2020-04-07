# Kling-Gupta Efficiency
#'
#' @description The function to plot a hyetograph over a hydrograph based on ggplot functions.
#' @name hydrograph
#' @title Hydrograph plus Hyetograph
#' @param data tbl or Data frame object.
#' @param time Vector name of timestamp variable in `data`.
#' @param runoff Vector name of Runoff variable in `data`.
#' @param rain Vector name of Rainfall variable in `data`.
#' @param agg.time_unit Integrer. Default "3 hours". Time lag unit to aggregate
#'     rainfall bars by sum function. This will be showed as light blue line of
#'     hyetrograph.
#' @param agg.window Integrer. Default 8. Value of number o values from rainfall
#'     data to aggregate by rolling mean function. This will be showed as dark
#'     blue line of hyetrograph.
#' @param vsep numeric. Default 0.5. Value between 0 and 1 to scale the vertical
#'     separation of hyetograph from hydrograph. 1 means 100% of difference
#'     between min and max runoff.
#' @keywords hydrograph, rainfall-runoff plot, plot
#' @author Gabriel Gaona
#' @import dplyr ggplot2 rlang RcppRoll scales
#' @details ggplot2 object is generated. Other functions passed to \code{theme()}
#'     guides() or labs() can be used to customize the hydrograph.
#'
#' @export
#' @examples
#' # ploting soil moisture
#' data("soil_moisture")
#' #
#' hydrograph(data = soil_moisture,
#'     time = timestamp,
#'     runoff = runoff,
#'     rain = rain,
#'     agg.time_unit = "day",
#'     agg.window = 7) +
#'     ggplot2::theme_light()+
#'     ggplot2::theme(legend.position = "bottom")
#' @note Based on the request at \url{https://stackoverflow.com/q/42057832}

kge <- function(x, y){
  #Calculates Kling-Gupta efficiency
  #   Detailed explanation goes here

  sd_mod = sd(x)   #Standard deviation modelled
  sd_obs = sd(y)   #Standard deviation observations

  mn_mod = mean(x)  #Mean modeled
  mn_obs = mean(y)  #Mean observed

  r = cor(x, y)  #Correlation coefficient modeled and observed
  a = sd_mod / sd_obs    #ratios of standard deviations
  b = mn_mod / mn_obs    #ratios of means

  1 - sqrt(((r - 1) ^ 2) + ((a - 1) ^ 2) + ((b - 1) ^ 2)) #Kling-Gupta efficiency
}
