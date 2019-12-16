#' Rainfall Scale and window aggregation
#' @name rainfall_agg
#'
#' @description The function to scale and agregate rainfall data for hydrograph
#' @title Rainfall: Scale and window aggregation
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

StatRainfall <- ggplot2::ggproto(
  "StatRainfall",
  ggplot2::Stat,
  required_aes = c("x", "y"),
  setup_params = function(data, params){
    params$xend <- params$x
    params$yend <- params$ytop
  },

  compute_group = function(data,
                           scales,
                           scale.rng = NULL,
                           agg.time_unit = NULL,
                           agg.window = NULL,
                           ytop = 0.5) {

    if(!rlang::is_null(agg.time_unit)) {
      data <- data %>%
        dplyr::group_by(
          y = lubridate::floor_date(y, unit = agg.time_unit)) %>%
        dplyr::summarise(rain = sum(y, na.rm = TRUE))
    }

    if(rlang::is_null(scale.rng)) {
      scale.rng <- range(data$y, na.rm = TRUE)
    }

    data$rain_scl <- scales::rescale(data$y, to = c(ytop, min(scale.rng)))

    if(!rlang::is_null(agg.window)){
      data$rain_rm <- RcppRoll::roll_mean(data$rain_scl, n = agg.window,
                                      align = "center", fill = ytop,
                                      na.rm = FALSE)
    }
    return(data)
  }
)

#' @import ggplot2
#' @name stat_rainfall
#' @export
stat_rainfall <- function(mapping = NULL, data = NULL,
                          geom = "Segment",
                          position = "identity",
                          show.legend = NA,
                          inherit.aes = TRUE,
                          scale.rng = NULL,
                          agg.time_unit = NULL,
                          agg.window = NULL,
                          ytop = 0.05, ...) {
  ggplot2::layer(
    stat = StatRainfall, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(scale.rng = scale.rng,
                  agg.time_unit = agg.time_unit,
                  agg.window = agg.window,
                  ytop = ytop, ...)
  )
}

rainfall_agg <- function(data,
                     time,
                     rain = y,
                     scale.rng = NULL,
                     agg.time_unit = NULL,
                     agg.window = NULL,
                     vsep = 0.5) {

  if(!rlang::is_null(agg.time_unit)) {
  data <- data %>%
    dplyr::group_by(
      {{ time }} := lubridate::floor_date({{ time }}, unit = agg.time_unit)) %>%
    dplyr::summarise({{ rain }} := sum({{ rain }}, na.rm = TRUE))
  }

  if(rlang::is_null(scale.rng)) {
    scale.rng <- range(data[[rlang::ensym(rain)]], na.rm = TRUE)
  }
  lim <- max(scale.rng) + diff(scale.rng) * vsep

  data <- dplyr::mutate(
    data,
    rain_scl = scales::rescale({{ rain }}, to = c(lim, min(scale.rng))))

  if(!rlang::is_null(agg.window)){
    dplyr::mutate(
      data,
      rain_rm = RcppRoll::roll_mean(rain_scl, n = agg.window,
                                    align = "center", fill = lim,
                                    na.rm = FALSE)

    )
  }
  return(data)
}
