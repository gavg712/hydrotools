# hydrograph plus hyetograph
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

hydrograph <- function(data,
                       time = time,
                       runoff = Runoff,
                       rain = Rain,
                       agg.time_unit = "3 hours",
                       agg.window = 8,
                       vsep = 0.5) {

  lim <- max(dplyr::pull(data, {{ runoff }}), na.rm = TRUE) +
    {range(
      dplyr::pull(data, {{ runoff }}), na.rm = TRUE) %>%
       diff()} * vsep

  rain_data <- data %>%
    dplyr::group_by({{ time }} := lubridate::floor_date({{ time }},
                                                        unit = agg.time_unit)) %>%
    dplyr::summarise({{ rain }} := sum({{ rain }}, na.rm = TRUE)) %>%
    dplyr::mutate(rain_scl = scales::rescale({{ rain }},
                                      to = c(lim,
                                             min(
                                               dplyr::pull(data,
                                                           {{ runoff }}),
                                                 na.rm = TRUE))))

  # base plot object
  gg <- ggplot2::ggplot(data = rain_data, aes(x = {{ time }}))

  # layers for hyetograph
  gg <- gg +
    ggplot2::geom_segment(mapping =
                              ggplot2::aes(xend = {{ time }},
                                           y = rain_scl, yend = lim,
                                           color = "Rainfall"),
                            alpha = 0.5) +
    ggplot2::geom_line(aes(y = RcppRoll::roll_mean(rain_scl, n = agg.window,
                                align = "center", fill = lim),
                  color = paste0("Running mean Rainfall [", agg.time_unit, "]")))

  # layers for hydrogrpah
  gg <- gg +
    ggplot2::geom_line(data = data, mapping = aes(y = {{ runoff }}, color = "Runoff"))

  # Adding a second axis
  gg <- gg +
    ggplot2::scale_y_continuous(expression(
      paste("Runoff [m"^3, "/s]")),
      sec.axis = ggplot2::sec_axis(
        trans = ~scales::rescale(.,
                                 to = base::rev(
                                   base::range(
                                     dplyr::pull(rain_data,
                                                 {{ rain }}))),
                                 from = c(0, lim)),
        name = paste0("Rainfall [mm/", agg.time_unit, "]")))
  gg <- gg +
    ggplot2::scale_color_manual(
      values = rlang::set_names(
        c("#3f88c5","#19364e", "#f97171"),
        c("Rainfall", paste0("Running mean Rainfall [", agg.time_unit, "]"),
          "Runoff"))
      )

  gg <- gg +
    guides(color = guide_legend(title = ""))

  gg
}
