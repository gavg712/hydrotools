#' @importFrom ggplot2 ggproto Stat
#' @importFrom lubridate floor_date
#' @importFrom dplyr group_by summarise mutate
#' @import rlang
#' @export
StatRainfall <- ggproto("StatRainfall",
                        Stat,
                        compute_group = function(data,
                                                 scales,
                                                 agg.time.unit = "hour",
                                                 ysep = 0.5
                        ) {
                          lim <- max(ref.range) + diff(ref.range) * ysep

                          data %>%
                            dplyr::group_by(
                              x = lubridate::floor_date(x, by = agg.time.unit)) %>%
                            dplyr::summarise(rain = sum(rain, na.rm = na.rm)) %>%
                            dplyr::mutate(rain_scl = scales::rescale(
                              rain, to = c(lim, min(ref.range)))) %>%
                            dplyr::mutate(rain_rm = RcppRoll::roll_mean(
                              rain_scl, n = agg.window, align = "center",
                              fill = lim, na.rm = TRUE))
                        },

                        required_aes = c("x", "rain", "ref.range")
)

stat_rainfall <- function(mapping = NULL, data = NULL, geom = "segment",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, ...) {
    layer(
      stat = StatChull, data = data, mapping = mapping, geom = geom,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(na.rm = na.rm, ...)
    )
  }
