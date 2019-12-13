#' Pipe operator
#'
#' Like dplyr, hydroutils also uses the pipe function, \code{\%>\%} to turn
#' function composition into a series of imperative statements.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs A chain process to concatenate
#' @examples
#' # Instead of
#' ggplot(mtcars, aes(mpg, wt)) + geom_points()
#' # you can write
#' mtcars %>% ggplot(aes(mpg, wt)) + geom_points()
NULL
