## Level (1-a) Wilson confidence interval for proportion
## WILSON, E. B. 1927. Probable inference, the law of succession,
## and statistical inference. Journal of the American Statistical
## Association 22: 209-212.
#adapted from https://www.r-bloggers.com/reference-chart-for-precision-of-wilson-binomial-proportion-confidence-interval/
# n is the absolute number of paricipants, p is the proportion

WilsonBinCI <-  function(n, p, a=0.05) {
  z <- qnorm(1-a/2,lower.tail=FALSE)
  l <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 + 
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
  u <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 -
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
list(lower=l, upper=u, width=u-l)
}

# this version requires you to provide a vector x, and specify which values indicate success
WilsonBinCIVec <-  function(x, markerForSuccess, a=0.05) {
  z <- qnorm(1-a/2,lower.tail=FALSE)
  n <- length(x)
  p <- mean(x == markerForSuccess)
  l <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 + 
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
  u <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 -
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
  list(lower=l, upper=u, width=u-l)
}


# somewhat hackish solution to:
# https://twitter.com/EamonCaddigan/status/646759751242620928
# based mostly on copy/pasting from ggplot2 geom_violin source:
# https://github.com/hadley/ggplot2/blob/master/R/geom-violin.r

library(ggplot2)
library(dplyr)


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )
