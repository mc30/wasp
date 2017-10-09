###############################################
# Auxillary functions
###############################################


#' @title Get colors for numeric vector values
#' @description Function to convert numbers into colors.
#'
#' @param x A numeric vector of values.
#' @param colors A vector of colors for the palette.
#' @param range A range of values for color picking.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
getColorsFromNumbers <- function(x, colors = c("blue", "white", "red"), range = range(x)) {
  mColors <- colorRamp(colors)((x - min(range)) / (max(range) - min(range)))
  return(rgb(mColors[, 1] / 255, mColors[, 2] / 255, mColors[, 3] / 255))
}