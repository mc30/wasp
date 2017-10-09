###############################################
# Phase analysis
###############################################


#' @title Calculate average phase angles
#' @description Creates a matrix with mean phase differences between regions for a moving time period.
#'
#' @param phaseMatrix A matrix of phase angles.
#' @param window The size of the moving average window.
#' @param bSort Logical; if \code{TRUE} sorts the output matrix according to the average phase (over time).
#' 
#' @return Returns a matrix with average phase angle differences from other locations.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
createSlidingPhaseMatrix <- function(phaseMatrix, window = 12, bSort = TRUE) {
  print("Calculating sliding phase matrix...")
  
  time <- 1:ncol(phaseMatrix)
  phaseSlideMatrix <- matrix(rep(NA, each = nrow(phaseMatrix), ncol(phaseMatrix)), nrow(phaseMatrix), byrow = TRUE)
  
  # colnames(phaseSlideMatrix) <- row.names(phaseMatrix)
  row.names(phaseSlideMatrix) <- row.names(phaseMatrix)
  
  prog.bar <- txtProgressBar(min = 0, ncol(phaseMatrix), style = 3)
  for (t in 1:ncol(phaseMatrix)) {
    for (i in 1:nrow(phaseMatrix)) {
      vals <- c()
      for (j in 1:nrow(phaseMatrix)) {
        if (i == j) {
          vals <- c(vals, 0)
          next
        }
        
        time <- max(1, t - window / 2):min(t + window / 2, ncol(phaseMatrix))
        
        # vec <- phaseMatrix[i, time] - phaseMatrix[j, time]
        vec <- phaseMatrix[j, time] - phaseMatrix[i, time]
        val <- mean(correctPhase(vec), na.rm = TRUE)
        vals <- c(vals, val)
      }
      
      phaseSlideMatrix[i, t] <- mean(vals, na.rm = TRUE)
    }
    
    setTxtProgressBar(prog.bar, t)
  }
  close(prog.bar)
  
  if (bSort)
    phaseSlideMatrix <- phaseSlideMatrix[order(rowMeans(phaseSlideMatrix), decreasing = T), ] # Sort regions according to their mean phase lags
  
  return(phaseSlideMatrix)
}


#' @title Plot mean phase lags from other locations
#' @description Plot mean phase lags from other locations.
#'
#' @param slideMat A matrix of phase angles.
#' @param sel Vector of locations to use for the plot.
#' @param colors Vector of colors to create a palette for different locations.
#' @param legendSize Size of the legend. 
#' @param legendNcol Number of columns in the legend. 
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
plotMeanPhaseLagOverTime <- function(slideMat, sel = 1:nrow(slideMat), colors = c("red", "yellow", "green", "blue"),
                                     legendSize = 1, legendNcol = 1) {
  plot(1:ncol(slideMat), rep(0, ncol(slideMat)), cex = 0, ylim = c(-pi, pi), xlim = c(1, ncol(slideMat)),
       xlab = "Time", ylab = "Phase lag from other regions")
  
  meanPhases <- rowMeans(slideMat)
  pal <- getColorsFromNumbers(meanPhases, colors, range(meanPhases))
  
  for (i in sel)
    lines(1:ncol(slideMat), slideMat[i, ], col = pal[i])
  
  legend("bottomright", rownames(slideMat)[order(meanPhases, decreasing = T)][sel], lty = 1, 
         col = pal[order(meanPhases, decreasing = T)][sel], text.col = pal[order(meanPhases, decreasing = T)][sel], 
         cex = legendSize, ncol = legendNcol)
}


#' @title Plot mean phase lags averaged over time
#' @description Plot mean phase lags averaged over time.
#'
#' @param slideMat A matrix of phase angles.
#' @param perc Percentage envelope for CI around the mean.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
plotAveragePhaseLag <- function(slideMat, perc = 0.95) {
  sds <- c()
  for (i in 1:nrow(slideMat))
    sds <- c(sds, sd(slideMat[i, ]))
  
  meanPhases <- rowMeans(slideMat)
  
  sds <- sds[order(meanPhases)]
  meanPhases <- sort(meanPhases)
  
  plot(meanPhases, ylim = c(-pi, pi), ylab = "Average phase lag from other locations", xlab = "", xaxt = "n", las = 1)
  
  ci <- qnorm(perc + (1 - perc) / 2) / sqrt(ncol(slideMat)) # Coefficient to get CI from SD (normal distribution)
  ci <- qt(perc + (1 - perc) / 2, df = ncol(slideMat) - 1) / sqrt(ncol(slideMat)) # Coefficient to get CI from SD (t-distribution)
  
  
  points(1:length(meanPhases), meanPhases + ci * sds, col = "red", pch = 19, cex = 0.5)
  points(1:length(meanPhases), meanPhases - ci * sds, col = "red", pch = 19, cex = 0.5)
  
  for (i in 1:length(meanPhases))
    lines(c(i, i), c(meanPhases[i] - ci * sds[i], meanPhases[i] + ci * sds[i]))
  
  mtext(names(meanPhases), 1, at = 1:length(meanPhases), las = 2, line = 0.2, cex = 0.5)
  
  return(meanPhases)
}