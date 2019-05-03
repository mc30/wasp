###############################################
# Analysis of peaks
###############################################


#' @title Get matrix of peaks
#' @description Creates a matrix with peak timing for each season.
#'
#' @param cases A matrix of case series.
#' @param start The start timepoint of each season.
#' @param period The period of seasons (e.g. 12 for monthly data, and 52 for weekly).
#' 
#' @return Returns a matrix with peak timing for each season and location.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
getPeaksMatrix <- function(cases, start, period = 12) {
  for (loc in 1:nrow(cases)) {
    peaks <- c()
    time <- c(start:(start + period))
    
    while(max(time) < ncol(cases)) {
      index <- which.max(cases[loc, time]) - 1
      peaks <- c(peaks, index)
      time <- time + period
    }
    
    if (loc == 1)
      pMat <- peaks
    else
      pMat <- rbind(pMat, peaks)
  }
  rownames(pMat) <- rownames(cases)
  colnames(pMat) <- colnames(cases)[seq(start, ncol(cases), length.out = ncol(pMat))]
  
  return(pMat)
}


#' @title Plot annual curves
#' @description Plots seasonal curves for each location with peaks and mean peak timing.
#'
#' @param cases A matrix of case series.
#' @param start The start timepoint of each season.
#' @param period The period of seasons (e.g. 12 for monthly data, and 52 for weekly).
#' @param perc CI range for the mean peak timing.
#' 
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
plotAnnualCurves <- function(cases, start, period = 12, perc = 0.95) {
  for (loc in 1:nrow(cases)) {
    timeStep <- period
    time <- c(start:(start + timeStep))
    
    if (period == 12)
      labels <- month.abb[c(start:period, 1:(start - 1))]
    else
      labels <- c(start:period, 1:(start - 1))
    
    plot(0, 0, xlim = c(0, period), ylim = c(0, max(cases[loc, ])), 
         main = row.names(cases)[loc], xlab = ifelse(period == 12, "Month", "Week"), ylab = "Cases", xaxt = "n",
         cex = 0)
    axis(1, at = 0:(period - 1), labels = labels)
    
    pal <- topo.colors(ncol(cases))
    
    peaks <- c()
    while(max(time) < ncol(cases)) {
      lines(time - min(time), cases[loc, time], col = pal[min(time)])
      index <- which.max(cases[loc, time]) - 1
      peaks <- c(peaks, index)
      points(index, cases[loc, index + min(time)], pch = 19, cex = 1, col = pal[min(time)])
      time <- time + timeStep
    }
    
    legend("topleft", legend = colnames(cases)[seq(start, ncol(cases), timeStep)], title = "Seasons' starts",
           col = pal[seq(start, ncol(cases), timeStep)], pch = 19, lty = 1, cex = 0.15)
    
    # ci <- qnorm(perc + (1 - perc) / 2) / sqrt(length(peaks)) # Coefficient to get CI from SD (normal distribution)
    ci <- qt(perc + (1 - perc) / 2, df = length(peaks) - 1) / sqrt(length(peaks)) # Coefficient to get CI from SD (t-distribution)
    
    abline(v = mean(peaks), col = "red", lty = 2)
    
    polygon(x = c(rep(mean(peaks) - ci * sd(peaks), 2), rep(mean(peaks) + ci * sd(peaks), 2)), 
            y = c(-max(cases), max(cases) * 2, max(cases) * 2, -max(cases)), border = NA, col = rgb(1, 0, 0, 0.1))
  }
}


#' @title Plot mean peak timing
#' @description Plot mean peak timing for all locations.
#'
#' @param pMat A matrix with peaks for each location and season.
#' @param start The start timepoint of each season.
#' @param period The period of seasons (e.g. 12 for monthly data, and 52 for weekly).
#' @param perc CI range for the mean peak timing.
#' @param bBootstrap Logical; if \code{TRUE} bootstrap will be used.
#' 
#' @return Returns a named vector of mean peaks for each location.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
plotMeanPeakTiming <- function(pMat, start, period = 12, perc = 0.95, bBootstrap = FALSE) {
  
  # Ordinary mean
  meanPeaks <- sapply(1:nrow(pMat), function(i) mean(pMat[i, ]))
  
  # ci <- qnorm(perc + (1 - perc) / 2) / sqrt(length(peaks)) # Coefficient to get CI from SD (normal distribution)
  ci <- qt(perc + (1 - perc) / 2, df = ncol(pMat) - 1) / sqrt(ncol(pMat)) # Coefficient to get CI from SD (t-distribution)
  
  sds <- sapply(1:nrow(pMat), function(i) sd(pMat[i, ]))
  
  upper <- meanPeaks + ci * sds
  lower <- meanPeaks - ci * sds
  
  # Bootstrapping
  if (bBootstrap) {
    mboot <- sapply(1:nrow(pMat), function(i) boot(pMat[i, ], function(x, k) mean(x[k]), R = 10000)$t[, 1])
    meanPeaks <- sapply(1:nrow(pMat), function(i) mean(mboot[, i]))
    
    lower <- sapply(1:nrow(pMat), function(i) quantile(mboot[, i], probs = (1 - perc) / 2))
    upper <- sapply(1:nrow(pMat), function(i) quantile(mboot[, i], probs = perc + (1 - perc) / 2))
  }
  
  names(meanPeaks) <- rownames(pMat)
  
  old.par <- par(mar = c(5.1 + 2.5, 4.1, 4.1, 2.1))
  
  upper <- upper[order(meanPeaks)]
  lower <- lower[order(meanPeaks)]
  meanPeaks <- sort(meanPeaks)
  
  plot(c(0, 0), range(upper, lower), cex = 0, 
       xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", main = paste0("Mean peak timing", ifelse(bBootstrap, " (bootstrapped)", "")), 
       xlim = c(-0.5, 0.5))
  
  if (period == 12)
    axis(2, at = 0:(period - 1), month.abb[c(start:period, 1:(start - 1))], las = 1, cex.axis = 0.8)
  else
    axis(2, at = 0:(period - 1), c(start:period, 1:(start - 1)), las = 1, cex.axis = 0.8)
  abline(h = 0:(period - 1), lty = 3, col = "lightgray")
  
  x <- seq(-0.5, 0.5, length.out = length(meanPeaks))
  points(x, meanPeaks, col = "red", pch = 19)
  for (i in order(row.names(cases)))
    lines(c(x[i], x[i]), c(lower[i], upper[i]), col = "red")
  axis(side = 1, at = x, labels = names(meanPeaks), las = 2, cex.axis = 0.6)
  par(old.par)
  
  return(meanPeaks)
}