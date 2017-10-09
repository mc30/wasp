###############################################
# Wavelets: phase angles and phase differences
###############################################


# @encoding UTF-8
#' @title Calculate wavelet spectra and phase angle matrix for time series
#' @description Computes wavelet spectra and creates a matrix of phase angles.
#'
#' @param bwn A matrix of values for which the wavelet spectrum is desired.
#' @param locs A vector of rows (i.e. locations) to be used in the analysis.
#' @param units Number of time units in one period (e.g. months in a year) (X axis)
#' @param periodStep Resolution of time component (i.e. periods) (Y axis)
# @param \dots Additional arguments (currently ignored)
#' 
#' @return Returns a matrix with phase angles.
#' @details Rows with uniform zero incidence will be ignored and excluded from the resulting list.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#'
#' @export
calculateWaveletSpectrum <- function(bwn, locs = 1:nrow(bwn), units = 12, periodStep = 1 / 100) {
  print("Calculating wavelet power spectrum...")
  
  wts <- vector("list", length(locs))
  k <- 0
  for (i in locs) {
    print(rownames(bwn)[i])
    
    if (length(which(bwn[i, ] == 0)) == length(bwn[i, ])) { # if all the case numbers are zeros
      k <- k + 1
      next
    }
    
    # Construct object for wavelet input
    tsfoc <- cbind(1:ncol(bwn), bwn[i, ])
    
    ##show the focal wavelet plus ts plot
    wtfoc <- wt(tsfoc, dj = periodStep)
    wtfoc$t <- wtfoc$t / units
    wtfoc$period <- wtfoc$period / units
    wtfoc$coi <- wtfoc$coi / units
    
    k <- k + 1
    wts[[k]] <- wtfoc 
  }
  names(wts) <- rownames(bwn)[locs]
  
  # Remove NULLs
  if (length(which(sapply(wts, is.null))) > 0)
    wts <- wts[-which(sapply(wts, is.null))]
  
  print("... Done!")
  
  return(wts)
}


#' @title Plot wavelet spectrum
#' @description Plots wavelet spectrum
#'
#' @param wtfoc Wavelet spectrum object.
#' @param units Number of time units in one period (e.g. months in a year) (X axis)
#' @param period The period we are interested in.
#' @param periodRange The range of periods to plot (Y axis).
#' @param title Main title for the plot.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' @export
plotWaveletSpectrum <- function(wtfoc, units = 12, period = 1, periodRange = c(0.2, 4), title = "") {
  plot(wtfoc, plot.cb = FALSE, plot.phase = FALSE, main = title,
       xlab = "Time", xaxt = "n", las = 1, ylim = periodRange)
  abline(h = log(period, 2), lty = 2) # line for the period of interest
  seq <- seq(1, ncol(wtfoc$power), units)
  axis(1, at = (1:ncol(wtfoc$power) / units)[seq], labels = names(wtfoc$t)[seq])
}


#' @title Plot wavelet spectrum summary graph
#' @description Plots a summary graph with three subgraphs: timeline of case series, wavelet spectrum and 
#'     overall wavelet power.
#'
#' @param wts A list of wavelet spectrum objects.
#' @param cases A matrix with case series.
#' @param i A row number (i.e. location) for plotting.
#' @param period The period we are interested in.
#' @param periodRange The range of periods to plot.
#' @param bLog Logical; if \code{TRUE} adds a line of logged cases.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
plotWaveletGraph <- function(wts, cases, i, period = 1, periodRange = c(0.2, 4), bLog = FALSE) {
  old.par <- par(mar = c(2, 4.1, 4.1, 0))
  
  m <- rbind(c(1, 4), c(2, 3))
  layout(m, widths = c(3, 1), heights = c(1, 2))
  
  
  ### First graph: normalised cases
  
  maxCases <- max(cases[i, ])
  plot(1:ncol(cases), 10 * cases[i, ] / maxCases, type = "l", xaxs = "i", xaxt = "n", yaxt = "n", las = 1, xlab = " ", 
       ylab = "Cases", ylim = c(-1, 10))
  interval <- 0:5 * 2
  axis(2, at = interval, labels = maxCases * interval, las = 1, cex.axis = 0.6)
  
  # Add log cases
  if (bLog) {
    lines(1:ncol(cases), log(cases[i, ] + 1), lty = 3, col = "blue")
    labels <- (log(maxCases * interval + 1) - mean(log(maxCases * interval + 1))) / sd (log(maxCases * interval + 1))
    axis(4, at = 0:5 * 2, labels = signif(labels, 2), las = 1, cex.axis = 0.6, col.axis = "blue")
    legend("topleft", legend = c("Original cases", "Log-scaled and normalised"), lty = c(1, 3), col = c("black", "blue"), cex = 0.5, bty = "n")
  }
  title(main = paste("Cases in ", row.names(cases)[i], sep = ""))
  
  
  ### Second graph: wavelet spectrum
  par(mar = c(5.1, 4.1, 0, 0))
  plotWaveletSpectrum(wts[[i]])
  
  
  ### Third graph: overall wavelet power
  par(mar = c(5.1, 2, 0, 2.1))
  wtfoc <- wts[[i]]
  
  plot(rowSums(wtfoc$power), nrow(wtfoc$power):1,
       type = "l", xlab = "Power", ylab = "", las = 1, 
       yaxt = "n",
       # yaxs = "i", 
       ylim = c(nrow(wtfoc$power) - which.min(abs(wtfoc$period - periodRange[2])), nrow(wtfoc$power) - which.min(abs(wtfoc$period - periodRange[1])))
       # ylim = c(0.2, 4)
  )
  # seq <- seq(nrow(wtfoc$power) - which.min(abs(wtfoc$period - 4)), nrow(wtfoc$power) - which.min(abs(wtfoc$period - 0.2)), length.out = 5)
  # axis(2, at = seq, wtfoc$period[seq], las = 1)
  
  abline(h = nrow(wtfoc$power) - which.min(abs(wtfoc$period - period)), lty = 2) # line for the period of interest
  # axis(2, at = 1:ncol(bwn) / units, labels = colnames(bwn))
  
  par(mfrow = c(1, 1))
  par(old.par)
}


#' @title Create a matrix of phase angles
#' @description Creates a matrix of phase angles for a specified period
#'
#' @param wts List of wavelet transform outputs.
#' @param period The period of interest.
#' @param significantPowerThreshold A threshold for power, e.g. if equals to 5, phase angles with corresponding power < 5 will be NAs.
#' 
#' @return Returns a matrix with phase angles.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
createPhaseMatrix <- function(wts, period = 1, significantPowerThreshold = 0) {
  ncol <- ncol(wts[[1]]$power)
  
  # Matrix with phases
  phaseMatrix <- matrix(rep(NA, each = ncol * length(wts)), length(wts), ncol, byrow = TRUE)
  rownames(phaseMatrix) <- names(wts)
  colnames(phaseMatrix) <- colnames(wts[[1]]$power)
  
  for (i in 1:length(wts)) {
    index <- which.min(abs(wts[[i]]$period - period)) # Index of closest period
    
    # Exclude non-significant timepoints (those that have weak signal of the selected period)
    nonSignificant <- which(wts[[i]]$power[index, ] < significantPowerThreshold)
    wts[[i]]$phase[index, nonSignificant] <- NA
    
    phaseMatrix[i, ]  <- wts[[i]]$phase[index, ] # Store phase angles in phaseMatrix
  }
  
  return(phaseMatrix)
}


#' @title Create power lines
#' @description Creates power lines 
#'
#' @param wts List of wavelet transform outputs.
#' 
#' @return Returns a matrix where each row is a power line.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#'
#' @export
createPowerLines <- function(wts) {
  powerLines <- matrix(rep(NA, each = length(wts)), length(wts), nrow(wts[[1]]$power), byrow = TRUE)
  colnames(powerLines) <- wts[[1]]$period
  rownames(powerLines) <- names(wts)
  
  for (i in 1:length(wts))
    powerLines[i, ] <- rowMeans(wts[[i]]$power, na.rm = TRUE)
  
  return(powerLines)
}


#' @title Plot power lines
#' @description Plots power lines
#'
#' @param powerLines Matrix with wavelet power vectors for each location.
#' 
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#'
#' @export
plotPowerLines <- function(powerLines) {
  plot(0, 0, cex = 0, xlim = range(as.numeric(colnames(powerLines))), ylim = c(0, max(powerLines)), las = 1,
       xlab = "Period", ylab = "Power")
  for (i in 1:nrow(powerLines))
    lines(as.numeric(colnames(powerLines)), powerLines[i, ])
}


#' @title Correct phase angles
#' @description Function to correct for spikes while computing phase difference.
#'
#' @param vector A vector of phase angles.
#' 
#' @return A vector with phase angles that fall within [-pi, pi] range.
#' @details While perfrming arithmetical operations with phase angles, one could end up with values outside [-pi, pi] range.
#'     This function adds or subtracts additional \code{2 * pi} for correction.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
correctPhase <- function(vector) {
  for (i in 1:length(vector)) {
    if (is.na(vector[i])) {
      vector[i] <- NA
      next
    }
    
    while (vector[i] > pi)
      vector[i] <- vector[i] - 2 * pi
    
    while (vector[i] < -pi)
      vector[i] <- vector[i] + 2 * pi
  }
  return(vector)
}


# doCrossWavelet <- function(bwn, i, j) {
#   # pick out two time-series to focus on
#   tsfoc <- cbind(1:ncol(bwn), bwn[i, ])
#   tscomp <- cbind(1:ncol(bwn), bwn[j, ])
#   
#   
#   ##cross wavelet 
#   #computation (correct for time step)
#   xwt.compfoc <- xwt(tscomp, tsfoc)
#   # xwt.compfoc$t <- xwt.compfoc$t / 12
#   xwt.compfoc$period <- xwt.compfoc$period / 12
#   xwt.compfoc$coi <- xwt.compfoc$coi / 12
#   
#   #function to return vector on 0-1
#   sc01 <- function(x) { 
#     x <- x - min(x, na.rm=TRUE)
#     x <- x / max(x) 
#     return(x)
#   }
#   
#   #plotting - showing signifcant coherence
#   par(fig = c(0, 1, 0, 0.8))
#   plot(xwt.compfoc, xlab = "Years", plot.cb = F, plot.phase = TRUE)
#   
#   par(fig = c(0, 1, 0.6, 1), new = TRUE)
#   plot(tsfoc[, 1], -sc01(tscomp[, 2]), type = "l", xaxs = "i", xaxt = "n", col = 1, 
#        main = paste(rownames(bwn)[i], " - ", rownames(bwn)[j], sep = ""),
#        ylim = c(-1, 1), xlab = " ")
#   abline(h = 0)
#   lines(tsfoc[, 1], sc01(tsfoc[, 2]), col = 2)
# }


# crossWavelet <- function(bwn) {
#   for (i in 1:nrow(bwn)) {
#     for (j in i:nrow(bwn)) {
#       doCrossWavelet(bwn, i, j)
#     }
#   }
# }


#' @title Calculate phase differences from phase matrix
#' @description Calculates pairwise differences between phase angles.
#' 
#' @param phaseMatrix A matrix with phase angles over time for considered locations.
#' @param refLoc A row number for reference location.
#' @param bPlot Logical; if \code{TRUE} plots are generated for each pair of locations.
#' 
#' @return  A matrix with pairwise phase differences.
#' @details Phase differences are calculated between the reference location and all other locations. 
#'     The row for the reference locations will contain NAs.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#'
#' @export
calculatePhaseDiff <- function(phaseMatrix, refLoc, bPlot = TRUE) {
  phaseDiffMatrix <- matrix(rep(NA, each = ncol(phaseMatrix)), nrow(phaseMatrix), ncol(phaseMatrix), byrow = TRUE)
  rownames(phaseDiffMatrix) <- rownames(phaseMatrix)
  
  old.par <- par(mfrow = c(2, 1))
  for (i in refLoc) {
    for (j in 1:nrow(phaseMatrix)) {
      if (i == j)
        next
      
      locI <- row.names(phaseMatrix)[i]
      locJ <- row.names(phaseMatrix)[j]
      
      if (bPlot) {
        plot(phaseMatrix[i,], type = "l", xlab = "Month", ylab = "Phase (radians)", las = 1, col = "red")
        lines(phaseMatrix[j,], col = "blue")
        abline(h = 0, lty = 2)
        
        title(bquote(paste("Phase angles for ", phantom(.(locI)), " and ", phantom(.(locJ)))), col.main = "black")
        title(bquote(paste(phantom("Phase angles for "), .(locI), phantom(paste(" and ", .(locJ))))), col.main = "red")
        title(bquote(paste(phantom(paste("Phase angles for ", .(locI), " and ")), .(locJ))), col.main = "blue")
      }
      
      vec <- phaseMatrix[i, ] - phaseMatrix[j, ]
      phaseDiffMatrix[j, ] <- correctPhase(vec)
      
      if (bPlot) {
        plot(phaseDiffMatrix[j, ], type = "l", las = 1,
             ylim = c(-pi, pi), col = "black", xlab = "Month", ylab = "Phase (radians)",
             main = paste("Phase difference", sep = ""))
        abline(h = 0, lty = 2)
      }
    }
  }
  par(old.par)
  
  return(phaseDiffMatrix)
}


#' @title Perform wavelet clustering
#' @description Compares wavelet spectra and generates an hierarchical tree for locations.
#'
#' @param bwn A matrix of time series.
#' 
#' @return A matrix of pairwise distances between wavelet spectra.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
waveletClustering <- function(bwn) {
  t <- cbind(1:ncol(bwn), bwn[1, ])
  wt.t1 <- wt(t)
  
  ## Store all wavelet spectra into array
  w.arr <- array(dim = c(nrow(bwn), NROW(wt.t1$wave), NCOL(wt.t1$wave)))
  
  for (i in 1:nrow(bwn)) {
    t <- cbind(1:ncol(bwn), bwn[i, ])
    spectrum <- wt(t)
    w.arr[i, , ] <- spectrum$wave
  }
  
  ## Compute dissimilarity and distance matrices
  w.arr.dis <- wclust(w.arr)
  
  # Distance between wavelet spectra
  # wdist(w.arr[1, ,], w.arr[2, ,])
  
  # Matrix of distances
  # w.arr.dis$dist.mat
  
  return(w.arr.dis)
}


#' @title Plot wavelet cluster dendrogram
#' @description Generates and plots an hierarchical tree for locations based on pairwise distances between their
#'     wavelet spectra.
#'
#' @param w.arr.dis A matrix with distances between wavelet spectra.
#' @param k A number of clusters of locations.
#' @param labels A vector of labels for the tree tips (i.e. locations).
#' @param title A title for the generated graph.
#' 
#' @return An object of class \code{hclust} which describes the resulting tree.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#'
#' @import grDevices
#' @import stats
#' @import graphics
#' 
#' @export
plotWaveletClusterDendrogram <- function(w.arr.dis, k = 4, labels = rownames(w.arr.dis$dist.mat), title = "Wavelet cluster dendrogram") {
  hc <- hclust(w.arr.dis$dist.mat, method = "ward.D")
  
  clust <- cutree(hc, k)
  palette <- rainbow(k)
  
  plot(hc, 
       # type = "fan", 
       # tip.color = palette[clust], 
       labels = labels,
       main = title,
       xlab = "",
       # label.offset = 1, 
       cex = 0.7)
  
  return(hc)
}
