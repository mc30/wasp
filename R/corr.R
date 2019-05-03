###############################################
# Correlation analysis
###############################################


# 

#' @title Spatial correlation function as in 'ncf' package
#' @description .
#'
#' @param mat A matrix.
#' @param coords .
#' @param resamp .
#' @param main .
#' @param \dots Additional graphical parameters.
#' 
#' @return .
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@pasteur.fr}).
#' 
#' @export
spatialCorrelationFunction <- function(mat, coords, resamp = 100, main = "", ...) {
  fit1 <- Sncf(x = coords[, 1], y = coords[, 2], z = mat, resamp = resamp, na.rm = TRUE)
  # plot.Sncf(fit1)
  # s <- summary.Sncf(fit1)
  
  plot(fit1$real$predicted$x, fit1$real$predicted$y, 
       # fit1$boot$boot.summary$predicted$x, fit1$boot$boot.summary$predicted$y[6, ],
       ylim = c(0, 1), 
       type = "l", lty = 1, lwd = 2, col = "blue",
       xlab = "Distance (km)", ylab = "Correlation", main = main, ...)
  # lines(fit1$real$predicted$x, fit1$real$predicted$y)
  
  lines(fit1$boot$boot.summary$predicted$x, fit1$boot$boot.summary$predicted$y[6, ], lty = 2) # 0.5
  lines(fit1$boot$boot.summary$predicted$x, fit1$boot$boot.summary$predicted$y[2, ], lty = 2) # 0.025
  lines(fit1$boot$boot.summary$predicted$x, fit1$boot$boot.summary$predicted$y[10, ], lty = 2) # 0.975
  
  # i <- 4
  # lines(fit1$boot$boot.summary$predicted$x, fit1$boot$boot.summary$predicted$y[i,]) #
  # lines(fit1$boot$boot.summary$predicted$x, fit1$boot$boot.summary$predicted$y[12 - i,]) #
  
  abline(h = fit1$real$cbar, col = "red")
  
  # text(0, fit1$real$cbar - 0.05, "Country wide correlation", pos = 4, col = "red", cex = 1)
  
  # # threshold <- fit1$real$predicted$x[which.min(abs(fit1$real$predicted$y - fit1$real$cbar))]
  threshold <- fit1$real$cbar.intercept
  
  # # abline(h = fit1$boot$boot.summary$cbar, col = "red", lty = 3)
  # abline(h = fit1$boot$boot.summary$cbar[2], col = "red")
  # 
  # abline(v = fit1$boot$boot.summary$cbar.intercept, lty = 3)
  # threshold <- fit1$boot$boot.summary$cbar.intercept[2]
  
  abline(v = threshold, lty = 1)
  
  mtext(signif(threshold, 3), at = threshold)
  
  legend("topright", legend = c("Countrywide correlation", "Estimated correlation from data", "Bootstrapped correlation (95% envelope)"), 
         lty = c(1, 1, 2), col = c("red", "blue", "black"), cex = 0.7)
}