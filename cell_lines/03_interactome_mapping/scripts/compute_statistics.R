# =============================================================================
# Function: computeStatistics
# =============================================================================
# Description: Computes z-score and one-tailed p-value for an observed network
#              statistic against a null distribution. Optionally plots the
#              distribution with the observed value marked.
#
# Args:
#   distribution : numeric vector — null distribution from permutations
#   observation  : numeric — observed value
#   xlab         : string — x-axis label for the plot
#   plot         : logical — whether to generate a density plot (default: TRUE)
#
# Returns:
#   data.frame with columns: m (mean), sd, z (z-score), pval
#
# Author: Federica D'Annunzio
# =============================================================================

computeStatistics <- function(distribution, observation, xlab, plot = TRUE) {

  m   <- mean(distribution, na.rm = TRUE)
  sd  <- sd(distribution, na.rm = TRUE)
  z   <- (observation - m) / sd
  pval <- pnorm(z, lower.tail = FALSE)  # one-tailed: H1 = greater than expected

  if (plot) {
    xmin <- min(m - 3 * sd, observation)
    xmax <- max(m + 3 * sd, observation)

    fun  <- function(x) dnorm(x, m, sd)

    par(xpd = FALSE)  # clip plotting to plot region
    curve_data <- curve(fun, xlim = c(xmin, xmax),
                        xlab = xlab, ylab = "Probability density", n = 1000)
    polygon(curve_data, col = "#00FF08")

    abline(v = observation, lty = 2, lwd = 2, col = "#FF0095")

    legend("topright",
           paste("Observation =", format(observation, digits = 3), "\n",
                 "p-value =", format(pval, digits = 1)),
           text.col = "#FF0095", text.font = 2, bty = "o",
           xjust = 0.5, yjust = 0.5)
  }

  data.frame(m = m, sd = sd, z = z, pval = pval)
}
