#' Fast solution to obtain extrema of RMPW estimator in causal decomposition
#' problems
#'
#' @param G Indicator of unmodifiable group or characteristic
#' @param Y Outcome
#' @param gamma Sensitivity parameter (log of MSM Lambda)
#' @param w fitted RMPW weights
#' @param estimand Target estimand, either "point", "reduction" or "residual"
#' @param stab Logical indicating whether to use stabilized weights
#' @param RD logical indicating whether to use risk diff or risk ratio (RR)
#'
#' @return Extrema (an interval) of the RMPW estimator.
#'
#' @export

getExtrema <- function(G, Y, gamma = 0, w, estimand = "point", stab = TRUE, RD = TRUE) {
  estimand <- match.arg(estimand, c("point", "reduction", "residual"))
  # assume observed rmpw weights already computed
  mu1 <- mean(Y[G == 1])
  mu0 <- mean(Y[G == 0])
  if (stab) {
    mu_10 <- sum(w * Y * G) / sum(w * G)
  } else {
    mu_10 <- mean(w * Y * G)
  }
  # message("mu1: ", round(mu1, 2))
  # message("mu0: ", round(mu0, 2))

  state <- switch(estimand, point = round(mu_10, 3),
                  reduction = round(mu1 - mu_10, 3),
                  residual = round(mu_10 - mu0, 3))

  message("Parameter of interest: ", estimand)
  message("Estimate: ", state)
  # prints the value of the estimator before bootstrapping

  w <- w[G != 0]
  Y <- Y[G != 0]
  G <- G[G != 0]

  w <- w[order(-Y)]
  Y <- Y[order(-Y)]
  G <- G[order(-Y)]

  ## maximization
  num.each.low <- G * Y * (exp(-gamma) * w)
  num.each.up <- G * Y * (exp(gamma) * w)
  num <- c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)

  if (stab) {
    den.each.low <- G * (exp(-gamma) * w)
    den.each.up <- G * (exp(gamma) * w)
  } else {
    den.each.low <- G * (exp(-gamma))
    den.each.up <- G * (exp(gamma))
  }
  den <- c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)

  maximum <- max(num / den)
  ## print(den[which.max(num/den)] / n)

  ## minimization
  num <- c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
  den <- c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
  minimum <- min(num / den)
  ## print(den[which.min(num/den)] / n)

  out <- c(minimum, maximum)
  if (RD) {
    return(switch(estimand,
           point = out,
           reduction = mu1 - rev(out),
           residual = out - mu0))
  } else {
    return(switch(estimand,
           point = out,
           reduction = mu1 / rev(out),
           residual = out / mu0))
  }
}
