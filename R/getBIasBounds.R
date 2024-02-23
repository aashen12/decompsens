#' Function to obtain bias bounds for amplification
#' @export

getBiasBounds <- function(G, Z, XA, XN, Y, Lambda, estimand, trim = 0.05,
                          allowable = FALSE, stab = TRUE) {
  ####################################################################
  # Get bounds on bias = \mu^*_10 - \hat{mu}_10
  # lower bound = inf h \hat{mu}_10^{(h)} - \hat{mu}_10
  # upper bound = inf h \hat{mu}_10^{(h)} - \hat{mu}_10
  #
  # input:
  #   Z: treatment vector
  #   X: covariate matrix
  #   Y: outcome vector
  #   Lambda: sensitivity parameter
  #
  # output:
  #   bounds: out[[1]]
  #     bounds[1]: upper bound = sup h \hat{mu}_0^{(h)} - \hat{mu}_0
  #     bounds[2]: lower bound = inf h \hat{mu}_0^{(h)} - \hat{mu}_0
  #   mu_0_hat: out[[2]]. estimate of E[Y(0)|Z==1]
  #   eb.out$w: out[[3]]. Estimated balancing weights
  #
  ####################################################################

  mu1 <- mean(Y[G == 1])
  mu0 <- mean(Y[G == 0])

  w_rmpw <- estimateRMPW(G = G, Z = Z, Y = Y, XA = XA, XN = XN,
                         trim = trim, allowable = allowable)

  if (stab) {
    mu_10 <- sum(w * Y * G) / sum(w * G)
  } else {
    mu_10 <- mean(w * Y * G)
  }

  # compute bounds
  bounds <- switch(estimand,
                   point = getExtrema(G = G, Y = Y, gamma = log(Lambda), w = w,
                                      estimand = "point", stab = stab),
                   reduction = mean(Y[G == 1]) - rev(getExtrema(G = G, Y = Y, gamma = log(Lambda), w = w,
                                                                estimand = "point", stab = stab)),
                   residual = getExtrema(G = G, Y = Y, gamma = log(Lambda), w = w,
                                         estimand = "point", stab = stab) - mean(Y[G == 0]))

  # return bounds and mu_0_hat
  out <- list(bounds, mu_0_hat, weights_x)

  return(out)
}
