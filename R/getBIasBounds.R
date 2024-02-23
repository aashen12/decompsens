#' Function to obtain bias bounds for amplification
#' @export

getBiasBounds <- function(G, Z, XA, XN, Y, Lambda, trim = 0.05,
                          allowable = FALSE, stab = TRUE) {
  ####################################################################
  # Get bounds on bias = \mu^*_10 - \hat{mu}_10
  # lower bound = inf h \hat{mu}_10^{(h)} - \hat{mu}_10
  # upper bound = inf h \hat{mu}_10^{(h)} - \hat{mu}_10
  # there is no reduction or residual since it cancels out
  #
  ####################################################################

  mu1 <- mean(Y[G == 1])
  mu0 <- mean(Y[G == 0])
  w <- estimateRMPW(G = G, Z = Z, Y = Y, XA = XA, XN = XN,
                         trim = trim, allowable = allowable)
  if (stab) {
    mu_10 <- sum(w * Y * G) / sum(w * G)
  } else {
    mu_10 <- mean(w * Y * G)
  }

  # compute bounds
  bounds <- getExtrema(G = G, Y = Y, gamma = log(Lambda), w = w,
                       estimand = "point", stab = stab) - mu_10
  # even if we care about red or res, we use point bc the lambda already takes
  # into account whether we care about red or res

  # return bounds and mu_0_hat
  out <- list(bounds, mu_10, w[G == 1])

  return(out)
}
