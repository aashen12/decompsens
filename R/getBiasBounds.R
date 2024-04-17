#' Function to obtain bias bounds for amplification
#' @param G Indicator of unmodifiable group or characteristic
#' @param XA Allowable covariates
#' @param XN Non-allowable covariates
#' @param Z Treatment
#' @param Y Outcome
#' @param w RMPW weights
#' @param mu_10 Point estimate of the target estimand
#' @param Lambda Sensitivity parameter (The MSM Lambda)
#' @param trim Trimming proportion
#' @param allowable Logical indicating whether to use allowability framework
#'
#' @returns A list of 3 items: the lower and upper bounds: mu_10(h) - mu_10,
#' the point estimate mu_10, the RMPW weights
#'
#'
#' @export

getBiasBounds <- function(G, Z, XA, XN, Y, w, mu_10, Lambda, trim = 0.05,
                          allowable = TRUE, stab = TRUE) {
  ####################################################################
  # Get bounds on bias = \mu^*_10 - \hat{mu}_10
  # lower bound = inf h \hat{mu}_10^{(h)} - \hat{mu}_10
  # upper bound = inf h \hat{mu}_10^{(h)} - \hat{mu}_10
  # there is no reduction or residual since it cancels out
  # compute wrt point estimate
  ####################################################################

  # mu1 <- mean(Y[G == 1])
  # mu0 <- mean(Y[G == 0])
  # w <- estimateRMPW(G = G, Z = Z, Y = Y, XA = XA, XN = XN,
  #                        trim = trim, allowable = allowable)
  # if (stab) {
  #   mu_10 <- sum(w * Y * G) / sum(w * G)
  # } else {
  #   mu_10 <- mean(w * Y * G)
  # }

  # compute bounds
  bounds <- getExtrema(G = G, Y = Y, gamma = log(Lambda), w = w,
                       estimand = "point", stab = stab) - mu_10
  # even if we care about red or res, we use point bc the lambda already takes
  # into account whether we care about red or res

  # return bounds and mu_10_hat

  return(bounds)
}
