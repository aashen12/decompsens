#' Helper function to estimate RMPW weights
#' @param G Indicator of unmodifiable group or characteristic
#' @param Z Treatment
#' @param Y Outcome
#' @param XA Allowable covariates WITHOUT INTERCEPT
#' @param XN Non-allowable covariates WITHOUT INTERCEPT
#' @param trim Trimming proportion
#' @param allowable Logical indicating whether to use allowability framework
#'
#' @return RMPW weights
#'
#' @export


estimateRMPW <- function(G, Z, Y, XA, XN, trim = 0.05, allowable = FALSE) {
  X <- cbind(XA, XN)
  if (allowable) {
    e0 <- glm(Z ~ XA, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
    e1 <- glm(Z ~ X, family = binomial, weights = G, na.action = na.exclude)$fitted.values
  } else {
    e0 <- glm(Z ~ X, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
    e1 <- glm(Z ~ X, family = binomial, weights = G, na.action = na.exclude)$fitted.values
  }

  e0 <- pmax(pmin(e0, 1 - trim), trim)
  e1 <- pmax(pmin(e1, 1 - trim), trim)

  w1 = e0 / e1
  w0 = (1 - e0) / (1 - e1)
  w_rmpw = w1 * Z + w0 * (1 - Z)

  return(w_rmpw)
}
