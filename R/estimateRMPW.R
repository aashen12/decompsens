#' Helper function to estimate RMPW weights
#' @param G Indicator of unmodifiable group or characteristic
#' @param Z Treatment
#' @param Y Outcome
#' @param XA Allowable covariates WITHOUT INTERCEPT
#' @param XN Non-allowable covariates WITHOUT INTERCEPT
#' @param trim0 Trimming proportion for e_0
#' @param trim1 Trimming proportion for e_1
#' @param allowable Logical indicating whether to use allowability framework
#'
#' @return w_rmpw RMPW weights, e_g group propensity scores
#'
#' @export


estimateRMPW <- function(G, Z, Y, XA, XN, trim0 = 0.05, trim1 = 0.05, allowable = TRUE) {
  X <- cbind(XA, XN)
  # model matrix
  if (allowable) {
    e0 <- glm(Z ~ ., data = XA, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
    e1 <- glm(Z ~ ., data = X, family = binomial, weights = G, na.action = na.exclude)$fitted.values
  } else {
    e0 <- glm(Z ~ ., data = X, family = binomial, weights = 1 - G, na.action = na.exclude)$fitted.values
    e1 <- glm(Z ~ ., data = X, family = binomial, weights = G, na.action = na.exclude)$fitted.values
  }

  e0 <- pmax(pmin(e0, 1 - trim0), trim0)
  e1 <- pmax(pmin(e1, 1 - trim1), trim1)

  w1 = e0 / e1
  w0 = (1 - e0) / (1 - e1)
  w_rmpw = w1 * Z + w0 * (1 - Z)

  return(list(w_rmpw = w_rmpw, e0 = e0, e1 = e1))
}
