#' Procedure to obtain bootstrap CI of RMPW estimator in causal decomposition
#' problems
#'
#' @param G Indicator of unmodifiable group or characteristic
#' @param Y Outcome
#' @param gamma Sensitivity parameter (log of MSM Lambda)
#' @param w fitted RMPW weights
#' @param alpha Significance level
#' @param estimand Target estimand, either "point", "reduction" or "residual"
#' @param parallel Logical indicating whether to use parallel processing
#' @param B Number of bootstrap samples
#'
#'
#' @return Extrema (an interval) of the RMPW estimator.
#'
#' @import parallel
#'
#' @export


boostrapCI <- function(G, Y, gamma = 0, w, alpha = 0.05, estimand = "point",
                       parallel = TRUE, B = 1000) {
  estimand <- match.arg(estimand, c("point", "reduction", "residual"))
  # assume observed rmpw weights already computed
  mu1 <- mean(Y[G == 1])
  mu0 <- mean(Y[G == 0])
  mu_10 <- sum(w * Y * G) / sum(w * G)
  # message("mu1: ", round(mu1, 2))
  # message("mu0: ", round(mu0, 2))

  state <- switch(estimand, point = round(mu_10, 3),
                  reduction = round(mu1 - mu_10, 3),
                  residual = round(mu_10 - mu0, 3))

  message("Parameter of interest: ", estimand)
  message("Estimate: ", state)

  no.cores <- max(1, ifelse(parallel, detectCores(), 1))
  n <- length(G)

  # if (warm.start) {
  #   start <- glm(A ~ X, family = "binomial")$coefs
  # } else {
  #   start <- NULL
  # }

  out <- mclapply(1:B, function(iter) {
    s <- sample(1:n, n, TRUE);
    res <- tryCatch(getExtrema(G[s], Y[s], gamma, w[s], estimand),
                    error = function(e) {print(e)});
    res
  }, mc.cores = no.cores)

  out <- do.call(rbind, out)

  c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))
}
