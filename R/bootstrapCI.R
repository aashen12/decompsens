#' Procedure to obtain bootstrap CI of RMPW estimator in causal decomposition
#' problems
#'
#' @param G Indicator of unmodifiable group or characteristic
#' @param XA Allowable covariates
#' @param XN Non-allowable covariates
#' @param Z Treatment
#' @param Y Outcome
#' @param gamma Sensitivity parameter (log of MSM Lambda)
#' @param w fitted RMPW weights
#' @param alpha Significance level
#' @param estimand Target estimand, either "point", "reduction" or "residual"
#' @param parallel Logical indicating whether to use parallel processing
#' @param B Number of bootstrap samples
#' @param trim Trimming proportion
#' @param allowable Logical indicating whether to use allowability framework
#'
#' @return Bootstrap CI of the desired estimator.
#'
#' @import parallel
#'
#' @export


bootstrapCI <- function(G, Z, Y, XA, XN, gamma = 0, w, alpha = 0.05, estimand = "point",
                       parallel = TRUE, B = 1000, stab = TRUE,
                       allowable = FALSE, trim = 0.05) {
  estimand <- match.arg(estimand, c("point", "reduction", "residual"))

#   state <- switch(estimand, point = round(mu_10, 3),
#                   reduction = round(mu1 - mu_10, 3),
#                   residual = round(mu_10 - mu0, 3))

  message("Parameter of interest: ", estimand)
  #message("Estimate: ", state)

  no.cores <- max(1, ifelse(parallel, detectCores(), 1))
  n <- length(G)

  # if (warm.start) {
  #   start <- glm(Z ~ X, family = "binomial")$coefs
  # } else {
  #   start <- NULL
  # }

  out <- mclapply(1:B, function(iter) {
    s <- sample(1:n, n, TRUE);
    Ys <- Y[s]
    Gs <- G[s]
    XAs <- XA[s,]
    XNs <- XN[s,]
    Zs <- Z[s]

    # re-estimate weights
    w_rmpw_boot <- estimateRMPW(G = Gs, Z = Zs, Y = Ys, XA = XAs, XN = XNs,
                                trim = trim, allowable = allowable)

    mu10_B <- tryCatch(getExtrema(G[s], Y[s], gamma, w_rmpw_boot, estimand = "point", stab),
                    error = function(e) {print(e)});
    switch(estimand,
           point = mu10_B,
           reduction = mean(Ys[Gs == 1]) - rev(mu10_B),
           residual = mu10_B - mean(Ys[Gs == 0]))
  }, mc.cores = no.cores)

  out <- do.call(rbind, out)

  c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))
}
