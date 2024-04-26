#' Function to perform amplification of sensitivity parameter
#' @param G Indicator of unmodifiable group or characteristic
#' @param XA Allowable covariates WITHOUT INTERCEPT
#' @param XN Non-allowable covariates WITHOUT INTERCEPT
#' @param Z Treatment
#' @param Y Outcome
#' @param w RMPW weights
#' @param mu_10 Point estimate of the target estimand
#' @param Lambda Sensitivity parameter (The MSM Lambda)
#' @param e1 Propensity score for G = 1
#' @param e0 Propensity score for G = 0
#' @param trim Trimming proportion
#' @param allowable Logical indicating whether to use allowability framework
#'
#' @returns \code{data.frame} with each covariate's imbalance and beta_u term: informal calibration
#'
#' @import dplyr
#'
#' @export


decompAmplify <- function(G, Z, XA, XN, Y, mu_10, Lambda, e1, e0, trim = 0.01, allowable = TRUE, stab = TRUE) {

  bounds <- decompsens::getBiasBounds(G, Z, XA, XN, Y, w, mu10, Lambda = Lam, trim = 0.01, allowable = TRUE)
  maxbias <- max(abs(bounds)) # max{|inf mu_10^h - mu_10|, |sup mu_10^h - mu_10|}
  message("Max bias: ", round(maxbias, 3))
  X <- cbind(XA, XN)
  # standardize X
  X_stnd <- apply(X, MARGIN = 2, FUN = function(x) {(x - mean(x))/sd(x)})

  # standardize X for group G = 1
  X_G1 <- X[G == 1, ] # [, -1] for intercept
  ZG1 <- Z[G == 1]
  X_G1_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {(x - mean(x))/sd(x)})

  ## Compute \beta_u ##
  mod_matrix_y <- data.frame(y = Y[G==1], model.matrix(~ . - 1, data = data.frame(X_G1_stnd))) %>%
    select(-sexM)
  beta1 <- lm(y ~ ., data = mod_matrix_y, weights = ZG1)$coef[-1]
  beta0 <- lm(y ~ ., data = mod_matrix_y, weights = 1 - ZG1)$coef[-1]
  coeffs <- map_dbl(1:length(beta1), function(i) {
    out1 <- beta1[i] * e0
    out0 <- beta0[i] * e1 * (1-e0) / (1-e1)
    name1 <- names(beta1)[i]
    name0 <- names(beta0)[i]
    out <- abs(mean(out1 - out0))
    if (name1 == name0) {
      return(out)
    } else {
      stop("Names of coeffs do not match!")
    }
  }); if (all(names(beta1) == names(beta0))) {names(coeffs) <- names(beta1)} else stop("Names of coeffs do not match!")
  max_betau <- max(coeffs)


  ## Compute standardized imbalance in U: \delta_u ##

  ## Imbalance before weighting
  imbal_stnd <- colMeans(X_G1_stnd) - colMeans(X_G1_stnd[ZG1 == 1, ])
  max_imbal_stnd <- max(abs(imbal_stnd), na.rm = TRUE)

  # Post-weighting imbalance
  wg1 <- w[G == 1]
  Xw_stnd <- apply(X_G1_stnd, MARGIN = 2, FUN = function(x) {x * wg1 / sum(wg1)})

  imbal_stnd_weight <- colMeans(X_G1_stnd) - colSums(Xw_stnd[ZG1 == 1, ]) # sum is reweighted
  max_imbal_stnd_wt <- max(abs(imbal_stnd_weight), na.rm = TRUE)

  # Get coordinates for strongest observed covariates to plot
  coeff_df <- data.frame(
    covar = gsub("X_stnd[G == 1, ]", "", names(coeffs)),
    coeff = abs(as.numeric(coeffs)))

  imbal_df <- data.frame(
    covar = names(imbal_stnd),
    imbal = abs(as.numeric(imbal_stnd)),
    imbal_wt = abs(as.numeric(imbal_stnd_weight)))

  # merge coefficients and imbalance, arrange from largest to smallest by imbal * beta_u
  strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, by = "covar") %>%
    dplyr::arrange(desc(coeff*imbal))

  return(list(df = strongest_cov_df, max_beta = max_betau, max_imbal_stnd = max_imbal_stnd,
              max_imbal_stnd_wt = max_imbal_stnd_wt, maxbias = maxbias))
}
