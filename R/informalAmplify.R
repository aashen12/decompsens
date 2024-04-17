#' Function to perform amplification of sensitivity parameter
#' @param G Indicator of unmodifiable group or characteristic
#' @param XA Allowable covariates WITHOUT INTERCEPT
#' @param XN Non-allowable covariates WITHOUT INTERCEPT
#' @param Z Treatment
#' @param Y Outcome
#' @param w RMPW weights
#' @param mu_10 Point estimate of the target estimand
#' @param Lambda Sensitivity parameter (The MSM Lambda)
#' @param trim Trimming proportion
#' @param allowable Logical indicating whether to use allowability framework
#'
#' @returns \code{data.frame} with each covariate's imbalance and beta_u term: informal calibration
#'
#' @import dplyr
#'
#' @export


informalAmplify <- function(G, Z, XA, XN, Y, mu_10, Lambda, trim = 0.01, allowable = TRUE, stab = TRUE) {

  bounds <- getExtrema(G=G, Y=Y, w=w, gamma = log(Lam), estimand = "point", RD = TRUE, verbose = FALSE)
  maxbias <- max(abs(bounds)) # max{|inf mu_10^h - mu_10|, |sup mu_10^h - mu_10|}
  message("Max bias: ", round(maxbias, 3))
  X <- cbind(XA, XN)
  # standardize X
  X_stnd <- apply(X, MARGIN = 2, FUN = function(x) {(x - mean(x))/sd(x)})

  # standardize X for group G = 1
  X_G1 <- X[G == 1, ] # [, -1] for intercept
  X_G1_stnd <- apply(X_G1, MARGIN = 2, FUN = function(x) {(x - mean(x))/sd(x)})

  ## Compute \beta_u ##
  mod_matrix_y <- data.frame(y = Y[G==1], model.matrix(~ . - 1, data = data.frame(X_G1_stnd))) %>% select(-sexM)
  coeffs <- lm(y ~ ., data = mod_matrix_y)$coef[-1]
  max_betau <- max(abs(coeffs), na.rm = TRUE)

  ## Compute standardized imbalance in U: \delta_u ##

  ZG1 <- Z[G == 1]

  ## Imbalance before weighting
  imbal_stnd <- colMeans(X_G1_stnd[ZG1 == 1, ]) - colMeans(X_G1_stnd[ZG1 == 0, ])
  max_imbal_stnd <- max(abs(imbal_stnd), na.rm = TRUE)

  # Post-weighting imbalance
  wg1 <- w[G == 1]
  Xw_stnd <- apply(X_G1_stnd, MARGIN = 2, FUN = function(x) {x * wg1 / sum(wg1)})

  imbal_stnd_weight <- colSums(Xw_stnd[ZG1 == 1, ]) - colSums(Xw_stnd[ZG1 == 0, ]) # sum is reweighted
  max_imbal_stnd_wt <- max(abs(imbal_stnd_weight), na.rm = TRUE)

  # Get coordinates for strongest observed covariates to plot
  coeff_df <- data.frame(
    covar = gsub("X_stnd[G == 1, ]", "", names(coeffs)),
    coeff = abs(as.numeric(coeffs)))

  # get imbal
  imbal_df <- data.frame(
    covar = names(imbal_stnd),
    imbal = abs(as.numeric(imbal_stnd)),
    imbal_wt = abs(as.numeric(imbal_stnd_weight)))

  # merge coefficients and imbalance, arrange from largest to smallest by imbal * beta_u
  strongest_cov_df <- dplyr::inner_join(coeff_df, imbal_df, by = "covar") %>%
    dplyr::arrange(desc(coeff*imbal))

  return(list(df = strongest_cov_df, max_beta = max_betau, max_imbal_stnd = max_imbal_stnd,
               max_imbal_stnd_wt = max_imbal_stnd_wt))
}
