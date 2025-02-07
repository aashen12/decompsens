# R Package for performing sensitivity analysis for causal decompositions

## Installation

`devtools::install_github("https://github.com/aashen12/decompsens")`

`devtools::install_github("https://github.com/aashen12/decompsens", ref = "main", auth_token = Sys.getenv("GITHUBTOKEN"))`

## Functions

- `bootstrapCI.R`: uses percentile bootstrap to create confidence intervals for disparity reduction/residual under unmeasured confounding (Equation 18, Theorem 1)
- `estimateRMPW.R`: estimates RMPW weights used to compute the weighted average in Equation 7
- `getExtrema.R` and `getBiasBounds.R`: computes extrema of point estimate (Equations 17)
- `decompAmplify.R` and `informalAmplify.R`: perform amplification
