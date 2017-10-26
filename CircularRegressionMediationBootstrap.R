library(boot)
library(circular)

# Function to compute a single circular mediation effect.
# data is a data.frame with named columns.
# outcome is a string giving the name of the outcome
# predictor is a string giving the name of the predictor
# mediators is a string vector giving the name of the mediator
computeCircMediation <- function(data, inds, outcome, predictor, mediators) {

  # This iteration's dataset.
  td <- data[inds, ]

  # Run the two required circular regression models.
  total_model    <- circular:::lm.circular.cl(y    = as.matrix(td[, outcome]),
                                              x    = as.matrix(td[, predictor]),
                                              init = rep(0, length(predictor)))

  mediator_model <- circular:::lm.circular.cl(y    = as.matrix(td[, outcome]),
                                              x    = as.matrix(td[, c(predictor, mediators)]),
                                              init = rep(0, length(predictor) + length(mediators)))



  # Direct effect
  a_tilde <- total_model$coefficients[1]

  # Total effect
  a       <- mediator_model$coefficients[1]

  # Return the total, direct, and indirect effect.
  c(total_effect    = a_tilde,
    direct_effect   = a,
    indirect_effect = a_tilde - a)
}

# Function to compute a single circular mediation effect.
# data is a data.frame with named columns.
# outcome is a string giving the name of the outcome
# predictor is a string giving the name of the predictor
# mediators is a string vector giving the name of the mediator
# R is the number of bootstrap iterations
# probs is numeric vector giving the proportions at which quantiles should be given
circMediationBootstrap <- function(data, outcome, predictor, mediators, R = 1000, probs = c(.025, .975), ...) {

  # Compute the bootstrap using the boot package.
  circmedboot <- boot(data, computeCircMediation, R = R, outcome = outcome, predictor = predictor, mediators = mediators, ...)

  # Obtain the p-values for t0 > 0, which will be taken as 1 - p for t0 < 0.
  ps <- apply(circmedboot$t, 2L, function(x) mean(x <= 0))

  # Obtain result
  cbind(original     = circmedboot$t0,
        bias         = apply(circmedboot$t, 2L, mean, na.rm = TRUE) - circmedboot$t0,
        "std. error" = sqrt(apply(circmedboot$t, 2L, function(t.st) var(t.st[!is.na(t.st)]))),
        t(apply(circmedboot$t, 2L, function(x) quantile(x, probs = probs))),
        "p-value"    = ifelse(circmedboot$t0 > 0, ps, 1 - ps))
}

