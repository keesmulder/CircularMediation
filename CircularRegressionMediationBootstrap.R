library(boot)
library(circular)

# Function to compute
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

# Convenience function to only obtain one on the results provided.


boot(data, computeCircMediation, R = 10)
