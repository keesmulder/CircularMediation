library(boot)


computeCircIndirectEffect <- function(data) {

  # Direct effect
  a     <- cg_mod_log1$coefficients[1]

  # Total effect
  a_til <- c_mod_log1$coefficients

  a
  a_til

}


boot()
