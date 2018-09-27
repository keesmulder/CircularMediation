## ----setup, include=FALSE------------------------------------------------


library(boot)
library(ggplot2)

knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
n <- 100
th <- rnorm(n, 1, .5)

cs_cor <- Vectorize(function(rot) cor(sin(th - rot), cos(th - rot)))

# Currently, the cosine and sine are correlated.
cs_cor(0)

# Compute the first eigenvector. 
eigenvec2 <- eigen(cov(cbind(cos(th), sin(th))))$vectors[, 2]

# Direction of the second eigenvector. 
rot <- atan2(eigenvec2[2], eigenvec2[1])

# Now the components are uncorrelated. 
cs_cor(rot)

# There are four solutions.
curve(cs_cor, 0, 2*pi, ylab = "Camel height")
abline(h = 0, col = "skyblue")
abline(v = rot %% (2*pi), col = "red")


th_s <- th - rot
cor(cos(th_s), sin(th_s))

## ------------------------------------------------------------------------
# Calculate the standardized path coefficient between x and y given z.
path_coef <- function(y, x, z) unname(coef(lm(scale(y) ~ scale(x) + z))[2])

## ------------------------------------------------------------------------
# Euclidean norm for later use. 
l2norm <- function(x) sqrt(x[1]^2 + x[2]^2)

# Compute the circular-linear correlation model, with linear predictor x,
# circular outcome th and linear mediator z.
clcor_mediation <- function(th, x, z) {
  
  # Rotation by second eigenvector.
  evec2 <- eigen(cov(cbind(cos(th),  sin(th))))$vectors[, 2]
  rot   <- atan2(evec2[2], evec2[1]) 
  
  # Vector of angles with uncorrelated sine and cosine. 
  th_rot <- th - rot
  
  # First get the results as (cosine, sine) vectors. 
  # Total effect.
  tv <- c(cor(cos(th_rot), x), cor(sin(th_rot), x))

  # Direct effect
  dv <- c(path_coef(cos(th_rot), x, z), path_coef(sin(th_rot), x, z))

  # Indirect effect.
  a <- cor(x, z)
  p_cz <- path_coef(cos(th_rot), z, x)
  p_sz <- path_coef(sin(th_rot), z, x)
  iv <- c(a * p_cz, a * p_sz)
  
  list(effects = c(total    = l2norm(tv), 
                   direct   = l2norm(dv), 
                   indirect = l2norm(iv)), 
       vectors = cbind(total    = tv,
                       direct   = dv,
                       indirect = iv))
}

## ------------------------------------------------------------------------
set.seed(100)

x  <- scale(rnorm(100))
z  <- scale(x + rnorm(100))
th <- 3 + x + z + rnorm(100) %% (2*pi)

cl_med <- clcor_mediation(th, x, z)

## ---- echo = FALSE-------------------------------------------------------
# The results.
knitr::kable(t(cl_med$effects), caption = "Effects in terms of linear-circular correlation.")
knitr::kable(cl_med$vectors, caption = "Effects as vectors of correlation with cosine, sine of outcomee.")

## ------------------------------------------------------------------------
total    <- cl_med$vectors[,1]
direct   <- cl_med$vectors[,2]
indirect <- cl_med$vectors[,3]

# Demonstrate the vector addition formula.
l2norm(total - direct)
l2norm(indirect)

## ------------------------------------------------------------------------

computeClcorMediation <- function(data, inds, outcome, predictor, mediators) {
  listres <- clcor_mediation(th = data[inds, outcome], 
                             z  = data[inds, mediators], 
                             x  = data[inds, predictor])
  res        <- unlist(listres)
  tdi        <- c("Total", "Direct", "Indirect")
  names(res) <- c(tdi, paste0(rep(tdi, each = 2), c("_cos", "_sin")))
  res
}


clcorMediationBootstrap <- function(data, outcome, predictor, mediators, 
                                    R = 1000, probs = c(.025, .975), ...) {
  
  
  # Compute the bootstrap using the boot package.
  suppressWarnings({
    clcormedboot <- boot(data, computeClcorMediation, R = R, 
                         outcome   = outcome, 
                         predictor = predictor, 
                         mediators = mediators, ...)
  })

  # Obtain the p-values for t0 > 0, which will be taken as 1 - p for t0 < 0.
  ps <- apply(clcormedboot$t, 2L, function(x) mean(x <= 0))

  # Obtain result
  res <- list(tab = cbind(original     = clcormedboot$t0,
                         bias         = apply(clcormedboot$t, 2L, mean, na.rm = TRUE) - clcormedboot$t0,
                         "std. error" = sqrt(apply(clcormedboot$t, 2L, function(t.st) var(t.st[!is.na(t.st)]))),
                         t(apply(clcormedboot$t, 2L, function(x) quantile(x, probs = probs))),
                         "one-sided p-value"    = ifelse(clcormedboot$t0 > 0, ps, 1 - ps),
                         "two-sided p-value"    = ifelse(clcormedboot$t0 > 0, ps*2, (1 - ps) * 2 )),
              bootsam = clcormedboot$t,
              bootobj = clcormedboot)

  class(res) <- c("clcorMedBoot", class(res))
  res
}


print.clcorMedBoot <- function(object, ...) print(object$tab, ...)


## ------------------------------------------------------------------------

n <- 1000
x  <- scale(rnorm(n))
z  <- scale(x + rnorm(n))
th <- 3 + x + z + rnorm(n) %% (2*pi)

mydat <- data.frame(th, x, z)

bt <- clcorMediationBootstrap(mydat, "th", "x", "z")

knitr::kable(bt$tab)

## ------------------------------------------------------------------------
bmat <- data.frame(bt$bootsam)
colnames(bmat) <- rownames(bt$tab)

alph <- .4

ggplot(bmat) + 
  geom_point(aes(Total_cos, Total_sin      , col = "Total"   ), alpha = alph) +
  geom_point(aes(Direct_cos, Direct_sin    , col = "Direct"  ), alpha = alph) +
  geom_point(aes(Indirect_cos, Indirect_sin, col = "Indirect"), alpha = alph) +
  ggthemes::theme_fivethirtyeight() + coord_fixed()

## ------------------------------------------------------------------------

ggplot(bmat) + 
  geom_histogram(aes(Total      , fill = "Total"   ), alpha = alph) +
  geom_histogram(aes(Direct     , fill = "Direct"  ), alpha = alph) +
  geom_histogram(aes(Indirect   , fill = "Indirect"), alpha = alph) +
  ggthemes::theme_fivethirtyeight()

