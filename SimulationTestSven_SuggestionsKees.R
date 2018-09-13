library(circglmbayes)

#data simulation
n <- 1000
a <- 0.5 # Relation x to m
b <- 0.5 # Relation m to y
c <- 0 # Relation x to y (direct only)

x <- rnorm(n, 0, 5)
m <- rnorm(n, a*x, 1)
beta <- c(c, b)
pred <- cbind(x,m)
linkfun <- function(x) 2 * atan(x)

y_pred <- linkfun(pred %*% as.matrix(beta))

err <- rvmc(n, 0, 15)
y <- y_pred + err
y <- as.circular(y)

plot(x, y)

plot(x, y %% (2*pi))

plot(data.frame(x, m, y))
