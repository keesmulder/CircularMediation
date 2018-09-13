
# Calculate the partial correlation between x and y given z.
part_cor <- function(x, y, z) cor(residuals(lm(x ~ z)), residuals(lm(y ~ z)))

# Get circular linear association between a circular th and linear x.
rho <- function(th, x) {
  r_cx <- cor(cos(th), x)
  r_sx <- cor(sin(th), x)
  r_cs <- cor(cos(th), sin(th))
  sqrt( (r_cx^2 + r_sx^2 - 2 * r_cx * r_sx * r_cs) / (1 - r_cs^2) )
}

# Get partial circular linear association between a circular th and linear x.
rho_partial <- function(th, x, z) {
  r_cx <- part_cor(cos(th), x, z)
  r_sx <- part_cor(sin(th), x, z)
  r_cs <- part_cor(cos(th), sin(th), z)
  sqrt( (r_cx^2 + r_sx^2 - 2 * r_cx * r_sx * r_cs) / (1 - r_cs^2) )
}

x  <- scale(rnorm(100))
z  <- scale(x + rnorm(100))
th <- 3 + x + z + rnorm(100) %% (2*pi)


# Total, direct, indirect effects:

total <- rho(th, x)
direct <- rho_partial(th, x, z)

#### Indirect

# product
a <- cor(x, z)
b <- rho_partial(th, z, x)
a * b

total - direct
