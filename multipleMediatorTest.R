x <- scale(rnorm(100))
y <- scale(x + rnorm(100))
cor(x, y)
lm(y ~ x)
lm(x ~ y)

rnddig <- 6

x <- scale(rnorm(100))
z1 <- scale(x + rnorm(100))
z2 <- scale(x + rnorm(100))
y <- scale(x/2 + z1/2 + z2/2 + rnorm(100))

coef_med <- round(c(z1 = coef(lm(z1~x))[2], z2 = coef(lm(z2~x))[2]), rnddig)
coef_y   <- round(coef(lm(y ~ x + z1 + z2))[-1], rnddig)
coef_tot <- round(coef(lm(y ~ x))[-1], rnddig)

# Total
total <- coef_tot

# Direct
direct <- coef_y[1]

# Indirect effect separated
indirect <- coef_med * coef_y[-1] # (a_1 * b_1, a_2 * b_2)

# Indirect effect combined, two ways!
total - direct
sum(indirect)


