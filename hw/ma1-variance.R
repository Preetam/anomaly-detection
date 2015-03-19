### Variance of EWMA smoothed MA(1) ###

# Import code
source("hw.R")

# Testing helper
assertWithin <- function(x, expected, bound) {
  withinBound <-  abs(x - expected) < bound
  print(sprintf("expecting %0.3g, got %0.3g (off by %0.3g). %s.",
                expected, x, abs(x-expected),
                ifelse(withinBound, "OK", "FAIL")))
  stopifnot( withinBound )
}

#
# s[t] = (1 - α) s[t-1] + α x[t]
#      = (1 - α)^t s[0] + α sum( (1 - α)^(t-i) x[i], {i, 1, t})
#

alpha <- 0.1
t <- 10

x <- arima.sim(list(ma=c(0.1)), n=t+1)

x0 <- x[1]
x <- x[-1]

hw.obj <- hw.new(alpha = alpha, beta = 0)
hw.obj <- hw.step(hw.obj, x0)
for (i in 1:t) {
  hw.obj <- hw.step(hw.obj, x[i])
}

# This is the iterative result.
s_10 <- hw.obj$level.value

# Calculating the closed-form result.
s_10_cf <- 0
for (i in 1:t) {
  s_10_cf <- s_10_cf + (1 - alpha)^(t - i) * x[i]
}
s_10_cf <- s_10_cf * alpha
s_10_cf <- s_10_cf + (1-alpha)^t * x0

c(s_10, s_10_cf)
assertWithin(s_10_cf, s_10, 0.001)

####
alpha <- 0.25
sigma <- 1 # standard deviation of MA(1)
theta <- 0.1
t <- 10 # variance of t-th value.
n <- 20000 # number of iterations

s_t <- rep(0, n)
x_t <- rep(0, n)
for (iteration in 1:n) {
  x <- arima.sim(list(ma=c(theta), start.innov=c(0)), n=t)
  x0 <- 0
  hw.obj <- hw.new(alpha = alpha, beta = 0)
  hw.obj <- hw.step(hw.obj, x0)
  for (i in 1:t) {
    hw.obj <- hw.step(hw.obj, x[i])
  }
  
  s_t[iteration] <- hw.obj$level.value
  x_t[iteration] <- x[t]
}

sum_1 <- 0
for (i in 1:t) {
  sum_1 <- sum_1 + alpha^2 * (1 - alpha)^(2*t - 2*i)
}

sum_2 <- 0
if (t > 1) {
  for (i in (1:(t-1))) {
    print(i)
    sum_2 <- sum_2 + (1 - alpha)^(2*t - (2*i+1))
  }
}

var_s_t <- sigma^2 * (theta^2 + 1) * sum_1 + 2 * alpha^2 * sigma^2 * theta * sum_2

assertWithin(var_s_t, var(s_t), 0.01)
