### Variance of EWMA smoothed MA(2) ###

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

####
alpha <- 0.25
sigma <- 1 # standard deviation of MA(1)
theta1 <- 0.1
theta2 <- 0.1
t <- 50 # variance of t-th value.
n <- 10000 # number of iterations

s_t <- rep(0, n)
x_t <- rep(0, n)
for (iteration in 1:n) {
  x <- arima.sim(list(ma=c(theta1, theta2), start.innov=c(0)), n=t)
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
    sum_2 <- sum_2 + alpha^2 * (1 - alpha)^(2*t - (2*i+1))
  }
}

sum_3 <- 0
if (t > 2) {
  for (i in (1:(t-2))) {
    sum_3 <- sum_3 + alpha^2 * (1 - alpha)^(2*t - (2*i+2))
  }
}

var_s_t <- sigma^2 * (1 + theta1^2 + theta2^2) * sum_1 + 2 * sigma^2 * (theta1 + theta1*theta2) * sum_2 + 2 * sigma^2 * theta2 * sum_3

assertWithin(var_s_t, var(s_t), 0.01)

