### Variance of EWMA smoothed MA(n) ###

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

sample.ma.covariance <- function(thetas, t, n) {
  x1 <- rep(0, n + 1)
  xt <- rep(0, n + 1)
  
  for (i in 1:n) {
    x <- arima.sim(list(ma=thetas, start.innov=c(0)), n=t+1)
    x1[i] <- x[1]
    xt[i] <- x[t+1]
  }
  
  cov(x1,xt)
}

smoothed.ma.variance <- function(alpha, thetas, t) {
  sum_1 <- 0
  for (i in 1:t) {
    sum_1 = sum_1 + ( alpha^2 * (1 - alpha)^(2*t - 2*i) )
  }
  
  sum_2 <- 0
  for (h in 1:(length(thetas))) {
    cov.lag.h <- sample.ma.covariance(thetas, h, 20000)
    
    sum_3 <- 0
    if (t > 1) {
      for (i in 1 : (t - h)) {
        sum_3 = sum_3 + ( alpha^2 * (1 - alpha)^(2 * t - (2*i + h)) )
      }
    }
    
    sum_2 = sum_2 + cov.lag.h * sum_3
  }
  
  sample.ma.covariance(thetas, 0, 20000) * sum_1 + 2 * sum_2
}

alpha <- 0.25
t <- 10 # lag
n <- 10000 # number of iterations
thetas <- c(0.3, 0.2, 0.1, 0.2, 0.3)
s_t <- rep(0, n)
x_t <- rep(0, n)

for (iteration in 1:n) {
  x <- arima.sim(list(ma=thetas, start.innov=c(0)), n=t)
  x0 <- 0
  hw.obj <- hw.new(alpha = alpha, beta = 0)
  hw.obj <- hw.step(hw.obj, x0)
  for (i in 1:t) {
    hw.obj <- hw.step(hw.obj, x[i])
  }
  
  s_t[iteration] <- hw.obj$level.value
  x_t[iteration] <- x[t]
}

var(s_t)
smoothed.ma.variance(alpha, thetas, t)
