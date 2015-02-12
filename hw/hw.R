###
# Holt-Winters and AR(1) Forecasting
###
# Copyright (c) 2015, Preetam Jinka
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   
#   * Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#          SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###

# hw.new creates a new hw.obj.
hw.new <- function(alpha = 0.1, beta = 0.1) {
  hw.obj <- list()
  hw.obj$alpha <- alpha
  hw.obj$beta <- beta
  hw.obj$level <- 0
  hw.obj$level.value <- NULL
  hw.obj$trend.value <- NULL
  hw.obj$error <- 1

  # return the new hw.obj
  hw.obj
}

# hw.step updates the given hw.obj and returns
# a new hw.obj with an updated state.
hw.step <- function(hw.obj, new.value) {
  if (is.null(hw.obj)) {
    stop("hw.obj is NULL")
  }

  if (is.null(hw.obj$level.value)) {
    # First value
    hw.obj$level.value <- new.value
    return(hw.obj)
  }

  if (is.null(hw.obj$trend.value)) {
    # Second value
    hw.obj$trend.value <- 0
    return(hw.obj)
  }

  # Calculate new smoothed level value
  s <- hw.obj$alpha * new.value + (1 - hw.obj$alpha) * (hw.obj$level.value + hw.obj$trend.value)

  # Calculate new smoothed trend value
  b <- hw.obj$beta * (s - hw.obj$level.value) + (1 - hw.obj$beta) * hw.obj$trend.value

  # Calculate error
  e <- (hw.obj$level.value - new.value)^2
  hw.obj$error <- hw.obj$alpha * e + (1 - hw.obj$alpha) * hw.obj$error

  # Update level and trend components
  hw.obj$level.value <- s
  hw.obj$trend.value <- b

  hw.obj
}

# arma.new creates a new arma.obj.
arma.new <- function(phi) {
  arma.obj <- list()
  arma.obj$phi <- phi
  arma.obj$predicted.value <- NULL

  arma.obj
}

# arma.step updates the given arma.obj and returns
# a new arma.obj with an updated state.
arma.step <- function(arma.obj, new.value) {
  arma.obj$predicted.value <- arma.obj$phi * new.value

  arma.obj
}

###
# Test
###

plot.simulation <- function(N = 200, phi = 0.1, alpha = 0.01, beta = 0) {
  # Simulate AR(1) with parameter phi for N points.
  sim <- arima.sim(list(ar=c(phi)), n=N)

  arma.obj <- arma.new(phi)
  hw.obj <- hw.new(alpha=alpha, beta=beta)

  ###
  # Vectors to hold time series points.
  ###

  # EWMA series
  smoothed <- rep(NA, N)
  low <- rep(NA, N)
  high <- rep(NA, N)

  # ARMA series
  arma.prediction <- rep(NA, N)
  arma.low <- rep(NA, N)
  arma.high <- rep(NA, N)

  # Iteration through N steps
  for (i in seq(1,N)) {
    hw.obj <- hw.step(hw.obj, sim[i])
    arma.obj <- arma.step(arma.obj, sim[i])

    smoothed[i] <- hw.obj$level.value
    arma.prediction[i] <- arma.obj$predicted.value

    if (is.null(hw.obj$error) == F) {
      low[i] <- smoothed[i] - hw.obj$error * 1.96
      high[i] <- smoothed[i] + hw.obj$error * 1.96

      arma.low[i] <- arma.prediction[i] - 1.96
      arma.high[i] <- arma.prediction[i] + 1.96
    }
  }

  ###
  # Plotting
  ###
  plot.ts(sim, lwd=2, ylim=c(-5,5), ylab='Value',
          main='AR(1) with Exponential Smoothing',
          sub=paste('phi = ', phi, ", alpha = ", alpha, ", beta = ", beta, sep=''))
  lines(smoothed,        lwd=3, lty=1, col='purple1')
  lines(arma.prediction, lwd=2, lty=6, col='chartreuse3')
  lines(low,             lwd=.5, lty=1, col='purple1')
  lines(high,            lwd=.5, lty=1, col='purple1')
  lines(arma.low,        lwd=1, lty=3, col='chartreuse4')
  lines(arma.high,       lwd=1, lty=3, col='chartreuse4')

  legend('topright', c("AR series", "EWMA prediction", "AR prediction", "EWMA PI", "AR PI"),
         col=c('black', 'purple1', 'chartreuse3', 'purple1', 'chartreuse4'),
         lwd=c(2,3,2,.5,1),
         lty=c(1,1,3,1,3))
}

plot.simulation()
plot.simulation(phi=0.5)
plot.simulation(phi=0.2, alpha=0.1)
plot.simulation(phi=0.2, alpha=0.01, beta=0.2)
plot.simulation(phi=0.2, alpha=0.01, beta=0.01, N=300)
