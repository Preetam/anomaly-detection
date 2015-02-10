###
# Holt-Winters
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

hw.new <- function(alpha = 0.1, beta = 0.1) {
  hw.obj <- list()
  hw.obj$alpha <- alpha
  hw.obj$beta <- beta
  hw.obj$level <- 0
  hw.obj$level.value <- NULL
  hw.obj$trend.value <- NULL
  hw.obj$error <- NULL
  
  hw.obj
}

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
    hw.obj$error <- (hw.obj$level.value - new.value)^2
    return(hw.obj)
  }
  
  # Calculate new smoothed level value
  s <- hw.obj$alpha * new.value + (1 - hw.obj$alpha) * (hw.obj$level.value + hw.obj$trend.value)
  
  # Calculate new smoothed trend value
  b <- hw.obj$beta * (s - hw.obj$level.value) + (1 - hw.obj$beta) * hw.obj$trend.value
  
  # Calculate error
  e <- (hw.obj$level.value - new.value)^2
  hw.obj$error <- hw.obj$alpha * e + (1 - hw.obj$alpha) * hw.obj$error
  
  hw.obj$level.value <- s
  hw.obj$trend.value <- b
  
  hw.obj
}

###
# Test
###
# sim <- arima.sim(list(ar=c(0.2)), n=2000)
# hw.obj <- hw.step(hw.new(alpha=0.001), 0)
# smoothed <- rep(NA, 2000)
# low <- rep(NA, 2000)
# high <- rep(NA, 2000)
# for (i in seq(1,2000)) {
#   hw.obj <- hw.step(hw.obj, sim[i])
#   smoothed[i] <- hw.obj$level.value
#   print(hw.obj$error)
#   if (is.null(hw.obj$error) == F) {
#     low[i] <- smoothed[i] - hw.obj$error * 3
#     high[i] <- smoothed[i] + hw.obj$error * 3
#   }
# }
# plot.ts(sim)
# lines(smoothed, lwd=5)
# lines(low, lwd=4, col='blue')
# lines(high, lwd=4, col='blue')
###