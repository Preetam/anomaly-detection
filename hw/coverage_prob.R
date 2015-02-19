### Coverage probabilities ###

# Import code
source("hw.R")

ar.bounds <- function(arma.obj) {
  list(low = arma.obj$predicted.value - 1.959964,
       high = arma.obj$predicted.value + 1.959964)
}

hw.bounds <- function(hw.obj) {
  list(low = hw.obj$level.value - sqrt(hw.obj$error) * 1.959964,
       high = hw.obj$level.value + sqrt(hw.obj$error) * 1.959964)
}

assertWithin <- function(x, expected, bound) {
  print(sprintf("expecting %0.3g, got %0.3g (off by %0.3g)", expected, x, abs(x-expected)))
  stopifnot( abs(x - expected) < bound )
}

####

exceeded <- 0
iterations <- 10000
for (iteration in 1:iterations) {
  sim <- arima.sim(list(ar=c(0.1)), n=11)
  
  arma.obj <- arma.new(0.1)
  for (i in 1:10) {
    arma.obj <- arma.step(arma.obj, sim[i])
  }
  
  bounds <- ar.bounds(arma.obj)
  if (isTRUE(sim[11] > bounds$high | sim[11] < bounds$low)) {
    exceeded <- exceeded + 1
  }
}

assertWithin(exceeded/iterations, 0.05, 0.01)

exceeded <- 0
iterations <- 10000
for (iteration in 1:iterations) {
  sim <- arima.sim(list(ar=c(0.1)), n=201)
  
  hw.obj <- hw.new(alpha=0.05, beta=0.05)
  for (i in 1:200) {
    hw.obj <- hw.step(hw.obj, sim[i])
  }
  
  bounds <- hw.bounds(hw.obj)
  if (isTRUE(sim[201] > bounds$high | sim[201] < bounds$low)) {
    exceeded <- exceeded + 1
  }
}

assertWithin(exceeded/iterations, 0.05, 0.01)

