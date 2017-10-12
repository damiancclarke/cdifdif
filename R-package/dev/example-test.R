library(broom)
library(magrittr)

data <- spilloverDGP
formula <- y1 ~ time + treat
dist <- spilloverDGP$dist
maxDist <- 30
delta <- 1
alpha <- 0.05
k <- 10
verbose <- TRUE

step <- 5
