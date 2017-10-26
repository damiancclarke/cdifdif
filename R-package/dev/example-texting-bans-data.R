# #delimit ;
# cdifdif `yvar' `cont2' handheldban `wt', distance(closeDistanceHH) `se'
# h(11) regtype(areg) tlimit(1.64) maxdist(100) ;
#
# eststo: cdifdif `yvar' handheldban `cont2' `wt', distance(closeDistanceHH) `se'
# maxdist(`mdist') delta(`delta') regtype(areg) optimal crossfold kfold(10)
#
# eststo: cdifdif `yvar' strongban `cont2' `wt', distance(closeDistanceSt) `se'
# maxdist(`mdist') delta(`delta') regtype(areg) optimal crossfold kfold(10)
# tlimit(1.64) verbose stub(strongclose);
# graph export "$OUT/RMSE_strong.eps", replace;
# local strongdist = e(distmax);
#
# eststo: cdifdif `yvar' weakban `cont2' `wt', distance(closeDistanceWeak) `se'
#                         maxdist(`mdist') delta(`delta') regtype(areg) optimal crossfold kfold(10)
# tlimit(1.64) verbose stub(weakclose);
rm(list = ls())
library(haven)
library(tidyverse)
library(broom)
library(viridis)


data <- read_dta("../example/estdata_munic.dta.zip")

saveRDS(data, "data-raw/estdata_munic.rds")

# remove state dummy variables
data <- data %>%
  select(-matches("st\\d{1,2}"), -matches("t\\d{1,2}"), -matches("mon\\d{1,2}"))

mdist <- 40
delta <- 2

select(data, laccidentsvso2, strongban, weakban)

data$closeDistanceAll

data %>% count(treatedt) %>% tail()
data %>% count(treated)
data %>% count(time)


dat <- data %>%
  select(
    county1, time, month, year, laccidentsvso2,
    starts_with("treat"), contains("Distance"), ends_with("ban")
    ) %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x))

dat

mod <- cdifdif(laccidentsvso2 ~ time + treated, data = dat, dist = data$closeDistanceStrong,
               maxDist = mdist, delta = delta, alpha = 0.05, k = 20, verbose = TRUE)


plot(mod$cvs, type = "l", main = "CV RMSEs")

nmod <- which.min(mod$cvs)
nmod

mod$cvs[[nmod]]

tidy(mod$mods[[nmod]])

step <- mod$steps[[nmod]]
step

cuts <- c(seq(0, 5), Inf) * step
cuts


