data <-

#' Spillover Robust Diff-in-Diff Estimation
#'
#' Implements Spillover Robust Diff-in-Diff Estimates (Clarke, 2017)(Clarke, 2017)
#' @param y a N-by-1 dependent variable.
#' @param X a N-by-k baseline independent variables
#' @param dist a N-by-1 distance to treatment.
#' @param maxDist a maximum spillover bandwidth to consider.
#' @param delta a step-size for bandwidth search (based on dist variable).
#' @param tlimit a minimum t-stat to consider marginal spillover to be significant.
#' @param CVtype a type of Cross-Validation (must be either 'kfoldcv' or 'loocv').
#' @param kfolds a number of folds for k-fold Cross-Validation.
#'
#' @export
cdifdif <- function(y , X, dist,
                    maxDist = quantile(dist[dist!=0], 0.75),
                    delta   = quantile(dist[dist!=0], 0.025),
                    alpha = 0.05,
                    CVtype, kfolds,
                    verbose = TRUE) {

  ff <- log(Volume) ~ log(Height) + log(Girth)

  mat <- model.matrix(y ~ ., data)
  mat <- model.matrix(y ~ . -1, data)
  head(mat)

  # linspace(delta,maxDist,round(maxDist/delta))
  # http://mathesaurus.sourceforge.net/octave-r.html
  steps <- seq(from = delta, to = maxDist, length = round(maxDist/delta))

  # mods <- purrr::map(steps, marginal_dist, data = data, dist = dist, alpha = alpha, verbose = verbose)
  mods <- lapply(steps, marginal_dist, data = data, dist = dist, alpha = alpha, verbose = verbose)

  cvs <- mods %>%
    lapply(function(x) x[["mod"]]) %>%
    lapply(cvTools::repCV, K = 1000)

  cvs[[1]]
  lapply(cvs, function(x){ x[["cv"]]}) %>%
    unlist() %>%
    plot(type = "l", ylim = c(0, max(.)*1.1))

  cvs <- mods %>%
    purrr::map(1) %>%
    purrr::map(cvTools::repCV, K = 1000)



}

data("spilloverDGP")
spilloverDGP
data <- spilloverDGP %>% select(y = y2, time, treat)
step <- 5
step <- 52.375
dist <- spilloverDGP$dist

library(dplyr)
library(purrr)
library(broom)
library(lmvar)

marginal_dist <- function(data, dist, step, alpha = 0.05, verbose = TRUE) {

  dl    <- 0
  t     <- 1
  pval  <- 0

  daux  <- data

  mod_aux  <- lm(y ~ ., data = daux)
  mod_lst <- list()

  while(pval < alpha) {

    if(verbose) message("step: ", step, " iteration: ", t)

    dnew <- as.numeric(dl < dist & dist <= dl + step)
    daux <- dplyr::bind_cols(daux, purrr::set_names(data_frame(dnew), paste0("d", t)))

    mod_old <- mod_aux
    mod_aux <- lm(y ~ ., data = daux)

    modsummary <- broom::tidy(mod_aux) # summary(mod_aux)

    if(!paste0("d", t) %in% modsummary[["term"]]) break

    mod_lst <- c(mod_lst, list(modsummary))

    # updates
    dl   <- dl + step
    pval <- abs(modsummary[nrow(modsummary), "p.value"])
    t    <- t + 1

  }

  list(mod = mod_old, summaries = mod_lst)

}

library(boot)

k_fold_rsq <- function(lmfit, ngroup=10) {
  # assumes library(bootstrap)
  # adapted from http://www.statmethods.net/stats/regression.html
  mydata <- lmfit$model
  outcome <- names(lmfit$model)[1]
  predictors <- names(lmfit$model)[-1]

  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef}
  X <- as.matrix(mydata[predictors])
  y <- as.matrix(mydata[outcome])

  results <- crossval(X,y,theta.fit,theta.predict,ngroup=ngroup)
  raw_rsq <- cor(y, lmfit$fitted.values)**2 # raw R2
  cv_rsq <- cor(y,results$cv.fit)**2 # cross-validated R2

  c(raw_rsq=raw_rsq, cv_rsq=cv_rsq)
}
