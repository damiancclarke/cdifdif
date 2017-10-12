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
#' @param k a number of folds for k-fold Cross-Validation.
#' @param verbose Verbose
#'
#' @export
cdifdif <- function(y , X, dist,
                    maxDist = 30, #quantile(dist[dist!=0], 0.75),
                    delta   = 1,  #quantile(dist[dist!=0], 0.025),
                    alpha = 0.05,
                    k = 5,
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

  cvs1 <- mods %>%
    lapply(function(x) x[["mod"]]) %>%
    lapply(cvTools::repCV, K = k)

  cvs2 <- mods %>%
    lapply(function(x) x[["mod"]]) %>%
    lapply(get_cv_rmse_from_mod_k, k = k)


  microbenchmark::microbenchmark(
    cvs1 <- mods %>%
      lapply(function(x) x[["mod"]]) %>%
      lapply(cvTools::repCV, K = k)
    ,
    cvs2 <- mods %>%
      lapply(function(x) x[["mod"]]) %>%
      lapply(get_cv_rmse_from_mod_k, k = k)
    , times = 5
  )


  cvs1v <- lapply(cvs1, function(x){ x[["cv"]]}) %>%
    unlist() %>%
    as.vector()
  cvs2v <- unlist(cvs2)

  cvs1v
  cvs2v


}


# library(dplyr)
# library(purrr)
# library(broom)
#
# data("spilloverDGP")
# spilloverDGP
# data    <- spilloverDGP %>% select(y = y1, time, treat)
# maxDist <- 30
# k       <- 5
# delta   <- 1
# verbose <- TRUE
# dist <- spilloverDGP$dist

marginal_dist <- function(data, dist, step, alpha = 0.05, verbose = TRUE) {

  dl    <- 0
  t     <- 1
  pval  <- 0

  daux  <- data

  mod_aux  <- lm(y ~ ., data = daux)
  mod_lst <- list()

  while(pval < alpha) {

    if(verbose) message("iteration: ", t, " step: ", step)

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

get_cv_rmse_from_mod_k <- function(mod, k) {

  d <- mod$model

  folds <- createFolds(seq_len(nrow(mod$model)), k = k)

  ss_errs <- lapply(folds, function(f) {

    # print(f)

    mod <- lm(y ~ ., data = d[-f,])
    yht <- predict(mod, newdata = d[f, ])

    err <- d[f, "y"] - yht

    ss_err <-  sum(err^2)

    ss_err

  })

  rmse <- sqrt( sum(unlist(ss_errs)) / nrow(d) )

  rmse

}
