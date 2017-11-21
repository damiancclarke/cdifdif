#' Spillover Robust Diff-in-Diff Estimation
#'
#' Implements Spillover Robust Diff-in-Diff Estimates (Clarke, 2017)(Clarke, 2017)
#' @param formula formula.
#' @param data data.
#' @param dist a N-by-1 distance to treatment.
#' @param maxDist a maximum spillover bandwidth to consider.
#' @param delta a step-size for bandwidth search (based on dist variable).
#' @param alpha a minimum t-stat to consider marginal spillover to be significant.
#' @param weights (only for weighted fits) the specified weights.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param k a number of folds for k-fold Cross-Validation.
#' @param verbose Verbose
#'
#' @examples
#'
#' data("spilloverDGP")
#' mod <- cdifdif(y1 ~ time + treat, data = spilloverDGP, dist = spilloverDGP$dist,
#'                maxDist = 30, delta = 1, alpha = 0.05, k = 10, verbose = FALSE)
#' mod
#'
#'
#' @importFrom stats model.frame quantile gaussian
#' @export
cdifdif <- function(formula, data, dist,
                    maxDist = 30, #quantile(dist[dist!=0], 0.75),
                    delta   = 1,  #quantile(dist[dist!=0], 0.025),
                    alpha = 0.05,
                    weights = NULL,
                    family = gaussian,
                    k = 1000,
                    verbose = TRUE) {

  data <- model.frame(formula, data)
  names(data) <- c("y", "t1", "t2")

  # head(data)

  # linspace(delta,maxDist,round(maxDist/delta))
  # http://mathesaurus.sourceforge.net/octave-r.html
  steps <- seq(from = delta, to = maxDist, length = round(maxDist/delta))

  mods <- lapply(steps, marginal_dist, data = data, dist = dist, alpha = alpha,
                 verbose = verbose, weights = weights, family = family)

  finalmods <- lapply(mods, function(x) x[["mod"]])

  cvs <- lapply(finalmods, get_cv_rmse_from_mod_k, k = k, verbose = verbose)
  cvs <- unlist(cvs)
  # plot(cvs)

  list(
    mods = finalmods,
    cvs = cvs,
    steps = steps
  )

}

#' @importFrom stats lm setNames glm
#' @importFrom broom tidy
#' @importFrom magrittr %>%
marginal_dist <- function(data, dist, step, alpha = 0.05, verbose = TRUE,
                          weights = NULL, family = gaussian) {

  dl    <- 0
  t     <- 1
  pval  <- 0

  daux  <- data

  mod_aux  <- lm(y ~ ., data = daux)
  mod_lst <- list()

  while(pval < alpha) {

    if(verbose) message("iteration: ", t, " step: ", step)

    dnew <- as.numeric(dl < dist & dist <= dl + step)
    daux <- cbind(daux, setNames(data.frame(dnew), paste0("d", t)))

    mod_old <- mod_aux

    # https://stats.stackexchange.com/a/181180/59417
    mod_aux <- glm(y ~ ., data = daux, weights = weights, family = family)

    modsummary <- tidy(mod_aux) # summary(mod_aux)

    if(!paste0("d", t) %in% modsummary[["term"]]) break

    mod_lst <- c(mod_lst, list(modsummary))

    # updates
    dl   <- dl + step
    pval <- abs(modsummary[nrow(modsummary), "p.value"])
    t    <- t + 1

  }

  list(mod = mod_old, summaries = mod_lst)

}

#' @importFrom stats predict
get_cv_rmse_from_mod_k <- function(mod, k, verbose = verbose) {

  d <- mod$model

  folds <- createFolds(seq_len(nrow(mod$model)), k = k)

  s_errs <- lapply(seq_along(folds), function(fold) {

    if(verbose) message("fold: ", fold)

    idx <- folds[[fold]]

    mod <- lm(y ~ ., data = d[-idx,])
    yht <- predict(mod, newdata = d[idx, ])

    err <- d[idx, "y"] - yht

    s_err <-  err^2

    s_err

  })

  s_errs <- unlist(s_errs)

  rmse <- sqrt(mean(s_errs))

  rmse

}
