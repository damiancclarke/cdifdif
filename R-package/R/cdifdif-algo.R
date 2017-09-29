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
cdifdif <- function(y , X, dist,maxDist, delta, tlimit, CVtype, kfolds) {

  TRUE

}

