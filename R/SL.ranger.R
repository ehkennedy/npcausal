#' @title Add Ranger wrapper for SuperLearner
#'
#' @description \code{SL.ranger} is a wrapper for \code{SuperLearner} that adds the fast random forests method \code{ranger}.
#'
#' @usage SL.ranger(Y, X, newX, family, ...)
#'
#' @param Y outcome vector.
#' @param X covariate dataframe for training.
#' @param newX covariate dataframe for predictions.
#' @param family link function (currently only supports "gaussian" identity link).
#'
#' @return Predictions and fits from \code{ranger}.
#'
#' @references Wright MN, Ziegler A (2016). ranger: A fast implementation of random forests for high dimensional data in C++ and R. \emph{Journal of Statistical Software}.
#'
SL.ranger <- function (Y, X, newX, family, ...) {
  require("ranger")
  fit.rf <- ranger::ranger(Y ~ ., data=X)
  pred <- predict(fit.rf,data=newX)$predictions
  fit <- list(object = fit.rf)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.ranger")
  return(out) }
