#' @title Estimating effects of incremental propensity score interventions
#'
#' @description \code{ipsi} is used to estimate effects of incremental
#'  propensity score interventions, i.e., estimates of mean outcomes if the odds
#'  of receiving treatment were multiplied by a factor delta.
#'
#' @usage ipsi(dat, x.trt, x.out, delta.seq, nsplits)
#'
#' @param y Outcome of interest measured at end of study.
#' @param a Binary treatment.
#' @param x.trt Covariate matrix for treatment regression.
#' @param x.out Covariate matrix for outcome regression.
#' @param time Measurement time.
#' @param id Subject identifier.
#' @param delta.seq Sequence of delta increment values for incremental
#'  propensity score intervention.
#' @param nsplits Integer number of sample splits for nuisance estimation. If
#'  \code{nsplits = 1}, sample splitting is not used, and nuisance functions are
#'  estimated n full sample (in which case validity of standard errors and
#'  confidence intervals requires empirical process conditions). Otherwise must
#'  have \code{nsplits > 1}.
#' @param ci_level A \code{numeric} value giving the level (1 - alpha) of the
#'  confidence interval to be computed around the point estimate.
#' @param progress_bar A \code{logical} value indicating whether to print a
#'  customized progress bar as various stages of computation reach completion.
#'  The default is \code{TRUE}, printing a progress bar to inform the user.
#' @param return_ifvals A \code{logical} indicating whether the estimated
#'  observation-level values of the influence function ought to be returned as
#'  part of the output object. The default is \code{FALSE} as these values are
#'  rarely of interest in standard usage.
#'
#' @section Details:
#' Treatment and covariates are expected to be time-varying and measured
#' throughout the course of the study. Therefore if \code{n} is the number of
#' subjects and \code{T} the number of timepoints, then \code{a}, \code{time},
#' and \code{id} should all be vectors of length \code{n}x\code{T}, and
#' \code{x.trt} and \code{x.out} should be matrices with \code{n}x\code{T} rows.
#' However \code{y} should be a vector of length \code{n} since it is only
#' measured at the end of the study. The subject ordering should be consistent
#' across function inputs, based on the ordering specified by \code{id}. See
#' example below for an illustration.
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs and uniform CIs for population means.}
#' \item{res.ptwise}{ estimates/SEs and pointwise CIs for population means.}
#' \item{calpha}{ multiplier bootstrap critical value.}
#'
#' @importFrom stats qnorm as.formula
#' @importFrom ranger ranger
#'
#' @examples
#' n <- 500; T <- 4
#'
#' time <- rep(1:T,n); id <- rep(1:n,rep(T,n))
#' x.trt <- matrix(rnorm(n*T*5),nrow=n*T)
#' x.out <- matrix(rnorm(n*T*5),nrow=n*T)
#' a <- rbinom(n*T,1,.5); y <- rnorm(n)
#'
#' d.seq <- seq(0.1,5,length.out=10)
#'
#' ipsi.res <- ipsi(y,a, x.trt,x.out, time,id, d.seq)
#'
#' @references Kennedy EH. Nonparametric causal effects based on incremental
#' propensity score interventions.
#' \href{https://arxiv.org/abs/1704.00211}{arxiv:1704.00211}
#
ipsi <- function(y, a, x.trt, x.out, time, id, delta.seq,
                 nsplits = 2, ci_level = 0.95,
                 progress_bar = TRUE, return_ifvals = FALSE) {
  # setup storage
  ntimes <- length(table(time))
  end <- max(time)
  n <- length(unique(id))
  ynew <- rep(NA, n * ntimes)
  ynew[time == end] <- y
  dat <- data.frame(time = time, id = id, y = ynew, a = a)
  k <- length(delta.seq)
  ifvals <- matrix(nrow = n, ncol = k)
  est.eff <- rep(NA, k)
  wt <- matrix(nrow = n * ntimes, ncol = k)
  cumwt <- matrix(nrow = n*ntimes, ncol = k)
  rt <- matrix(nrow = n * ntimes, ncol = k)
  vt <- matrix(nrow = n * ntimes, ncol = k)
  x.trt <- data.frame(x.trt)
  x.out <- data.frame(x.out)

  if (progress_bar) {
    pb <- txtProgressBar(min = 0, max = 2 * nsplits * length(delta.seq) + 3,
                         style = 3)
  }

  s <- sample(rep(seq_len(nsplits), ceiling(n / nsplits))[seq_len(n)])
  slong <- rep(s, rep(ntimes, n))

  if (progress_bar) {
    pbcount <- 0
  }

  for (split in seq_len(nsplits)) {
    if (progress_bar) {
      Sys.sleep(0.1);
      setTxtProgressBar(pb, pbcount);
      pbcount <- pbcount + 1
    }

    # fit treatment model
    trtmod <- ranger::ranger(stats::as.formula("a ~ ."),
                             dat = cbind(x.trt, a = dat$a)[slong != split, ])
    dat$ps <- predict(trtmod, data = x.trt)$predictions

    for (j in seq_len(k)) {
      if (progress_bar) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, pbcount)
        pbcount <- pbcount + 1
      }
      delta <- delta.seq[j]

      # compute weights
      wt[, j] <- (delta * dat$a + 1 - dat$a) / (delta * dat$ps + 1 - dat$ps)
      cumwt[, j] <- as.numeric(t(aggregate(wt[, j], by = list(dat$id),
                                           cumprod)[, -1]))
      vt[, j] <- (1 - delta) * (dat$a * (1 - dat$ps) -
                                (1 - dat$a) * delta * dat$ps) / delta

      # fit outcome models
      outmod <- vector("list", ntimes)
      rtp1 <- dat$y[dat$time == end]
      if (progress_bar) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, pbcount)
        pbcount <- pbcount + 1
      }
      for (i in seq_len(ntimes)) {
        t <- rev(unique(dat$time))[i]
        outmod[[i]] <-
          ranger::ranger(stats::as.formula("rtp1 ~ ."),
                         dat = cbind(x.out,
                                     rtp1)[dat$time == t & slong != split, ])

        # counterfactual case for treatment: A = 1
        newx1 <- x.out[dat$time == t, ]
        newx1$a <- 1
        m1 <- predict(outmod[[i]], data = newx1)$predictions

        # counterfactual case for no treatment: A = 0
        newx0 <- x.out[dat$time == t, ]
        newx0$a <- 0
        m0 <- predict(outmod[[i]], data = newx0)$predictions

        pi.t <- dat$ps[dat$time == t]
        rtp1 <- (delta * pi.t * m1 + (1 - pi.t) * m0) /
          (delta * pi.t + 1 - pi.t)
        rt[dat$time == t, j] <- rtp1
      }
      # compute influence function values
      ifvals[s == split, j] <- ((cumwt[, j] * dat$y)[dat$time == end] +
        aggregate(cumwt[, j] * vt[, j] * rt[, j],
                  by = list(dat$id), sum)[, -1])[s == split]
    }
  }

  # compute estimator
  for (j in seq_len(k)) {
    est.eff[j] <- mean(ifvals[, j])
  }

  # compute asymptotic variance
  sigma <- sqrt(apply(ifvals, 2, var))
  ci_norm_bounds <- abs(stats::qnorm(p = (1 - ci_level) / 2))
  eff.ll <- est.eff - ci_norm_bounds * sigma / sqrt(n)
  eff.ul <- est.eff + ci_norm_bounds * sigma / sqrt(n)

  # multiplier bootstrap
  if (progress_bar) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, pbcount)
    pbcount <- pbcount + 1
  }

  eff.mat <- matrix(rep(est.eff, n), nrow = n, byrow = TRUE)
  sig.mat <- matrix(rep(sigma, n), nrow = n, byrow = TRUE)
  ifvals2 <- (ifvals - eff.mat) / sig.mat
  nbs <- 10000
  mult <- matrix(2 * rbinom(n * nbs, 1, 0.5) - 1, nrow = n, ncol = nbs)
  maxvals <- sapply(seq_len(nbs), function(col) {
                      max(abs(apply(mult[, col] * ifvals2, 2, sum) / sqrt(n)))
                   })
  calpha <- quantile(maxvals, ci_level)
  eff.ll2 <- est.eff - calpha * sigma / sqrt(n)
  eff.ul2 <- est.eff + calpha * sigma / sqrt(n)

  if (progress_bar) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, pbcount)
    close(pb)
  }

  # conditionally return influence function values
  if (return_ifvals) {
    res <- data.frame(increment = delta.seq, est = est.eff, se = sigma,
                      ci.ll = eff.ll2, ci.ul = eff.ul2, ifvals = ifvals2)
    res2 <- data.frame(increment = delta.seq, est = est.eff, se = sigma,
                       ci.ll = eff.ll, ci.ul = eff.ul, ifvals = ifvals)
  } else {
    res <- data.frame(increment = delta.seq, est = est.eff, se = sigma,
                      ci.ll = eff.ll2, ci.ul = eff.ul2)
    res2 <- data.frame(increment = delta.seq, est = est.eff, se = sigma,
                       ci.ll = eff.ll, ci.ul = eff.ul)
  }

  # output
  return(invisible(list(res = res, res.ptwise = res2, calpha = calpha)))
}
