
#' @title Estimating effects of incremental propensity score interventions
#'
#' @description \code{ipsi} is used to estimate effects of incremental propensity score interventions, i.e., estimates of mean outcomes if the odds of receiving treatment were multiplied by a factor delta.
#'
#' @usage ipsi(dat, x.trt, x.out, delta.seq, nsplits)
#'
#' @param y outcome of interest measured at end of study.
#' @param a binary treatment.
#' @param x.trt covariate matrix for treatment regression.
#' @param x.out covariate matrix for outcome regression.
#' @param time measurement time.
#' @param id subject identifier.
#' @param delta.seq sequence of delta increment values.
#' @param nsplits integer number of sample splits for nuisance estimation.
#' If nsplits=1, sample splitting is not used, and nuisance functions are estimated
#' on full sample (in which case validity of SEs/CIs requires empirical
#' process conditions). Otherwise must have nsplits>1.
#'
#' @section Details:
#' Treatment and covariates are expected to be time-varying and measured throughout the course of the study. Therefore if \code{n} is the number of subjects and \code{T} the number of timepoints, then \code{a}, \code{time}, and \code{id} should all be vectors of length \code{n}x\code{T}, and \code{x.trt} and \code{x.out} should be matrices with \code{n}x\code{T} rows. However \code{y} should be a vector of length \code{n} since it is only measured at the end of the study. The subject ordering should be consistent across function inputs, based on the ordering specified by \code{id}. See example below for an illustration.
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs and uniform CIs for population means.}
#' \item{res.ptwise}{ estimates/SEs and pointwise CIs for population means.}
#' \item{calpha}{ multiplier bootstrap critical value.}
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
#' @references Kennedy EH. Nonparametric causal effects based on incremental propensity score interventions. \href{https://arxiv.org/abs/1704.00211}{arxiv:1704.00211}
#'
ipsi <- function(y,a, x.trt, x.out, time, id, delta.seq, nsplits=2){
require(ranger)

  # setup storage
  ntimes <- length(table(time)); end <- max(time); n <- length(unique(id))
  ynew <- rep(NA,n*ntimes); ynew[time==end] <- y
  dat <- data.frame(time=time,id=id,y=ynew,a=a)
  k <- length(delta.seq); ifvals <- matrix(nrow=n,ncol=k); est.eff <- rep(NA,k)
  wt <- matrix(nrow=n*ntimes,ncol=k); cumwt <- matrix(nrow=n*ntimes,ncol=k)
  rt <- matrix(nrow=n*ntimes,ncol=k); vt <- matrix(nrow=n*ntimes,ncol=k)
  x.trt <- data.frame(x.trt); x.out <- data.frame(x.out)

  pb <- txtProgressBar(min=0, max=2*nsplits*length(delta.seq)+3, style=3)

  s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
  slong <- rep(s,rep(ntimes,n))

  pbcount <- 0
  for (split in 1:nsplits){ Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

    # fit treatment model
    trtmod <- ranger(a ~ ., dat=cbind(x.trt,a=dat$a)[slong!=split,])
    dat$ps <- predict(trtmod, data=x.trt)$predictions

    for (j in 1:k){ Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
      delta <- delta.seq[j]

      # compute weights
      wt[,j] <- (delta*dat$a + 1-dat$a)/(delta*dat$ps + 1-dat$ps)
      cumwt[,j] <- as.numeric(t(aggregate(wt[,j],by=list(dat$id),cumprod)[,-1]))
      vt[,j] <- (delta-1)*(dat$a*(1-dat$ps) - (1-dat$a)*delta*dat$ps)/delta

      # fit outcome models
      outmod <- vector("list",ntimes); rtp1 <- dat$y[dat$time==end]
      Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
      for (i in 1:ntimes){
        t <- rev(unique(dat$time))[i]
        outmod[[i]] <- ranger(rtp1 ~ ., dat=cbind(x.out,rtp1)[dat$time==t & slong!=split,])
        newx1 <- x.out[dat$time==t,]; newx1$a <- 1
        m1 <- predict(outmod[[i]], data=newx1)$predictions
        newx0 <- x.out[dat$time==t,]; newx0$a <- 0
        m0 <- predict(outmod[[i]], data=newx0)$predictions
        pi.t <- dat$ps[dat$time==t]
        rtp1 <- (delta*pi.t*m1 + (1-pi.t)*m0) / (delta*pi.t + 1-pi.t)
        rt[dat$time==t,j] <- rtp1 }

      ifvals[s==split,j] <- ((cumwt[,j]*dat$y)[dat$time==end] +
        aggregate(cumwt[,j]*vt[,j]*rt[,j],by=list(dat$id),sum)[,-1])[s==split]

    } }

  # compute estimator
  for (j in 1:k){ est.eff[j] <- mean(ifvals[,j]) }

  # compute asymptotic variance
  sigma <- sqrt(apply(ifvals,2,var))
  eff.ll <- est.eff-1.96*sigma/sqrt(n); eff.ul <- est.eff+1.96*sigma/sqrt(n)

  # multiplier bootstrap
  Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
  eff.mat <- matrix(rep(est.eff,n),nrow=n,byrow=T)
  sig.mat <- matrix(rep(sigma,n),nrow=n,byrow=T)
  ifvals2 <- (ifvals-eff.mat)/sig.mat
  nbs <- 10000; mult <- matrix(2*rbinom(n*nbs,1,.5)-1,nrow=n,ncol=nbs)
  maxvals <- sapply(1:nbs, function(col){
    max(abs(apply(mult[,col]*ifvals2,2,sum)/sqrt(n))) } )
  calpha <- quantile(maxvals, 0.95)
  eff.ll2 <- est.eff-calpha*sigma/sqrt(n); eff.ul2 <- est.eff+calpha*sigma/sqrt(n)

  Sys.sleep(0.1)
  setTxtProgressBar(pb,pbcount)
  close(pb)

  res <- data.frame(increment=delta.seq, est=est.eff, se=sigma, ci.ll=eff.ll2, ci.ul=eff.ul2)
  res2 <- data.frame(increment=delta.seq, est=est.eff, se=sigma, ci.ll=eff.ll, ci.ul=eff.ul)

  return(invisible(list(res=res, res.ptwise=res2, calpha=calpha)))

}
