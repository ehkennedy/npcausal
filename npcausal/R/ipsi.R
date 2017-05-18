#' @title Estimating effects of incremental propensity score interventions
#'
#' @description \code{ipsi} is used to estimate effects of incremental propensity score interventions.
#'
#' @usage ipsi(dat, x.trt, x.out, delta.seq, nsplits)
#'
#' @param dat dataframe (in long not wide form if longitudinal) with columns ‘time’, ‘id’, outcome ‘y’, treatment ‘a’.
#' @param x.trt covariate matrix for treatment regression.
#' @param x.out covariate matrix for outcome regression.
#' @param delta.seq sequence of delta values.
#' @param nsplits number of sample splits.
#'
#' @return The sum of \code{x} and \code{y}.
#'
#' @references Kennedy EH. Nonparametric causal effects based on incremental propensity score interventions. \href{https://arxiv.org/abs/1704.00211}{arxiv:1704.00211}
#'
ipsi <- function(dat, x.trt, x.out, delta.seq, nsplits){
  # setup storage
  ntimes <- length(table(dat$time)); n <- length(unique(dat$id))
  k <- length(delta.seq); ifvals <- matrix(nrow=n,ncol=k); est.eff <- rep(NA,k)
  wt <- matrix(nrow=n*ntimes,ncol=k); cumwt <- matrix(nrow=n*ntimes,ncol=k)
  rt <- matrix(nrow=n*ntimes,ncol=k); vt <- matrix(nrow=n*ntimes,ncol=k)
  s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
  slong <- rep(s,rep(ntimes,n))
  for (split in 1:nsplits){ print(paste("split",split)); flush.console()
    # fit treatment model
    trtmod <- ranger(a ~ ., dat=cbind(x.trt,a=dat$a)[slong!=split,])
    dat$ps <- predict(trtmod, data=x.trt)$predictions
    for (j in 1:k){ print(paste("delta",j)); flush.console()
      delta <- delta.seq[j]
      # compute weights
      wt[,j] <- (delta*dat$a + 1-dat$a)/(delta*dat$ps + 1-dat$ps)
      cumwt[,j] <- as.numeric(t(aggregate(wt[,j],by=list(dat$id),cumprod)[,-1]))
      vt[,j] <- (delta-1)*(dat$a*(1-dat$ps) - (1-dat$a)*delta*dat$ps)/delta
      # fit outcome models
      outmod <- vector("list",ntimes); rtp1 <- dat$y[dat$time==end]
      print("fitting regressions"); flush.console()
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
  eff.mat <- matrix(rep(est.eff,n),nrow=n,byrow=T)
  sig.mat <- matrix(rep(sigma,n),nrow=n,byrow=T)
  ifvals2 <- (ifvals-eff.mat)/sig.mat
  nbs <- 10000; mult <- matrix(2*rbinom(n*nbs,1,.5)-1,nrow=n,ncol=nbs)
  maxvals <- sapply(1:nbs, function(col){
    max(abs(apply(mult[,col]*ifvals2,2,sum)/sqrt(n))) } )
  calpha <- quantile(maxvals, 0.95)
  eff.ll2 <- est.eff-calpha*sigma/sqrt(n); eff.ul2 <- est.eff+calpha*sigma/sqrt(n)
  return(list(est=est.eff, sigma=sigma, ll1=eff.ll,ul1=eff.ul,
              calpha=calpha, ll2=eff.ll2,ul2=eff.ul2))
}
