#' @title Doubly robust series estimation of counterfactual densities
#'
#' @description \code{cdensity} is used to estimate counterfactual densities,
#'  i.e., the density of the potential outcome in a population if everyone
#'  received given treatment levels, using doubly robust estimates of L2
#'  projections of the density onto a linear basis expansion. Nuisance functions
#'  are estimated with random forests. The L2 distance between the density of the
#'  counterfactuals is also estimated as a density-based treatment effect.
#'
#' @usage cdensity(y, a, x, kmax=5, l2 = TRUE,
#'  gridlen=20, nsplits=2, progress_updates = TRUE,
#'  makeplot=TRUE, kforplot=5, ylim=NULL)
#'
#' @param y outcome of interest.
#' @param a binary treatment (more than 2 levels are allowed, but only densities under
#'  A=1 and A=0 will be estimated).
#' @param x covariate matrix.
#' @param kmax Integer indicating maximum dimension of (cosine) basis
#'  expansion that should be used in series estimator.
#' @param l2 A \code{logical} value indicating whether an estimate of the L2 distance
#'  between counterfactual densities (under A=1 vs A=0) should be returned.
#' @param gridlen Integer number indicating length of grid for which the
#'  plug-in estimator of the marginal density is computed.
#' @param nsplits Integer number of sample splits for nuisance estimation. If
#'  \code{nsplits = 1}, sample splitting is not used, and nuisance functions are
#'  estimated n full sample (in which case validity of standard errors and
#'  confidence intervals requires empirical process conditions). Otherwise must
#'  have \code{nsplits > 1}.
#' @param progress_updates A \code{logical} value indicating whether to print a
#'  progress statement as various stages of computation reach completion.
#'  The default is \code{TRUE}, printing a progress bar to inform the user.
#' @param makeplot A \code{logical} value indicating whether to print a plot.
#' @param kforplot A vector of two integers indicating which k values to plot results for,
#'  with first argument for A=1 and second for A=0.
#' @param ylim Range of y values at which density should be plotted.
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for population means and relevant contrasts.}
#' \item{nuis}{ subject-specific estimates of nuisance functions (i.e., propensity score and outcome regression) }
#' \item{ifvals}{ matrix of estimated influence function values.}
#'
#' @importFrom stats qnorm as.formula
#' @importFrom ranger ranger
#'
#' @export
#'
#' @examples
#' n <- 100; x <- matrix(rnorm(n*5),nrow=n)
#' a <- sample(3,n,replace=TRUE)-2; y <- rnorm(n)
#'
#' cdens.res <- cdensity(y,a,x)
#'
#' @references Kennedy EH, Wasserman LA, Balakrishnan S. Semiparametric counterfactual
#' density estimation. \href{https://arxiv.org/abs/TBA}{arxiv:TBA}
#

cdensity <- function(y, a, x,
               kmax=5, l2=TRUE,
               gridlen=20, nsplits=2, progress_updates = TRUE,
               makeplot=TRUE, kforplot=c(5,5), ylim=NULL) {

  require("ranger")
  pos.part <- function(x){ x*(x>0)+0*(x<0) }

  ### preliminaries

  # rescale outcome in [0,1]
  n <- length(y)
  ymax <- max(y,na.rm=T); ymin <- min(y,na.rm=T)
  yorig <- y; y <- (yorig-ymin)/(ymax-ymin)
  # set treatment indicators
  a1 <- as.numeric(a==1); a0 <- as.numeric(a==0)
  if (sum(a1==0 & a0==0)>0){ y[a1==0 & a0==0] <- 0 }
  # set grid for plugin estimator
  ygrid <- quantile(y[a1==1 | a0==1],probs=seq(0,1,length.out=gridlen))

  # create splits & put into superlearner form
  v <- nsplits
  s <- sample(rep(1:v, ceiling(n/v))[1:n])
  splits <- vector(mode="list", length=v); names(splits) <- 1:v
  for (i in 1:v){ splits[[i]] <- (1:n)[s==i] }

  # initialize storage vectors
  pi1hat <- rep(NA,n); pi0hat <- rep(NA,n)
  mu1mat <- matrix(nrow=n,ncol=kmax)
  mu0mat <- matrix(nrow=n,ncol=kmax)
  eta1mat <- matrix(nrow=n,ncol=gridlen)
  eta0mat <- matrix(nrow=n,ncol=gridlen)
  p1hat <- rep(NA,gridlen); p0hat <- rep(NA,gridlen)
  p1yrega1hat <- rep(NA,n); p0yrega1hat <- rep(NA,n)
  p1yrega0hat <- rep(NA,n); p0yrega0hat <- rep(NA,n)

  # create basis functions & compute at obs y vals
  #bfun <- function(t,k){ (k%%2 == 1) *sqrt(2)*cos(2*ceiling(k/2)*pi*t) + (k%%2 == 0) *sqrt(2)*sin(2*ceiling(k/2)*pi*t) }
  bfun <- function(t,k){ sqrt(2)*cos(k*pi*t) }
  b <- Vectorize(bfun,vectorize.args="k")
  bmat <- NULL; for (k in 1:kmax){ bmat <- cbind(bmat, b(y,k)) }

  # loop through folds
  for (test in 1:v){
  if (progress_updates==T){ print(paste("fold:",test)); flush.console() }

  ### estimate nuisance functions

  # estimate propensity scores
  if (progress_updates==T){ print("    estimating propensity scores..."); flush.console() }
  pi1mod <- ranger(a1[s!=test] ~ ., data=data.frame(x[s!=test,]))
  pi0mod <- ranger(a0[s!=test] ~ ., data=data.frame(x[s!=test,]))
  pi1hat[s==test] <- predict(pi1mod, data=data.frame(x[s==test,]))$predictions
  pi0hat[s==test] <- predict(pi0mod, data=data.frame(x[s==test,]))$predictions

  # estimate basis-transformed outcome regressions
  if (progress_updates==T){ print("    estimating outcome regressions..."); flush.console() }
  for (i in 1:kmax){
    mu1mod <- ranger(bmat[a1==1 & s!=test,i] ~ ., data=data.frame(x[a1==1 & s!=test,]))
    mu0mod <- ranger(bmat[a0==1 & s!=test,i] ~ ., data=data.frame(x[a0==1 & s!=test,]))
    mu1mat[s==test,i] <- predict(mu1mod, data=data.frame(x[s==test,]))$predictions
    mu0mat[s==test,i] <- predict(mu0mod, data=data.frame(x[s==test,]))$predictions
  }

  # estimate conditional density
  if (progress_updates==T){ print("    estimating conditional densities..."); flush.console() }
  for (j in 1:gridlen){
    # regress kernel outcome w/silverman's rule on x w/ranger
    h <- bw.nrd0(y[a1==1 | a0==1])
    kern <- dnorm((y-ygrid[j])/h)/h
    eta1mod <- ranger(kern[a1==1 & s!=test] ~ ., data=data.frame(x[a1==1 & s!=test,]))
    eta0mod <- ranger(kern[a0==1 & s!=test] ~ ., data=data.frame(x[a0==1 & s!=test,]))
    eta1mat[s==test,j] <- predict(eta1mod, data=data.frame(x[s==test,]))$predictions
    eta0mat[s==test,j] <- predict(eta0mod, data=data.frame(x[s==test,]))$predictions
    # marginalize in training data for density regressions
    p1hat[j] <- mean(predict(eta1mod, data=data.frame(x[s!=test,]))$predictions)
    p0hat[j] <- mean(predict(eta0mod, data=data.frame(x[s!=test,]))$predictions)
  }

  # estimate regressions of p1y,p0y on x
  if (progress_updates==T){ print("    estimating density regressions..."); flush.console() }
  p1y.trn <- approx(ygrid,p1hat,xout=y,rule=2)$y
  p0y.trn <- approx(ygrid,p0hat,xout=y,rule=2)$y
  p1ymod1 <- ranger(p1y.trn[a1==1 & s!=test] ~ ., data=data.frame(x[a1==1 & s!=test,]))
  p0ymod1 <- ranger(p0y.trn[a1==1 & s!=test] ~ ., data=data.frame(x[a1==1 & s!=test,]))
  p1ymod0 <- ranger(p1y.trn[a0==1 & s!=test] ~ ., data=data.frame(x[a0==1 & s!=test,]))
  p0ymod0 <- ranger(p0y.trn[a0==1 & s!=test] ~ ., data=data.frame(x[a0==1 & s!=test,]))
  p1yrega1hat[s==test] <- predict(p1ymod1, data=data.frame(x[s==test,]))$predictions
  p0yrega1hat[s==test] <- predict(p0ymod1, data=data.frame(x[s==test,]))$predictions
  p1yrega0hat[s==test] <- predict(p1ymod0, data=data.frame(x[s==test,]))$predictions
  p0yrega0hat[s==test] <- predict(p0ymod0, data=data.frame(x[s==test,]))$predictions

  }

  ### plug-in estimator

  p1hatvals <- apply(eta1mat,2,mean)
  p0hatvals <- apply(eta0mat,2,mean)
  p1hatfn <- function(t){ approx(ygrid,p1hatvals,xout=t,rule=2)$y }
  p0hatfn <- function(t){ approx(ygrid,p0hatvals,xout=t,rule=2)$y }
  p1hat.const <- integrate(p1hatfn, lower=0,upper=1)$val
  p0hat.const <- integrate(p0hatfn, lower=0,upper=1)$val
  p1hat <- function(t){ approx(ygrid,p1hatvals,xout=t,rule=2)$y/p1hat.const }
  p0hat <- function(t){ approx(ygrid,p0hatvals,xout=t,rule=2)$y/p0hat.const }

  ###  l2 projection

  drmean1 <- apply((a1/pi1hat)*(bmat-mu1mat) + mu1mat,2,mean)
  drmean0 <- apply((a0/pi0hat)*(bmat-mu0mat) + mu0mat,2,mean)
  g1fun <- function(k,t){ pos.part(1 + b(t,1:k) %*% drmean1[1:k]) }
  g0fun <- function(k,t){ pos.part(1 + b(t,1:k) %*% drmean0[1:k]) }
  g1kconst <- function(k){ integrate(function(t){ g1fun(k,t) },lower=0,upper=1)$val }
  g0kconst <- function(k){ integrate(function(t){ g0fun(k,t) },lower=0,upper=1)$val }
  g1 <- function(k,t){ g1fun(k,t) / g1kconst(k) }
  g0 <- function(k,t){ g0fun(k,t) / g0kconst(k) }

  ### L2 distance

  if (l2==TRUE){
  l2.plugin <- integrate( function(t){ (p1hat(t)-p0hat(t))^2 }, lower=0,upper=1)$value

  p1yhat <- p1hat(y); p0yhat <- p0hat(y)
  psif.if <- 2 * ( (a1/pi1hat)*((p1yhat-p0yhat) - (p1yrega1hat-p0yrega1hat)) + (p1yrega1hat-p0yrega1hat) -
    ( (a0/pi0hat)*((p1yhat-p0yhat) - (p1yrega0hat-p0yrega0hat)) + (p1yrega0hat-p0yrega0hat) ) - l2.plugin)

  est <- l2.plugin + mean(psif.if)
  se <- sqrt(var(psif.if)/n)
  ci.ll <- est - 1.96*se; ci.ul <- est + 1.96*se

  ### results

  res <- data.frame(parameter="L2 distance", est,se,ci.ll,ci.ul) }

  ifvals <- list(a1 = (a1/pi1hat)*(bmat-mu1mat) + mu1mat, a0= (a0/pi0hat)*(bmat-mu0mat) + mu0mat)

  if (makeplot==T){
  se1 <- function(k,t){ sqrt(diag( b(t,1:k) %*% cov(ifvals$a1)[1:k,1:k] %*% t(b(t,1:k)) )/n) }
  se0 <- function(k,t){ sqrt(diag( b(t,1:k) %*% cov(ifvals$a0)[1:k,1:k] %*% t(b(t,1:k)) )/n) }

  par(mfrow=c(1,1))
  if (is.null(ylim)){ yrng <- c(0,1) }
  if (!is.null(ylim)){ yrng <- (ylim-ymin)/(ymax-ymin) }
  yseq <- seq(yrng[1],yrng[2],length.out=250)*(ymax-ymin)+ymin
  g1est <- g1(kforplot[1],(yseq-ymin)/(ymax-ymin))/(ymax-ymin); g0est <- g0(kforplot[2],(yseq-ymin)/(ymax-ymin))/(ymax-ymin)
  se1est <- se1(kforplot[1],(yseq-ymin)/(ymax-ymin))/(ymax-ymin); se0est <- se0(kforplot[2],(yseq-ymin)/(ymax-ymin))/(ymax-ymin)

  plot(yseq, g1est, type="l", lwd=2,col="darkblue",
     xlab="y",ylab="Estimate", main="Counterfactual density (& pointwise CI)",
     ylim=range(c(0,g0est+1.96*se0est,g1est+1.96*se1est)))
  polygon(c(yseq, rev(yseq)), c(g0est-1.96*se0est, rev(g0est+1.96*se0est)), col="lightsalmon",border=NA)
  polygon(c(yseq, rev(yseq)), c(g1est-1.96*se1est, rev(g1est+1.96*se1est)), col="lightblue3",border=NA)
  lines(yseq, g1est, lwd=2,col="darkblue")
  lines(yseq, g0est, lwd=2,col="darkred",lty=2)
  lgd <- legend("topright",legend=c(expression(Y^1),expression(Y^0)),lwd=2,lty=1:2,col=c("darkblue","darkred"),bg="transparent")
  points(rep((lgd$rect$left + lgd$text$x[1])/2,2), lgd$text$y, pch=15,cex=2.25,col=c("lightblue3","lightsalmon"))
  lgd <- legend("topright",legend=c(expression(Y^1),expression(Y^0)),lwd=2,lty=1:2,col=c("darkblue","darkred"),bg="transparent")
  }

  if (l2==T){ print(res) } else { res <- NULL }
  return(invisible(list(res=res, g1=g1, g0=g0, ifvals=ifvals)))

}
