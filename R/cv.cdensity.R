#' @title Cross-validation for doubly robust estimation of counterfactual densities
#'
#' @description \code{cv.cdensity} estimates counterfactual densities using
#'  linear cosine basis expansions at a sequence of dimensions, and then estimates
#'  the L2 pseudo-risk of each, which can be used for purposes of model selection.
#'  Nuisance functions are estimated with random forests.
#'
#' @usage cv.cdensity(y, a, x, kmax=5,
#'  gridlen=20, nsplits=2, progress_updates = TRUE)
#'
#' @param y outcome of interest.
#' @param a binary treatment (more than 2 levels are allowed, but only densities under
#'  A=1 and A=0 will be estimated).
#' @param x covariate matrix.
#' @param kmax Integer indicating maximum dimension of (cosine) basis
#'  expansion that should be used in series estimator.
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
#'
#' @return A plot of the pseudo L2 risk of candidate estimators for counterfactual
#'  densities, at each model dimension from 1 to kmax
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
#' cv.cdensity(y,a,x)
#'
#' @references Kennedy EH, Wasserman LA, Balakrishnan S. Semiparametric counterfactual
#' density estimation. \href{https://arxiv.org/abs/TBA}{arxiv:TBA}
#
cv.cdensity <- function(y, a, x,
                kmax=5, gridlen=20, nsplits=2,
                progress_updates = TRUE) {

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
  dat.all <- data.frame(y,a1,a0,x)

  # set grid for plugin estimator
  ygrid <- quantile(y[a1==1 | a0==1],probs=seq(0,1,length.out=gridlen))

  # select half for learning g / model selection
  train <- sample( rep(0:1,ceiling(n/2))[1:n] )

  for (set in 1:0){

  if (set==1 & progress_updates==T){ print("estimating candidate densities"); flush.console() }
  if (set==0 & progress_updates==T){ print("starting model selection"); flush.console() }

  # setup test data
  n <- sum(train==set); x <- dat.all[train==set,-(1:3)]
  y <- dat.all$y[train==set]; a1 <- dat.all$a1[train==set]; a0 <- dat.all$a0[train==set]

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

  if (set==1){
  betag1 <- apply((a1/pi1hat)*(bmat-mu1mat) + mu1mat,2,mean)
  betag0 <- apply((a0/pi0hat)*(bmat-mu0mat) + mu0mat,2,mean)
  g1fun <- function(k,t){ pos.part(1 + b(t,1:k) %*% betag1[1:k]) }
  g0fun <- function(k,t){ pos.part(1 + b(t,1:k) %*% betag0[1:k]) }
  g1kconst <- function(k){ integrate(function(t){ g1fun(k,t) },lower=0,upper=1)$val }
  g0kconst <- function(k){ integrate(function(t){ g0fun(k,t) },lower=0,upper=1)$val }
  g1 <- function(k,t){ g1fun(k,t) / g1kconst(k) }
  g0 <- function(k,t){ g0fun(k,t) / g0kconst(k) }
  }

  }

  ### estimate and plot pseudo L2 risk

  par(mfrow=c(1,2))

  for (trt in 1:0){

  if (trt==1){
    g <- g1; phat <- p1hat; mumat <- mu1mat; drmean <- apply((a1/pi1hat)*(bmat-mu1mat) + mu1mat,2,mean)
    aval <- a1; pihat <- pi1hat; title <- expression(Y^1," density") }
  if (trt==0){
    g <- g0; phat <- p0hat; mumat <- mu0mat; drmean <- apply((a0/pi0hat)*(bmat-mu0mat) + mu0mat,2,mean)
    aval <- a0; pihat <- pi0hat; title <- expression(Y^0," density") }

  # plugin
  delta.plugin <- rep(NA,kmax)
  for (j in 1:kmax){
    delta.plugin[j] <- integrate(function(t){ (g(j,t)-phat(t))^2 },lower=0,upper=1)$val }

  meth <- "shifted"

  # one-step estimator of L2
  if (meth=="L2"){
  delta.res <- data.frame(matrix(nrow=kmax,ncol=3))
  colnames(delta.res) <- c("est","ci.ll","ci.ul")
  for (j in 1:kmax){
    pyhat <- phat(y); gy <- g(j,y)
    gyregahat <- 1 + mumat %*% c(drmean[1:j],rep(0,kmax-j))
    delta.ifmean <- integrate(function(t){ (phat(t)-g(j,t))*phat(t) },lower=0,upper=1)$val
    delta.if <- 2 * ( (aval/pihat)*((pyhat-gy) - (pyregahat-gyregahat)) +
               (pyregahat-gyregahat) - delta.ifmean )
    delta.res$est[j] <- delta.plugin[j] + mean(delta.if)
    delta.res$ci.ll[j] <- delta.res$est[j] - 1.96*sqrt(var(delta.if)/n)
    delta.res$ci.ul[j] <- delta.res$est[j] + 1.96*sqrt(var(delta.if)/n)
  } }

  # shifted/pseudo L2
  if (meth=="shifted"){
  delta.res <- data.frame(matrix(nrow=kmax,ncol=3))
  colnames(delta.res) <- c("est","ci.ll","ci.ul")
  for (j in 1:kmax){
    gy <- g(j,y); gyregahat <- 1 + mumat %*% c(drmean[1:j],rep(0,kmax-j))
    gnorm <- integrate(function(t){ g(j,t)^2 },lower=0,upper=1)$val
    delta.if <- 2*( (aval/pihat)*(gy - gyregahat) + gyregahat )
    delta.res$est[j] <- gnorm - mean(delta.if)
    delta.res$ci.ll[j] <- delta.res$est[j] - 1.96*sqrt(var(delta.if)/n)
    delta.res$ci.ul[j] <- delta.res$est[j] + 1.96*sqrt(var(delta.if)/n)
  } }

  kseq <- 1:kmax
  plot(kseq,delta.res$est, ylim=range(delta.res),pch="",
     ylab="Pseudo L2 risk", xlab="Basis dimension",
     main=title)
  segments(kseq,delta.res$ci.ll,kseq,delta.res$ci.ul, col="lightgray")
  lines(kseq,delta.res$est); points(kseq,delta.res$est)

  }

  ### results


}
