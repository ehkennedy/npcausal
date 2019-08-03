
#' @title Estimating average effect curve for continuous treatment
#'
#' @description \code{ctseff} is used to estimate the mean outcomes in a population had all subjects received given levels of a continuous (unconfounded) treatment.
#'
#' @usage ctseff(y, a, x, bw.seq, sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet",
#'   "SL.glm.interaction","SL.mean","SL.ranger"))
#'

#' @param y outcome of interest.
#' @param a continuous treatment.
#' @param x covariate matrix.
#' @param bw.seq sequence of bandwidth values.
#' @param sl.lib algorithm library for SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", and "ranger".
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs for population means.}
#' \item{bw.risk}{ estimated risk at sequence of bandwidth values.}
#'
#' @examples
#' n <- 500; x <- matrix(rnorm(n*5),nrow=n)
#' a <- runif(n); y <- a + rnorm(n,sd=.5)
#'
#' ce.res <- ctseff(y,a,x, bw.seq=seq(.2,2,length.out=100))
#' plot.ctseff(ce.res)
#'
#' # check that bandwidth choice is minimizer
#' plot(ce.res$bw.risk$bw,ce.res$bw.risk$risk)
#'
#' @references Kennedy EH, Ma Z, McHugh MD, Small DS (2017). Nonparametric methods for doubly robust estimation of continuous treatment effects. \emph{Journal of the Royal Statistical Society, Series B}. \href{https://arxiv.org/abs/1507.00747}{arxiv:1507.00747}
#'
ctseff <- function(y,a,x, bw.seq, n.pts=100, a.rng=c(min(a),max(a)),
  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger")){

require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require(KernSmooth)
kern <- function(t){ dnorm(t) }

n <- dim(x)[1]

# set up evaluation points & matrices for predictions
a.min <- a.rng[1]; a.max <- a.rng[2]; a.vals <- seq(a.min,a.max,length.out=n.pts)
xa.new <- rbind(cbind(x,a), cbind( x[rep(1:n,length(a.vals)),],a=rep(a.vals,rep(n,length(a.vals))) ))
x.new <- xa.new[,-dim(xa.new)[2]]
x <- data.frame(x); x.new <- data.frame(x.new); colnames(x) <- colnames(x.new)
xa.new <- data.frame(xa.new)

# estimate nuisance functions via super learner
# note: other methods could be used here instead
pimod <- SuperLearner(Y=a, X=data.frame(x), SL.library=sl.lib, newX=x.new); pimod.vals <- pimod$SL.predict
pi2mod <- SuperLearner(Y=(a-pimod.vals[1:n])^2,X=x, SL.library=sl.lib, newX=x.new)
pi2mod.vals <- pi2mod$SL.predict
mumod <- SuperLearner(Y=y, X=cbind(x,a), SL.library=sl.lib,newX=xa.new)
muhat.vals <- mumod$SL.predict

# construct estimated pi/varpi and mu/m values
a.std <- (xa.new$a-pimod.vals)/sqrt(pi2mod.vals)
pihat.vals <- approx(density(a.std)$x,density(a.std[1:n])$y,xout=a.std)$y
pihat <- pihat.vals[1:n]; pihat.mat <- matrix(pihat.vals[-(1:n)], nrow=n,ncol=length(a.vals))
varpihat <- predict(smooth.spline(a.vals, apply(pihat.mat,2,mean)), x=a)$y
varpihat.mat <- matrix( rep(apply(pihat.mat,2,mean),n), byrow=T, nrow=n)
muhat <- muhat.vals[1:n]; muhat.mat <- matrix(muhat.vals[-(1:n)], nrow=n,ncol=length(a.vals))
mhat <- predict(smooth.spline(a.vals, apply(muhat.mat,2,mean)), x=a)$y
mhat.mat <- matrix( rep(apply(muhat.mat,2,mean),n), byrow=T, nrow=n)

# form adjusted/pseudo outcome xi
pseudo.out <- (y-muhat)/(pihat/varpihat) + mhat

# leave-one-out cross-validation to select bandwidth
w.fn <- function(bw){ w.avals <- NULL; for (a.val in a.vals){
  a.std <- (a-a.val)/bw; kern.std <- kern(a.std)/bw
  w.avals <- c(w.avals, mean(a.std^2*kern.std)*(kern(0)/bw) /
    (mean(kern.std)*mean(a.std^2*kern.std)-mean(a.std*kern.std)^2))
}; return(w.avals/n) }
hatvals <- function(bw){ approx(a.vals,w.fn(bw),xout=a)$y }
cts.eff.fn <- function(out,bw){
  approx(locpoly(a,out,bandwidth=bw),xout=a)$y }
# note: choice of bandwidth range depends on specific problem,
# make sure to inspect plot of risk as function of bandwidth
risk.fn <- function(h){ hats <- hatvals(h); mean( ((pseudo.out - cts.eff.fn(pseudo.out,bw=h))/(1-hats))^2) }
risk.est <- sapply(bw.seq,risk.fn); h.opt <- bw.seq[which.min(risk.est)]
bw.risk <- data.frame(bw=bw.seq, risk=risk.est)
# alternative approach:
#h.opt <- optimize(function(h){ hats <- hatvals(h); mean( ((pseudo.out-cts.eff.fn(pseudo.out,bw=h))/(1-hats))^2) } ,
#  bw.seq, tol=0.01)$minimum

# estimate effect curve with optimal bandwidth
est <- approx(locpoly(a,pseudo.out,bandwidth=h.opt),xout=a.vals)$y

# estimate pointwise confidence band
# note: other methods could also be used
se <- NULL; for (a.val in a.vals){
  a.std <- (a-a.val)/h.opt; kern.std <- kern(a.std)/h.opt
  beta <- coef(lm(pseudo.out ~ a.std, weights=kern.std))
  Dh <- matrix( c(mean(kern.std), mean(kern.std*a.std),
    mean(kern.std*a.std), mean(kern.std*a.std^2)), nrow=2)
    kern.mat <- matrix(rep(kern((a.vals-a.val)/h.opt)/h.opt,n), byrow=T,nrow=n)
  g2 <- matrix( rep((a.vals-a.val)/h.opt, n), byrow=T, nrow=n)
  intfn1.mat <- kern.mat*(muhat.mat - mhat.mat)*varpihat.mat
  intfn2.mat <- g2*kern.mat*(muhat.mat - mhat.mat)*varpihat.mat
  int1 <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)])/2,n),
    byrow=T,nrow=n)*intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)],1,sum)
  int2 <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)])/2,n),
    byrow=T,nrow=n)*intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)],1,sum)
  sigma <- cov(t(solve(Dh) %*%
    rbind( kern.std*(pseudo.out-beta[1]-beta[2]*a.std) + int1,
    a.std*kern.std*(pseudo.out-beta[1]-beta[2]*a.std) + int2 )))
se <- c(se, sqrt(sigma[1,1])) }

ci.ll <- est-1.96*se/sqrt(n); ci.ul <- est+1.96*se/sqrt(n)
res <- data.frame(a.vals, est,se,ci.ll,ci.ul)

return(invisible(list(res=res,bw.risk=bw.risk)))

}

