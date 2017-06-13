
#' @title Estimating average effect of treatment on the treated
#'
#' @description \code{att} is used to estimate the difference in mean outcome among treated subjects had a binary (unconfounded) treatment been withheld.
#'
#' @usage att(y, a, x, nsplits=2, sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet",
#'   "SL.glm.interaction","SL.mean","SL.ranger"))
#'
#' @param y outcome of interest.
#' @param a binary treatment.
#' @param x covariate matrix.
#' @param nsplits integer number of sample splits for nuisance estimation.
#' If nsplits=1, sample splitting is not used, and nuisance functions are estimated
#' on full sample (in which case validity of SEs/CIs requires empirical
#' process conditions). Otherwise must have nsplits>1.
#' @param sl.lib algorithm library if using SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", and "ranger".
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for treated means and contrast.}
#' \item{nuis}{ subject-specific estimates of nuisance functions (i.e., propensity score and outcome regression) }
#' \item{ifvals}{ vector of estimated influence function values.}
#'
#' @examples
#' n <- 1000; x <- matrix(rnorm(n*5),nrow=n)
#' a <- rbinom(n,1,.3); y <- rnorm(n)
#'
#' att.res <- att(y,a,x)
#'
#' @references (Also see references for function \code{ate})
#' @references Kennedy EH, Sjolander A, Small DS (2015). Semiparametric causal inference in matched cohort studies. \emph{Biometrika}.
#'
att <- function(y,a,x, nsplits=2,
  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean",
  "SL.ranger")){

require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require("rpart")

n <- dim(x)[1]
pb <- txtProgressBar(min=0, max=2*nsplits, style=3)

s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
pihat <- rep(NA,n); mu0hat <- rep(NA,n)

pbcount <- 0
Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
for (vfold in 1:nsplits){

train <- s!=vfold; test <- s==vfold
if (nsplits==1){ train <- test }

# estimate propensity score
pifit <- SuperLearner(a[train],as.data.frame(x[train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib, family=binomial)
pihat[test] <- pifit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate regression function
mu0fit <- SuperLearner(y[a==0 & train],
  as.data.frame(x[a==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
mu0hat[test] <- mu0fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
}

ey01hat <- mean((a/mean(a))*mu0hat + ((1-a)/mean(1-a))*(y-mu0hat)*pihat/(1-pihat))
psihat <- mean((a/mean(a))*(y-mu0hat) - ((1-a)/mean(1-a))*(y-mu0hat)*pihat/(1-pihat))

ifvals <- cbind( a*(y-mean(y))/mean(a),
  (a/mean(a))*(mu0hat - ey01hat) + ((1-a)/mean(1-a))*(y-mu0hat)*pihat/(1-pihat),
  (a/mean(a))*(y-mu0hat - psihat) - ((1-a)/mean(1-a))*(y-mu0hat)*pihat/(1-pihat) )

est <- c(mean(y[a==1]), ey01hat, psihat)
se <- apply(ifvals,2,sd)/sqrt(n)
ci.ll <- est-1.96*se; ci.ul <- est+1.96*se
pval <- round(2*(1-pnorm(abs(est/se))),3)
param <- c("E(Y|A=1)","E{Y(0)|A=1}", "E{Y-Y(0)|A=1}")
res <- data.frame(parameter=param, est,se,ci.ll,ci.ul,pval)
rownames(res) <- NULL

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount)
close(pb)

nuis <- data.frame(pi=pihat,mu0=mu0hat)

print(res)
return(invisible(list(res=res, nuis=nuis, ifvals=as.data.frame(ifvals) )))

}
