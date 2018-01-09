
#' @title Estimating bounds on treatment effects with instrumental variables
#'
#' @description \code{ivbds} is used to estimate bounds on various effects using instrumental variables.
#'
#' @usage ivbds(y, a, z, x, nsplits=2, sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet",
#'   "SL.glm.interaction", "SL.mean","SL.ranger","rpart"), project01=T)
#'
#' @param y outcome of interest.
#' @param a binary treatment.
#' @param z binary instrument.
#' @param x covariate matrix.
#' @param nsplits integer number of sample splits for nuisance estimation.
#' If nsplits=1, sample splitting is not used, and nuisance functions are estimated
#' on full sample (in which case validity of SEs/CIs requires empirical
#' process conditions). Otherwise must have nsplits>1.
#' @param sl.lib algorithm library for SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", "ranger", "rpart".
#' @param project01 should the estimated compliance score be projected to space respecting 0-1 bounds and monotonicity?
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for local average treatment effect E(Y(a=1)-Y(a=0)|A(z=1)>A(z=0)), as well as IV strength and sharpness.}
#' \item{nuis}{ subject-specific estimates of nuisance functions (i.e., IV propensity score and treatment/outcome regressions) }
#' \item{ifvals}{ matrix of estimated influence function values.}
#'
#' @examples
#' n <- 100; x <- matrix(rnorm(n*5),nrow=n)
#' z <- rbinom(n,1,0.5); a <- rbinom(n,1,0.6*z+0.2)
#' y <- rnorm(n)
#'
#' ivbds.res <- ivbds(y,a,z,x)
#'
#' @references (Also see references for function \code{ate})
#' @references Angrist JD, Imbens GW, Rubin DB (1996). Identification of causal effects using instrumental variables. \emph{Journal of the American Statistical Association}.
#' @references Abadie A (2003). Semiparametric instrumental variable estimation of treatment response models. \emph{Journal of Econometrics}.
#' @references Kennedy EH, Balakrishnan S, G'Sell M (2017). Complier classification with sharp instrumental variables. \emph{Working Paper}.
#'
ivbds <- function(y,a,z,x, nsplits=2,
  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger","SL.rpart"),
  project01=T){

require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require("rpart")

expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
project.01 <- function(x1,y1){
  x2 <- x1; y2 <- y1
  y2[y1>1 & 0<x1 & x1<1] <- 1
  x2[x1<0 & y1>1] <- 0; y2[x1<0 & y1>1] <- 1
  x2[x1<0 & 0<y1 & 11] <- 0
  x2[-x1>y1 & 10] <- 0; y2[-x1>y1 & y1<0] <- 0
  x2[x1>y1 & y1>-x1 & y1<2-x1] <- (x1+y1)[x1>y1 & y1>-x1 & y1<2-x1]/2
  y2[x1>y1 & y1>-x1 & y1<2-x1] <- (x1+y1)[x1>y1 & y1>-x1 & y1<2-x1]/2
  x2[y1>2-x1 & x1>1] <- 1; y2[y1>2-x1 & x1>1] <- 1
  return(cbind(x2,y2)) }

n <- dim(x)[1]
pb <- txtProgressBar(min=0, max=4*nsplits, style=3)

s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
pihat <- rep(NA,n); la1hat <- rep(NA,n); la0hat <- rep(NA,n)
vl1hat <- rep(NA,n); vl0hat <- rep(NA,n)
vu1hat <- rep(NA,n); vu0hat <- rep(NA,n)
mu1hat <- rep(NA,n); mu0hat <- rep(NA,n)
onesided <- sum(a[z==0])==0

pbcount <- 0
Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
for (vfold in 1:nsplits){

train <- s!=vfold; test <- s==vfold
if (nsplits==1){ train <- test }

# estimate iv propensity score
pifit <- SuperLearner(z[train],as.data.frame(x[train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib, family=binomial)
pihat[test] <- pifit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate treatment regressions
la1fit <- SuperLearner(a[z==1 & train],
  as.data.frame(x[z==1 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
la1hat[test] <- la1fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

if (!onesided){
la0fit <- SuperLearner(a[z==0 & train],
  as.data.frame(x[z==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
la0hat[test] <- la0fit$SL.predict }
if (onesided){ la0hat[test] <- 0 }

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate outcome regression
mu1fit <- SuperLearner(y[z==1 & train],
  as.data.frame(x[z==1 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
mu1hat[test] <- mu1fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

mu0fit <- SuperLearner(y[z==0 & train],
  as.data.frame(x[z==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
mu0hat[test] <- mu0fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate bound variable regressions

vu1 <- y*a + 1-a; vu0 <- y*(1-a)
vl1 <- y*a; vl0 <- y*(1-a)+a

vl1fit <- SuperLearner(vl1[z==1 & train],
  as.data.frame(x[z==1 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
vl1hat[test] <- vl1fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

vl0fit <- SuperLearner(vl0[z==0 & train],
  as.data.frame(x[z==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
vl0hat[test] <- vl0fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

vu1fit <- SuperLearner(vu1[z==1 & train],
  as.data.frame(x[z==1 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
vu1hat[test] <- vu1fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

vu0fit <- SuperLearner(vu0[z==0 & train],
  as.data.frame(x[z==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
vu0hat[test] <- vu0fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
}

# project onto 0-1 and la1hat>la0hat space
if (project01){
newla <- project.01(la0hat,la1hat)
la0hat <- newla[,1]; la1hat <- newla[,2] }

# compute if vals
ifvals.atel <- z*(vl1-vl1hat)/pihat - (1-z)*(vl0-vl0hat)/(1-pihat) + vl1hat-vl0hat
ifvals.ateu <- z*(vu1-vu1hat)/pihat - (1-z)*(vu0-vu0hat)/(1-pihat) + vu1hat-vu0hat
ifvals.trt <- z*(a-la1hat)/pihat - (1-z)*(a-la0hat)/(1-pihat) + la1hat-la0hat
ifvals.out <- z*(y-mu1hat)/pihat - (1-z)*(y-mu0hat)/(1-pihat) + mu1hat-mu0hat
muhat <- mean(ifvals.trt)
muhat2 <- muhat*(muhat>0.01 & muhat<0.99) + 0.01*(muhat<0.01) + 0.99*(muhat>0.99)
q <- quantile(la1hat-la0hat,probs=1-muhat2)
hq <- as.numeric( la1hat-la0hat > q)
ifvals.hql.num <- (z*(vl1-vl1hat)/pihat - (1-z)*(vl0-vl0hat)/(1-pihat) + vl1hat-vl0hat)*hq
ifvals.hqu.num <- (z*(vu1-vu1hat)/pihat - (1-z)*(vu0-vu0hat)/(1-pihat) + vu1hat-vu0hat)*hq
nulhat <- vl1hat-vl0hat; nuuhat <- vu1hat-vu0hat; gamhat <- la1hat-la0hat

latehat <- mean(ifvals.out)/mean(ifvals.trt)
bl.ate <- mean(ifvals.atel); bu.ate <- mean(ifvals.ateu)
bl.hq <- mean(ifvals.hql.num)/mean(ifvals.trt); bu.hq <- mean(ifvals.hqu.num)/mean(ifvals.trt)

if (bl.ate<= -1){ bl.ate <- -1 }; if (bl.ate> 1){ bl.ate <- 1 }
if (bu.ate<= -1){ bu.ate <- -1 }; if (bu.ate> 1){ bu.ate <- 1 }
if (bl.hq <= -1){ bl.hq <- -1 }; if (bl.hq > 1){ bl.hq <- 1 }
if (bu.hq <= -1){ bu.hq <- -1 }; if (bu.hq > 1){ bu.hq <- 1 }

ifvals.late <- (ifvals.out - latehat*ifvals.trt)/mean(ifvals.trt)
ifvals.hql <- (ifvals.hql.num - bl.hq*ifvals.trt)/mean(ifvals.trt)
ifvals.hqu <- (ifvals.hqu.num - bu.hq*ifvals.trt)/mean(ifvals.trt)

lb <- c(bl.ate,bl.hq); ub <- c(bu.ate,bu.hq)
se.l <- c(sd(ifvals.atel), sd(ifvals.hql))/sqrt(n)
se.u <- c(sd(ifvals.ateu), sd(ifvals.hqu))/sqrt(n)
cval <- rep(NA,2); for (j in 1:2){
crit <- function(cn){
  abs(pnorm(cn + (ub[j]-lb[j])/max(se.l[j],se.u[j])) - pnorm(-cn) - 0.95) }
(cval[j] <- optimize(crit,c(1,3))$min) }
ci.ll <- lb-cval*se.l; ci.ul <- ub+cval*se.u
lb <- c(lb,latehat); ub <- c(ub,latehat)
ci.ll <- c(ci.ll,latehat-1.96*sd(ifvals.late)/sqrt(n))
ci.ul <- c(ci.ul,latehat+1.96*sd(ifvals.late)/sqrt(n))
param <- c("ATE", "beta(h_q)", "LATE")
res <- data.frame(parameter=param, lb,ub,ci.ll,ci.ul)
rownames(res) <- NULL

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount)
close(pb)

nuis <- data.frame(pi=pihat,la1=la1hat,la0=la0hat,hq,gamhat=la1hat-la0hat)

print(res)
return(invisible(list(res=res, nuis=nuis,
  ifvals=data.frame(ifvals.atel,ifvals.ateu,ifvals.hql,ifvals.hqu) )))

}
