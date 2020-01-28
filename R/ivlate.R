
#' @title Estimating complier average effect of binary treatment using binary instrument
#'
#' @description \code{ivlate} is used to estimate the mean outcome among compliers (i.e., those encouraged by the instrument) had all subjects received treatment versus control.
#'
#' @usage ivlate(y, a, z, x, nsplits=2, sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet",
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
#' ivlate.res <- ivlate(y,a,z,x)
#'
#' @references (Also see references for function \code{ate})
#' @references Angrist JD, Imbens GW, Rubin DB (1996). Identification of causal effects using instrumental variables. \emph{Journal of the American Statistical Association}.
#' @references Abadie A (2003). Semiparametric instrumental variable estimation of treatment response models. \emph{Journal of Econometrics}.
#' @references Kennedy EH, Balakrishnan S, G'Sell M (2017). Complier classification with sharp instrumental variables. \emph{Working Paper}.
#'
ivlate <- function(y,a,z,x, nsplits=2,
  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger","SL.rpart"),
  project01=T){

require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require("rpart")

expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
logx <- function(x){ 2*((x-1)/(x+1) + ((x-1)/(x+1))^3/3 + ((x-1)/(x+1))^5/5
          + ((x-1)/(x+1))^7/7 + ((x-1)/(x+1))^9/9 + ((x-1)/(x+1))^11/11)  }
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

# estimate treatment regression
la1fit <- SuperLearner(a[z==1 & train],
  as.data.frame(x[z==1 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)#,family=binomial)
la1hat[test] <- la1fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

if (!onesided){
la0fit <- SuperLearner(a[z==0 & train],
  as.data.frame(x[z==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)#,family=binomial)
la0hat[test] <- la0fit$SL.predict }
if (onesided){ la0hat[test] <- 0 }

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate outcome regression
yfam <- "gaussian"; if (length(unique(y))==2){ yfam <- "binomial" }
mu1fit <- SuperLearner(y[z==1 & train],
  as.data.frame(x[z==1 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)#, family=yfam)
mu1hat[test] <- mu1fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

mu0fit <- SuperLearner(y[z==0 & train],
  as.data.frame(x[z==0 & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)#, family=yfam)
mu0hat[test] <- mu0fit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
}

# project onto 0-1 and la1hat>la0hat space
if (project01){
newla <- project.01(la0hat,la1hat)
la0hat <- newla[,1]; la1hat <- newla[,2] }

ifvals.out <- z*(y-mu1hat)/pihat - (1-z)*(y-mu0hat)/(1-pihat) + mu1hat-mu0hat
ifvals.trt <- z*(a-la1hat)/pihat - (1-z)*(a-la0hat)/(1-pihat) + la1hat-la0hat
ifvals.gam2 <- 2*(la1hat-la0hat)*(z*(a-la1hat)/pihat - (1-z)*(a-la0hat)/(1-pihat)) + (la1hat-la0hat)^2

psihat <- mean(ifvals.out)/mean(ifvals.trt)
ifvals <- (ifvals.out - psihat*ifvals.trt)/mean(ifvals.trt)
muhat <- mean(ifvals.trt); xihat <- mean(ifvals.gam2)

muhat2 <- muhat*(muhat>0.01 & muhat<0.99) + 0.01*(muhat<0.01) + 0.99*(muhat>0.99)
q <- quantile(la1hat-la0hat,probs=1-muhat2)
xihat2 <- mean(ifvals.trt*(la1hat-la0hat>q))
sharp2 <- (xihat2-muhat^2)/(muhat-muhat^2)
ifvals.sharp2 <- (ifvals.trt*((la1hat-la0hat>q)+q)-xihat2 - q*((la1hat-la0hat>q)-muhat2) )/(muhat-muhat^2) + (2*muhat*xihat2 - xihat2 - muhat^2)*(ifvals.trt - muhat)/((muhat-muhat^2)^2)
if (sharp2<0.001){ sharp2 <- 0.001 }; if (sharp2>0.999){ sharp2 <- 0.999 }

est <- c(psihat, muhat, sharp2 )
se <- c(sd(ifvals), sd(ifvals.trt), sd(ifvals.sharp2))/sqrt(n)
ci.ll <- est-1.96*se; ci.ul <- est+1.96*se
ci.ll[3] <- expit(logit(sharp2) - 1.96*sd(ifvals.sharp2/(sharp2-sharp2^2))/sqrt(n))
ci.ul[3] <- expit(logit(sharp2) + 1.96*sd(ifvals.sharp2/(sharp2-sharp2^2))/sqrt(n))
pval <- round(2*(1-pnorm(abs(est/se))),3); pval[2:3] <- NA
param <- c("LATE","Strength", "Sharpness")
res <- data.frame(parameter=param, est,se,ci.ll,ci.ul,pval)
rownames(res) <- NULL

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount)
close(pb)

nuis <- data.frame(pi=pihat,la1=la1hat,la0=la0hat,mu1=mu1hat,mu0=mu0hat)

print(res)
return(invisible(list(res=res, nuis=nuis, ifvals=as.data.frame(ifvals) )))

}
