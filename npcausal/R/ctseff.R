
#' @title Estimating average effect curve for continuous treatment
#'
#' @description \code{ctseff} is used to estimate the mean outcome in a population had all subjects received given levels of a continuous (unconfounded) treatment.
#'
#' @usage ctseff(x, a, y, nsplits=2,
#' sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet", "SL.glm.interaction","SL.mean","SL.ranger"))
#'
#' @param x covariate matrix.
#' @param a continuous treatment.
#' @param y outcome of interest.
#' @param nsplits integer number of sample splits for nuisance estimation.
#' If nsplits=1, sample splitting is not used, and nuisance functions are estimated
#' on full sample (in which case validity of SEs/CIs requires empirical
#' process conditions). Otherwise must have nsplits>1.
#' @param sl.lib algorithm library if using SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", and "ranger".
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for population means and relevant contrasts.}
#' \item{ifvals}{ matrix of estimated influence function values.}
#'
#' @examples
#' n <- 1000; x <- matrix(rnorm(n*4),nrow=n)
#' a <- runif(n); y <- rnorm(n)
#'
#' ctseff.res <- ctseff(x,a,y)
#'
#' @references Kennedy EH, Ma Z, McHugh MD, Small DS (2017). Nonparametric methods for doubly robust estimation of continuous treatment effects. \emph{Journal of the Royal Statistical Society, Series B}. \href{https://arxiv.org/abs/1507.00747}{arxiv:1507.00747}
#'
ctseff <- function(x, a, y, nsplits=2,
  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean",
  "SL.ranger")){

require("SuperLearner")
require("earth")
require("gam")
SL.ranger <- function (Y, X, newX, family, ...) { require("ranger")
  fit.rf <- ranger::ranger(Y ~ ., data=X); pred <- predict(fit.rf,data=newX)$predictions
  fit <- list(object = fit.rf); out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.ranger"); return(out) }

n <- dim(x)[1]
avals <- names(table(a))
n.avals <- length(avals)
pb <- txtProgressBar(min=0, max=2*nsplits*n.avals, style=3)

s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])

muhat <- as.data.frame(matrix(NA,nrow=n,ncol=n.avals))
colnames(muhat) <- paste("a",avals,sep="")
pihat <- muhat

pbcount <- 0
for (i in 1:n.avals){
  if (i==1){ Sys.sleep(0.1); setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
for (vfold in 1:nsplits){

train <- s!=vfold; test <- s==vfold
if (nsplits==1){ train <- test }

# estimate propensity score
if (i != n.avals){
pifit <- SuperLearner(as.numeric(a==avals[i])[train],as.data.frame(x[train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib, family=binomial)
pihat[test,i] <-pifit$SL.predict }

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate regression function
mufit <- SuperLearner(y[a==avals[i] & train],
  as.data.frame(x[a==avals[i] & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
muhat[test,i] <- mufit$SL.predict

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
}
if (i == n.avals){ pihat[,i] <- 1 - apply(pihat,1,sum, na.rm=T) }
}

amat <- matrix(rep(a,n.avals),nrow=n,byrow=F)
alevel <- matrix(rep(1:n.avals, rep(n,n.avals)),nrow=n,byrow=F)
ymat <- matrix(rep(y,n.avals),nrow=n,byrow=F)

ifvals <- as.matrix( (amat==alevel)*(ymat-muhat)/pihat + muhat )

est <- apply(ifvals,2,mean)
se <- apply(ifvals,2,sd)/sqrt(n)
ci.ll <- est-1.96*se; ci.ul <- est+1.96*se
pval <- round(2*(1-pnorm(abs(est/se))),3)
paste("E{Y(",avals,")}")
res1 <- data.frame(parameter=paste("E{Y(",avals,")}",sep=""), est,se,ci.ll,ci.ul,pval)

signdist <- function(x){ c(as.dist(outer(x,x,'-'))) }
ifvals2 <- t(apply(ifvals,1,signdist))
if (n.avals==2){ ifvals2 <- t(ifvals2)  }

tmp <- expand.grid(1:n.avals,1:n.avals)
tmp2 <- tmp[tmp[,1]>tmp[,2],]
contrasts <- apply(cbind(avals[tmp2[,1]],avals[tmp2[,2]]),1,paste,collapse=")-Y(")
contrasts <- paste("E{Y(",contrasts,")}",sep="")

est2 <- apply(ifvals2,2,mean)
se2 <- apply(ifvals2,2,sd)/sqrt(n)
ci.ll2 <- est2-1.96*se2; ci.ul2 <- est2+1.96*se2
pval2 <- round(2*(1-pnorm(abs(est2/se2))),3)
res2 <- data.frame(parameter=contrasts,est=est2,se=se2,ci.ll=ci.ll2,ci.ul=ci.ul2,pval=pval2)

res <- rbind(res1,res2); rownames(res) <- NULL

Sys.sleep(0.1)
setTxtProgressBar(pb,pbcount)
close(pb)

print(res)
return(invisible(list(res=res, ifvals=as.data.frame(ifvals))))

}
