
#' @title Doubly robust estimation of average treatment effect
#'
#' @description \code{ate} is used to estimate the mean outcome in a population had
#'  all subjects received given levels of a discrete (unconfounded) treatment, using
#'  doubly robust methods with ensembled nuisance estimation.
#'
#' @usage ate(y, a, x, nsplits=2, sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glmnet",
#'   "SL.glm.interaction", "SL.mean","SL.ranger","rpart"))
#'
#' @param y outcome of interest.
#' @param a discrete treatment.
#' @param x covariate matrix.
#' @param nsplits integer number of sample splits for nuisance estimation.
#' If nsplits=1, sample splitting is not used, and nuisance functions are estimated
#' on full sample (in which case validity of SEs/CIs requires empirical
#' process conditions). Otherwise must have nsplits>1.
#' @param sl.lib algorithm library for SuperLearner.
#' Default library includes "earth", "gam", "glm", "glmnet", "glm.interaction",
#' "mean", "ranger", "rpart.
#'
#' @return A list containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for population means and relevant contrasts.}
#' \item{nuis}{ subject-specific estimates of nuisance functions (i.e., propensity score and outcome regression) }
#' \item{ifvals}{ matrix of estimated influence function values.}
#'
#' @examples
#' n <- 100; x <- matrix(rnorm(n*5),nrow=n)
#' a <- sample(3,n,replace=TRUE); y <- rnorm(n,mean=x[,1])
#'
#' ate.res <- ate(y,a,x, sl.lib=c("SL.mean","SL.gam"))
#'
#' @references Robins JM, Rotnitzky A (1995). Semiparametric efficiency in multivariate regression models with missing data. \emph{Journal of the American Statistical Association}.
#' @references Hahn J (1998). On the role of the propensity score in efficient semiparametric estimation of average treatment effects. \emph{Econometrica}.
#' @references van der Laan MJ, Robins JM (2003). \emph{Unified Methods for Censored Longitudinal Data and Causality} (Springer).
#' @references Tsiatis AA (2006). \emph{Semiparametric Theory and Missing Data} (Springer).
#' @references Robins JM, Li L, Tchetgen Tchetgen ET, van der Vaart A (2008). Higher order influence functions and minimax estimation of nonlinear functionals. \emph{Probability and Statistics: Essays in Honor of David A. Freedman}.
#' @references Zheng W, van der Laan (2010). Asymptotic theory for cross-validated targeted maximum likelihood estimation \emph{UC Berkeley Division of Biostatistics Working Paper Series}.
#' @references Chernozhukov V, Chetverikov V, Demirer M, et al (2016). Double machine learning for treatment and causal parameters.
#'
ate <- function(y,a,x, nsplits=2,
  sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger","SL.rpart"),
  ps=NULL){

require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require("rpart")

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
  if (i==1){setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1 }
for (vfold in 1:nsplits){

train <- s!=vfold; test <- s==vfold
if (nsplits==1){ train <- test }

# estimate propensity score
if (i != n.avals & is.null(ps)){
pifit <- SuperLearner(as.numeric(a==avals[i])[train],as.data.frame(x[train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib, family=binomial)
pihat[test,i] <-pifit$SL.predict }


setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1

# estimate regression function
mufit <- SuperLearner(y[a==avals[i] & train],
  as.data.frame(x[a==avals[i] & train,]),
  newX=as.data.frame(x[test,]), SL.library=sl.lib)
muhat[test,i] <- mufit$SL.predict


setTxtProgressBar(pb,pbcount); pbcount <- pbcount+1
}
if (i == n.avals){ pihat[,i] <- 1 - apply(pihat,1,sum, na.rm=T) }
}

amat <- matrix(rep(a,n.avals),nrow=n,byrow=F)
alevel <- matrix(rep(as.numeric(avals), rep(n,n.avals)),nrow=n,byrow=F)
ymat <- matrix(rep(y,n.avals),nrow=n,byrow=F)

if (!is.null(ps)){ pihat <- amat*ps + (1-amat)*(1-ps) }
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


setTxtProgressBar(pb,pbcount)
close(pb)

nuis <- as.data.frame(cbind(pihat,muhat))
colnames(nuis) <- paste(rep(c("pi","mu"), rep(n.avals,2)), colnames(nuis), sep="_")

print(res)
return(invisible(list(res=res, nuis=nuis, ifvals=as.data.frame(ifvals) )))

}
