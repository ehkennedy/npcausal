% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.cdensity.R
\name{cv.cdensity}
\alias{cv.cdensity}
\title{Cross-validation for doubly robust estimation of counterfactual densities}
\usage{
cv.cdensity(y, a, x, kmax=5,
 gridlen=20, nsplits=2, progress_updates = TRUE)
}
\arguments{
\item{y}{outcome of interest.}

\item{a}{binary treatment (more than 2 levels are allowed, but only densities under
A=1 and A=0 will be estimated).}

\item{x}{covariate matrix.}

\item{kmax}{Integer indicating maximum dimension of (cosine) basis
expansion that should be used in series estimator.}

\item{gridlen}{Integer number indicating length of grid for which the
plug-in estimator of the marginal density is computed.}

\item{nsplits}{Integer number of sample splits for nuisance estimation. If
\code{nsplits = 1}, sample splitting is not used, and nuisance functions are
estimated n full sample (in which case validity of standard errors and
confidence intervals requires empirical process conditions). Otherwise must
have \code{nsplits > 1}.}

\item{progress_updates}{A \code{logical} value indicating whether to print a
progress statement as various stages of computation reach completion.
The default is \code{TRUE}, printing a progress bar to inform the user.}
}
\value{
A plot of the pseudo L2 risk of candidate estimators for counterfactual
 densities, at each model dimension from 1 to kmax
}
\description{
\code{cv.cdensity} estimates counterfactual densities using
 linear cosine basis expansions at a sequence of dimensions, and then estimates
 the L2 pseudo-risk of each, which can be used for purposes of model selection.
 Nuisance functions are estimated with random forests.
}
\examples{
n <- 100; x <- matrix(rnorm(n*5),nrow=n)
a <- sample(3,n,replace=TRUE)-2; y <- rnorm(n)

cv.cdensity(y,a,x)

}
\references{
Kennedy EH, Wasserman LA, Balakrishnan S. Semiparametric counterfactual
density estimation. \href{https://arxiv.org/abs/TBA}{arxiv:TBA}
}
