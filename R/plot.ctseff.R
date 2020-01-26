
#' @title Plot estimated average effect curve for continuous treatment
#'
#' @description \code{plot.ctseff} is used to plot results from \code{ctseff} fit.
#'
#' @usage plot.ctseff(ctseff.res)
#'
#' @param ctseff.res output from \code{ctseff} fit.
#'
#' @return A plot of estimated effect curve with pointwise confidence intervals.
#'
#' @examples
#' n <- 500
#' x <- matrix(rnorm(n * 5), nrow = n)
#' a <- runif(n)
#' y <- a + rnorm(n, sd = .5)
#' 
#' ce.res <- ctseff(y, a, x, bw.seq = seq(.2, 2, length.out = 100))
#' plot.ctseff(ce.res)
#' @references Kennedy EH, Ma Z, McHugh MD, Small DS (2017). Nonparametric methods for doubly robust estimation of continuous treatment effects. \emph{Journal of the Royal Statistical Society, Series B}. \href{https://arxiv.org/abs/1507.00747}{arxiv:1507.00747}
#'
plot.ctseff <- function(res, xlab = "Treatment level A=a", ylab = expression("E( Y"^"a" ~ ")")) {
  plot(res$res$a.vals, res$res$est,
    type = "l",
    ylim = c(min(res$res$ci.ll), max(res$res$ci.ul)), xlab = xlab, ylab = ylab
  )
  lines(res$res$a.vals, res$res$ci.ll, lty = 2)
  lines(res$res$a.vals, res$res$ci.ul, lty = 2)
}
