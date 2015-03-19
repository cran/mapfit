\name{phfit.density}
\alias{phfit.density}
\title{PH fitting with density function}

\description{
  estimates PH parameters from density function.
}
\usage{
phfit.density(ph, f, deformula = zero.to.inf,
  weight.zero = .Machine$double.eps,
  weight.reltol = sqrt(.Machine$double.eps),
  method = c("all", "increment"), lbound = 1, ubound = NULL,
  control = list(), verbose = list(), ...)
}
\arguments{
  \item{ph}{an object of S4 class for MAP. The estimation algorithm is selected depending on thie class.}
  \item{f}{a faunction object for a density function.}
  \item{deformula}{an object for formulas of numerical integration. It is not necessary to change it when the density function is defined on the positive domain [0,infinity).}
  \item{weight.zero}{a absolute value which is regarded as zero in numerical integration.}
  \item{weight.reltol}{a value for precision of numerical integration.}
  \item{method}{the name of estimation method for hyper Erlang (\code{\linkS4class{herlang}}).}
  \item{lbound}{a value for lower limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).}
  \item{ubound}{a value for upper limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).}
  \item{control}{a list of parameters for controlling the fitting process.}
  \item{verbose}{a list of parameters for displaying the fitting process.}
  \item{\dots}{further arguments for methods, which are also used to send the arguments to density function.}
}
\value{
  returns a list with components, which is an object of S3 class \code{phfit.result};
  \item{model}{an object for estimated PH class (\code{\linkS4class{ph}}, \code{\linkS4class{cf1}}, \code{\linkS4class{herlang}}).}
  \item{llf}{a value of the maximum log-likelihood (a netative value of the cross entropy).}
  \item{df}{a value of degrees of freedom of the model.}
  \item{aic}{a value of Akaike information criterion (this is not meaningless in this case).}
  \item{iter}{the number of iterations.}
  \item{convergence}{a logical value for the convergence of estimation algorithm.}
  \item{ctime}{computation time (user time).}
  \item{data}{an object for MAP data class}
  \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
  \item{rerror}{a value of relative error for llf at the last step of algorithm.}
  \item{control}{a list of the argument of \code{control}.}
  \item{verbose}{a list of the argument of \code{verbose}.}
  \item{call}{the matched call.}
}
\note{
  Any of density function can be applied to the argument \code{f}, where \code{f} should be defined \code{f <- function(x, ...)}. The first argument of \code{f} should be an integral parameter. The other parameters are set in the argument \code{...} of \code{phfit.density}. The truncated density function can also be used directly.
}
\seealso{
  \code{\link{phfit.point}}, \code{\link{phfit.group}}, \code{\linkS4class{ph}}, \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}
}
\examples{
####################
##### truncated density
####################

## PH fitting for general PH
(result1 <- phfit.density(ph=ph(2), f=dnorm,
                          mean=3, sd=1))

## PH fitting for CF1
(result2 <- phfit.density(ph=cf1(2), f=dnorm,
                          mean=3, sd=1))

## PH fitting for hyper Erlang
(result3 <- phfit.density(ph=herlang(3), f=dnorm,
                          mean=3, sd=1))

## mean
ph.mean(result1$model)
ph.mean(result2$model)
ph.mean(result3$model)

## variance
ph.var(result1$model)
ph.var(result2$model)
ph.var(result3$model)

## up to 5 moments 
ph.moment(5, result1$model)
ph.moment(5, result2$model)
ph.moment(5, result3$model)

}
