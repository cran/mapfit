\name{phfit.point}
\alias{phfit.point}
\title{PH fitting with point data}

\description{
  estimates PH parameters from point data.
}
\usage{
phfit.point(ph, x, weights, method = c("all", "increment"),
  lbound = 1, ubound = NULL, control = list(), verbose = list(), ...)
}
\arguments{
  \item{ph}{an object of S4 class for MAP. The estimation algorithm is selected depending on thie class.}
  \item{x}{a vector for point data.}
  \item{weights}{a vector of weights for points.}
  \item{method}{the name of estimation method for hyper Erlang (\code{\linkS4class{herlang}}).}
  \item{lbound}{a value for lower limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).}
  \item{ubound}{a value for upper limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).}
  \item{control}{a list of parameters for controlling the fitting process.}
  \item{verbose}{a list of parameters for displaying the fitting process.}
  \item{\dots}{further arguments for methods.}
}
\value{
  returns a list with components, which is an object of S3 class \code{phfit.result};
  \item{model}{an object for estimated PH class (\code{\linkS4class{ph}}, \code{\linkS4class{cf1}}, \code{\linkS4class{herlang}}).}
  \item{llf}{a value of the maximum log-likelihood.}
  \item{df}{a value of degrees of freedom of the model.}
  \item{aic}{a value of Akaike information criterion.}
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
\seealso{
  \code{\link{phfit.group}}, \code{\link{phfit.density}}, \code{\linkS4class{ph}}, \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}
}
\examples{
## make sample
wsample <- rweibull(n=100, shape=2, scale=1)

## PH fitting for general PH
(result1 <- phfit.point(ph=ph(2), x=wsample))

## PH fitting for CF1
(result2 <- phfit.point(ph=cf1(2), x=wsample))

## PH fitting for hyper Erlang
(result3 <- phfit.point(ph=herlang(3), x=wsample))

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
