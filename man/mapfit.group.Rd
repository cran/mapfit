\name{mapfit.group}
\alias{mapfit.group}
\title{MAP fitting with grouped data}

\description{
  estimates MAP parameters from grouped data.
}
\usage{
mapfit.group (map, counts, breaks, intervals, instant,
  stationary = TRUE, control = list(), verbose = list(), ...)
}
\arguments{
  \item{map}{S4 class for MAP. The estimation algorithm is selected depending on thie class.}
  \item{counts}{a vector for the number of arrivals in time interval.}
  \item{breaks}{a vector for time sequence to determine time interval. This is equivalent to \code{c(0,cumsum(intervals))}. If this is missing, it is assigned to \code{0:length(counts)}.}
  \item{intervals}{a vector for a sequence of time length for intervals. This is equivalent to \code{diff(breaks)}). If this is missing, it is assigned to \code{rep(1,length(counts))}.}
  \item{instant}{a vector of integer to indicate whether an arrival occurs at the last time of interval. If instant is 1, an arrival occurs at the last time of interval. If instant is 0, no arrival occurs at the last time of interval. By using instant, time point data can be expressed by grouped data class. If instant is missing, it is given by \code{rep(0,length(counts))}, i.e., there are no arrivals at the end of interval.}
  \item{stationary}{a logical value that determine whether initial probability is given by a stationary vector of underlying Markov process or not.}
  \item{control}{a list of parameters for controlling the fitting process.}
  \item{verbose}{a list of parameters for displaying the fitting process.}
  \item{\dots}{further arguments for methods.}
}
\value{
  returns a list with components, which is an object of S3 class \code{mapfit.result};
  \item{model}{an object for estimated MAP class (\code{\linkS4class{map}}, \code{\linkS4class{erhmm}}).}
  \item{llf}{a value of the maximum log-likelihood.}
  \item{df}{a value of degrees of freedom of the model.}
  \item{aic}{a value of Akaike information criterion.}
  \item{iter}{the number of iterations.}
  \item{convergence}{a logical value for the convergence of estimation algorithm.}
  \item{ctime}{computation time (user time).}
  \item{stationary}{a logical value for the argument \code{stationary}.}
  \item{data}{an object for MAP data class}
  \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
  \item{rerror}{a value of relative error for llf at the last step of algorithm.}
  \item{control}{a list of the argument of \code{control}.}
  \item{verbose}{a list of the argument of \code{verbose}.}
  \item{call}{the matched call.}
}
\seealso{
  \code{\link{mapfit.point}}, \code{\linkS4class{map}} and \code{\linkS4class{gmmpp}}
}
\examples{
## load trace data
data(BCpAug89)
BCpAug89s <- head(BCpAug89, 50)

## make grouped data
BCpAug89.group <- hist(cumsum(BCpAug89s),
                       breaks=seq(0, 0.15, 0.005),
                       plot=FALSE)

## MAP fitting for general MAP
(result1 <- mapfit.group(map=map(2),
                         counts=BCpAug89.group$counts,
                         breaks=BCpAug89.group$breaks))

## MAP fitting for MMPP
(result2 <- mapfit.group(map=mmpp(2),
                         counts=BCpAug89.group$counts,
                         breaks=BCpAug89.group$breaks))

## MAP fitting with approximate MMPP
(result3 <- mapfit.group(map=gmmpp(2),
                         counts=BCpAug89.group$counts,
                         breaks=BCpAug89.group$breaks))

## marginal moments for estimated MAP
map.mmoment(k=3, map=result1$model)
map.mmoment(k=3, map=result2$model)
map.mmoment(k=3, map=result3$model)

## joint moments for estimated MAP
map.jmoment(lag=1, map=result1$model)
map.jmoment(lag=1, map=result2$model)
map.jmoment(lag=1, map=result3$model)

## lag-k correlation
map.acf(map=result1$model)
map.acf(map=result2$model)
map.acf(map=result3$model)

}
