\name{herlang}
\alias{herlang}
\alias{dherlang}
\alias{pherlang}
\alias{rherlang}
\title{Hyper-Erlang Distribution}
\description{
Density function, distribution function and
random generation for the hyper-Erlang distribution, and
a function to generate an object of \code{\linkS4class{herlang}}.
}
\usage{
herlang(shape, mixrate = rep(1/length(shape), length(shape)),
	     rate = rep(1, length(shape)))
dherlang(x, herlang = herlang(shape = c(1)), log = FALSE)
pherlang(q, herlang = herlang(shape = c(1)), lower.tail = TRUE, log.p = FALSE)
rherlang(n, herlang = herlang(shape = c(1)))
}
\arguments{
\item{shape}{a integer vector of shape parameters of Erlang components.}
\item{mixrate}{a vector for the initial probabilities of hyper-Erlang distribution.}
\item{rate}{a vector of rate parameters of Erlang components.}
\item{x, q}{vectors of quantiles.}
\item{p}{a vector of probabilities.}
\item{n}{number of observations.}
\item{herlang}{an object of S4 class of hyper Erlang (\code{\linkS4class{herlang}}).}
\item{log}{logical; if \code{TRUE}, the log density is returned.}
\item{lower.tail}{logical; if \code{TRUE}, probabilities are P[X <= x], otherwise, P[X > x].}
\item{log.p}{logical; if \code{TRUE}, the log probability is returned.}
}
\details{
  The hyper-Erlang distribution with parameters \eqn{mixrate}, \eqn{shape} and \eqn{rate}:
  Cumulative probability function; \deqn{F(q) = \sum_i \int_0^q maxirate_i rate_i^shape_i x^(shape_i-1) e^(-rate_i x) / (shape_i - 1)! dx}
  Probability density function; \deqn{f(x) = \sum_i maxirate_i rate_i^shape_i x^(shape_i-1) e^(-rate_i x) / (shape_i - 1)!}
}
\note{
  \code{herlang} requires shape parameters. 
}
\value{
\code{herlang} gives an object of hyper-Erlang distribution.
\code{dherlang} gives the density function, \code{pherlang} gives the distribution function,
and \code{rherlang} generates random samples.
}

\seealso{\code{\link{ph}}, \code{\link{herlang}}}

\examples{
## create a hyper Erlang consisting of two Erlang
## with shape parameters 2 and 3.
(param1 <- herlang(c(2,3)))

## create a hyper Erlang consisting of two Erlang
## with shape parameters 2 and 3.
(param1 <- herlang(shape=c(2,3)))

## create a hyper Erlang with specific parameters
(param2 <- herlang(shape=c(2,3), mixrate=c(0.3,0.7),
	               rate=c(1.0,10.0)))

## convert to a general PH
as(param2, "ph")

## p.d.f. for 0, 0.1, ..., 1
(dherlang(x=seq(0, 1, 0.1), herlang=param2))

## c.d.f. for 0, 0.1, ..., 1
(pherlang(q=seq(0, 1, 0.1), herlang=param2))

## generate 10 samples
(rherlang(n=10, herlang=param2))

}
\keyword{distribution}

