\name{ph}
\alias{ph}
\alias{dph}
\alias{pph}
\alias{rph}
\title{Phase-Type (PH) Distribution}
\description{
Density function, distribution function and
random generation for the PH distribution, and
a function to generate an object of \code{\linkS4class{ph}}.
}
\usage{
ph(size, alpha, Q, xi, class = "CsparseMatrix")
dph(x, ph = ph(1), log = FALSE)
pph(q, ph = ph(1), lower.tail = TRUE, log.p = FALSE)
rph(n, ph = ph(1))
}
\arguments{
\item{size}{a value for the number of phases.}
\item{alpha}{a vector for the initial probabilities of PH distribution.}
\item{Q}{an object of Matrix class for the initesmal generator of PH distribution.}
\item{xi}{a vector for the exit rates of PH distribution.}
\item{class}{name of Matrix class for \code{Q}.}
\item{x, q}{vectors of quantiles.}
\item{p}{a vector of probabilities.}
\item{n}{number of observations.}
\item{ph}{an object of S4 class of PH (\code{\linkS4class{ph}}).}
\item{log}{logical; if \code{TRUE}, the log density is returned.}
\item{lower.tail}{logical; if \code{TRUE}, probabilities are P[X <= x], otherwise, P[X > x].}
\item{log.p}{logical; if \code{TRUE}, the log probability is returned.}
}
\details{
  The PH distribution with parameters \eqn{alpha}, \eqn{Q} and \eqn{xi}:
  Cumulative probability function; \deqn{F(q) = 1 - alpha \exp( Q q ) 1}
  Probability density function; \deqn{f(x) = alpha \exp( Q x ) xi}
}
\note{
  \code{ph} requires either \code{size} or (\code{alpha}, \code{Q}, \code{xi}).
  \code{rph} for \code{\linkS4class{ph}} is too slow. It is recommended to use
  \code{rph} for \code{\linkS4class{cf1}}.
}
\value{
\code{ph} gives an object of general PH distribution.
\code{dph} gives the density function, \code{pph}
gives the distribution function,
and \code{rph} generates random samples.
}

\seealso{\code{\link{cf1}}, \code{\link{herlang}}}

\examples{
## create a PH (full matrix) with 5 phases
(param1 <- ph(5))

## create a PH (full matrix) with 5 phases
(param1 <- ph(size=5))

## create a PH with specific parameters
(param2 <- ph(alpha=c(1,0,0), 
	          Q=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-1)),
	          xi=c(2,2,0))) 

## p.d.f. for 0, 0.1, ..., 1
(dph(x=seq(0, 1, 0.1), ph=param2))

## c.d.f. for 0, 0.1, ..., 1
(pph(q=seq(0, 1, 0.1), ph=param2))

## generate 10 samples
(rph(n=10, ph=param2))

}
\keyword{distribution}

