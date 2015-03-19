\name{map}
\alias{map}
\alias{mmpp}
\alias{gmmpp}
\title{Markovian Arrival Process (MAP)}
\description{
Functions to generate an object of \code{\linkS4class{map}}.
}
\usage{
map(size, alpha, D0, D1, class = "CsparseMatrix")
mmpp(size, class = "CsparseMatrix")
gmmpp(size, alpha, D0, D1, class = "dgeMatrix")
}
\arguments{
\item{size}{an integer for the number of phases.}
\item{alpha}{a vector of probabilities for determing an initial phase.}
\item{D0}{an object of Matrix class for the initesmal generator without arrivals.}
\item{D1}{an object of Matrix class for the initesmal generator with arrivals.}
\item{class}{name of Matrix class for \code{D0} and \code{D1}.}
}
\details{
  MAP parameters are \eqn{alpha}, \eqn{D0} and \eqn{D1}. \eqn{alpha} is the probability vector to determine an initial phase at time 0. \eqn{D0} is an infinitesimal generator of underlyinc continuous-time Markov chain (CTMC) without arrival. \eqn{D1} is an infinitesimal generator of CTMC with arrival. The infinitesimal generator of underlying CTMC becomes \eqn{D0+D1}. In the stationary case, \eqn{alpha} is often given by a stationary vector satisfying \eqn{alpha (D0+D1) = alpha}.

  \code{mmpp} generates an object of a specific MAP called MMPP.
  MMPP (Markov modulated Poisson process) is an MAP whose \eqn{D1} is given by a diagonal matrix. Unlike to general MAPs, MMPP never changes the phase at which an arrival occurs.

  \code{gmmpp} generates an object of \code{\linkS4class{gmmpp}}, which is exactly same as MMPP. In the estimation algorithm, \code{\linkS4class{gmmpp}} class uses an approximate method.
}
\note{
  \code{map} and \code{gmmpp} require either \code{size} or (\code{alpha}, \code{D0}, \code{D1}).
}
\value{
\code{map} gives an object of general MAP.
\code{mmpp} gives an object of MMPP with default parameters.
\code{gmmpp} gives an object of MMPP which uses an approximate estimation algorithm.
}

\seealso{
\code{\link{erhmm}}, \code{\link{map.mmoment}},
\code{\link{map.jmoment}}, \code{\link{map.acf}}
}

\examples{
## create an MAP (full matrix) with 5 phases
map(5)

## create an MAP (full matrix) with 5 phases
map(size=5)

## create an MMPP with 5 states
mmpp(5)

## create an MMPP with 5 states for approximate
## estimation
gmmpp(5)

## create an MAP with specific parameters
(param <- map(alpha=c(1,0,0),
	D0=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-4)),
	D1=rbind(c(1,1,0),c(1,0,1),c(2,0,1))))

## marginal moments of MAP
map.mmoment(k=3, map=param)

## joint moments of MAP
map.jmoment(lag=1, map=param)

## k-lag correlation
map.acf(map=param)

}

