\name{map.acf}
\alias{map.mmoment}
\alias{map.jmoment}
\alias{map.acf}
\title{Moments for Markovian arrival pcess (MAP)}
\description{
Moments for MAP.
}
\usage{
map.mmoment(k, map)
map.jmoment(lag, map)
map.acf(map)
}
\arguments{
\item{map}{an object of S4 class of MAP (\code{\linkS4class{map}}, \code{\linkS4class{gmmpp}}).}
\item{k}{an integer of dgrees of moments.}
\item{lag}{an integer of time lag for corrleation.}
}
\details{
  MAP parameters are \eqn{alpha}, \eqn{D0} and \eqn{D1};
  \deqn{P = (-D0)^(-1) D1},
  \deqn{s P = s}. Then the moments for MAP are
  marginal moment; \deqn{m_k = k! s (-D0)^{-k} 1}
  joint moment; \deqn{s_{ij}(lag) = i! j! s (-D0)^{-i} P^lag (-D0)^{-j} 1}
  k-lag correlation (autocorrelation); \deqn{rho(lag) = (s_{11}(lag) - m_1^2)/(m_2 - m_1^2)}
}
% \note{
% \code{map.mmoment} is a generic function for \code{\linkS4class{ph}} and \code{\linkS4class{herlang}}.
% }
\value{
\code{map.mmoment} gives a vector of up to k moments.
\code{map.jmoment} gives a matrix of \eqn{s_{ij}(lag), i=1,..,n, j=1,..,n} where n is the size of phases.
\code{map.acf} gives a vector of up to n-lag correlation, where n is the size of phases.
}

\seealso{\code{\link{map}}, \code{\link{gmmpp}}, \code{\link{erhmm}}}

\examples{
## create an MAP with specific parameters
(param1 <- map(alpha=c(1,0,0),
	D0=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-4)),
	D1=rbind(c(1,1,0),c(1,0,1),c(2,0,1))))

## create an ER-HMM with specific parameters
(param2 <- erhmm(shape=c(2,3), alpha=c(0.3,0.7),
                rate=c(1.0,10.0),
                P=rbind(c(0.3, 0.7), c(0.1, 0.9))))

## marginal moments of MAP
map.mmoment(k=3, map=param1)
map.mmoment(k=3, map=as(param2, "map"))

## joint moments of MAP
map.jmoment(lag=1, map=param1)
map.jmoment(lag=1, map=as(param2, "map"))

## k-lag correlation
map.acf(map=param1)
map.acf(map=as(param2, "map"))

}
\keyword{distribution}

