\name{ph.moment}
\alias{ph.mean}
\alias{ph.var}
\alias{ph.moment}
\alias{ph.moment-method}
\alias{ph.moment,ANY,ph-method}
\alias{ph.moment,ANY,herlang-method}
\title{Moments for Phase-Type (PH) Distribution}
\description{
Moments for PH distribution.
}
\usage{
ph.mean(ph)
ph.var(ph)
\S4method{ph.moment}{ANY,ph}(k, ph, ...)
\S4method{ph.moment}{ANY,herlang}(k, ph, ...)
}
\arguments{
\item{ph}{an object of S4 class of PH (\code{\linkS4class{ph}}) or Hyper-Erlang (\code{\linkS4class{herlang}}).}
\item{k}{an integer of dgrees of moments.}
\item{\dots}{further arguments for methods.}
}
\details{
  The PH distribution with parameters \eqn{alpha}, \eqn{Q} and \eqn{xi}:
  k-th moment; \deqn{(-1)^k k! alpha (-Q)^{-k} 1}
}
\note{
\code{ph.moment} is a generic function for \code{\linkS4class{ph}} and \code{\linkS4class{herlang}}.
}
\value{
\code{ph.mean} and \code{ph.var} give mean and variance of PH.
\code{ph.moment} gives a vector of up to k moments.
}

\seealso{\code{\link{ph}}, \code{\link{cf1}}, \code{\link{herlang}}}

\examples{
## create a PH with specific parameters
(param1 <- ph(alpha=c(1,0,0), 
	          Q=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-1)), 
	          xi=c(2,2,0)))

## create a CF1 with specific parameters
(param2 <- cf1(alpha=c(1,0,0), rate=c(1.0,2.0,3.0)))

## create a hyper Erlang with specific parameters
(param3 <- herlang(shape=c(2,3), mixrate=c(0.3,0.7),
	               rate=c(1.0,10.0)))

## mean
ph.mean(param1)
ph.mean(param2)
ph.mean(param3)

## variance
ph.var(param1)
ph.var(param2)
ph.var(param3)

## up to 5 moments 
ph.moment(5, param1)
ph.moment(5, param2)
ph.moment(5, param3)

}
\keyword{distribution}

