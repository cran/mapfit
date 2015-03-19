\name{cf1}
\alias{cf1}
\title{Canonical Form 1 for Phase-Type (PH) Distribution}
\description{
A function to generate an object of \code{\linkS4class{cf1}}.
}
\usage{
cf1(size, alpha, rate, class = "CsparseMatrix")
}
\arguments{
\item{size}{a value for the number of phases.}
\item{alpha}{a vector for the initial probabilities of PH distribution.}
\item{rate}{a vector for transition rates to next phase (diagonal elements of Q).}
\item{class}{name of Matrix class for \code{Q}.}
}
\details{
  The PH distribution with parameters \eqn{alpha}, \eqn{Q} and \eqn{xi}:
  Cumulative probability function; \deqn{F(q) = 1 - alpha \exp( Q q ) 1}
  Probability density function; \deqn{f(x) = alpha \exp( Q x ) xi},
  where \eqn{Q} is a bidiagonal matrix whose entries are sorted.
}
\note{
  \code{rph} is a generic function and is specified for \code{\linkS4class{cf1}}.
}
\value{
\code{cf1} gives an object of canonical form 1 that is a subclass of PH distribution.
}

\seealso{\code{\link{ph}}, \code{\link{herlang}}}

\examples{
## create a CF1 with 5 phases
(param1 <- cf1(5))

## create a CF1 with 5 phases
(param1 <- cf1(size=5))

## create a CF1 with specific parameters
(param2 <- cf1(alpha=c(1,0,0), rate=c(1.0,2.0,3.0)))

## p.d.f. for 0, 0.1, ..., 1
(dph(x=seq(0, 1, 0.1), ph=param2))

## c.d.f. for 0, 0.1, ..., 1
(pph(q=seq(0, 1, 0.1), ph=param2))

## generate 10 samples (this is quiker than rph with general ph)
(rph(n=10, ph=param2))

}
\keyword{distribution}

