\name{erhmm}
\alias{erhmm}
\title{ER-HMM (HMM with Erlang outputs)}
\description{
A function to generate an object of \code{\linkS4class{erhmm}}.
}
\usage{
erhmm(shape, alpha, rate, P, class = "CsparseMatrix")
}
\arguments{
\item{shape}{an integer vector of shape parameters of Erlang outputs.}
\item{alpha}{a vector for initial probabilities of HMM states.}
\item{rate}{a vector of rate parameters of Erlang outputs.}
\item{P}{an object of Matrix class for a transition probability matrix of HMM.}
\item{class}{name of Matrix class for \code{P}.}
}
\details{
  ER-HMM has parameters \eqn{alpha}, \eqn{shape}, \eqn{rate} and \eqn{P}. HMM state chages according to a discrete-time Markov chain with transition matrix \eqn{P}. At each HMM state, there is an inherent Erlang distriution as an output. This model can be converted to an MAP.
}
\note{
  \code{erhmm} requires shape parameters. Other parameters have default values.
}
\value{
\code{erhmm} gives an object of ER-HMM.
}

\seealso{
\code{\link{map}}, \code{\link{gmmpp}}, \code{\link{map.mmoment}},
\code{\link{map.jmoment}}, \code{\link{map.acf}}
}

\examples{

## create an ER-HMM consisting of two Erlang components with
## shape parameters 2 and 3.
erhmm(c(2,3))

## create an ER-HMM consisting of two Erlang components with
## shape parameters 2 and 3.
erhmm(shape=c(2,3))

## create an ER-HMM with specific parameters
(param <- erhmm(shape=c(2,3), alpha=c(0.3,0.7),
                rate=c(1.0,10.0),
                P=rbind(c(0.3, 0.7), c(0.1, 0.9))))

## convert to a general MAP
as(param, "map")

## marginal moments of MAP
map.mmoment(k=3, map=as(param, "map"))

## joint moments of MAP
map.jmoment(lag=1, map=as(param, "map"))

## k-lag correlation
map.acf(map=as(param, "map"))

}

