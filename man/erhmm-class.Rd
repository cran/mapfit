\name{erhmm-class}
\docType{class}
\alias{erhmm-class}

\title{Class of ER-HMM}
\description{Parameters for an ER-HMM (Hidden Markov Model with Erlang outputs).}
\section{Objects from the Class}{Objects are usually created by an \link{erhmm}.}
\section{Slots}{
  \describe{
    \item{\code{size}:}{the number of HMM states.}
    \item{\code{alpha}:}{a vector of initial probabilities for HMM states.}
    \item{\code{shape}:}{shape parameters for Erlang distributions. The sum of shape parameters is the number of phases of MAP.}
    \item{\code{rate}:}{rate parameters for Erlang distributions.}
    \item{\code{P}:}{an object of Matrix class for a transition probability matrix of HMM.}
  }
}
% \section{Methods}{
%   \describe{
%     \item{ph.moment}{\code{signature(ph = "herlang")}: ... }
%     % \item{emfit.init}{\code{signature(model = "herlang", data = "phdata.wtime")}: ... }
%     % \item{emfit.init}{\code{signature(model = "herlang", data = "phdata.group")}: ... }
%     % \item{emfit.estep}{\code{signature(model = "herlang", data = "phdata.wtime")}: ... }
%     % \item{emfit.estep}{\code{signature(model = "herlang", data = "phdata.group")}: ... }
%     % \item{emfit.mstep}{\code{signature(model = "herlang")}: ... }
%   }
% }
%\references{}
%\author{}
\note{
  This class can be converted to \code{\linkS4class{map}}.
}
\seealso{
  Classes \code{\linkS4class{map}} and \code{\linkS4class{gmmpp}}.

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
\keyword{classes}
