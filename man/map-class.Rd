\name{map-class}
\docType{class}
\alias{map-class}
\alias{gmmpp-class}

\title{Classes of MAP}
\description{
Parameters for MAP and MMPP.
}
\section{Objects from the Class}{Objects are usually created by \link{map}, \link{mmpp} or \link{gmmpp}.}
\section{Slots}{
  \describe{
    \item{\code{size}:}{the number of phases (internal states).}
    \item{\code{alpha}:}{a probability (row) vector to decide an initial phase.}
    \item{\code{D0}:}{a square matrix that means transition rates without arrivals.}
    \item{\code{D1}:}{a square matrix that means transition rates with arrivals. In the case of MMPP, D1 should be a diagonal matrix.}
    \item{\code{df}:}{the number of free parameters.}
  }
}
% \section{Methods}{
%   \describe{
%     \item{ph.moment}{\code{signature(ph = "ph")}: ... }
%     % \item{emfit.init}{\code{signature(model = "ph")}: ... }
%     % \item{emfit.estep}{\code{signature(model = "ph", data = "phdata.wtime")}: ... }
%     % \item{emfit.estep}{\code{signature(model = "ph", data = "phdata.group")}: ... }
%     % \item{emfit.mstep}{\code{signature(model = "ph")}: ... }
%   }
% }
%\references{}
%\author{}
%\note{}
\seealso{
  Classes \code{\linkS4class{erhmm}}.
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
\keyword{classes}
