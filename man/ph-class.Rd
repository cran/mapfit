\name{ph-class}
\docType{class}
\alias{ph-class}

\title{Class of general PH distributions}
\description{
Parameters for a general PH distribution.
}
\section{Objects from the Class}{Objects are usually created by a \link{ph}.}
\section{Slots}{
  \describe{
    \item{\code{size}:}{the number of phases (transient states).}
    \item{\code{alpha}:}{a probability (row) vector to decide an initial phase.}
    \item{\code{Q}:}{a square matrix that means transition rates between phases.}
    \item{\code{xi}:}{a column vector for exiting rates from phases to an absorbing state.}
    \item{\code{df}:}{the number of free parameters.}
  }
}
\section{Methods}{
  \describe{
    \item{ph.moment}{\code{signature(ph = "ph")}: ... }
    % \item{emfit.init}{\code{signature(model = "ph")}: ... }
    % \item{emfit.estep}{\code{signature(model = "ph", data = "phdata.wtime")}: ... }
    % \item{emfit.estep}{\code{signature(model = "ph", data = "phdata.group")}: ... }
    % \item{emfit.mstep}{\code{signature(model = "ph")}: ... }
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  Classes \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}.
}
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
\keyword{classes}
