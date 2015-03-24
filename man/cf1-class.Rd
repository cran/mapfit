\name{cf1-class}
\docType{class}
\alias{cf1-class}

\title{Class of canonical form 1 (PH distribution)}
\description{
Parameters for a canonical form 1 which is a subclass of PH.
This is extended from \code{\linkS4class{ph}}.}
\section{Objects from the Class}{Objects are usually created by a \link{cf1}.}
\section{Slots}{
  \describe{
    \item{\code{rate}:}{transition rates to the next phase.}
  }
  The following slots are inherited from \code{\linkS4class{ph}}:
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
    \item{\code{\link{ph.moment}}}{\code{signature(ph = "cf1")}: ... }
    % \item{emfit.init}{\code{signature(model = "cf1")}: ... }
    % \item{emfit.mstep}{\code{signature(model = "cf1")}: ... }
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  Classes \code{\linkS4class{ph}} and \code{\linkS4class{herlang}}.

}
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
\keyword{classes}
