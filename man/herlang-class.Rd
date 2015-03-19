\name{herlang-class}
\docType{class}
\alias{herlang-class}

\title{Class of hyper Erlang}
\description{Parameters for a hyper Erlang.}
\section{Objects from the Class}{Objects are usually created by a \link{herlang}.}
\section{Slots}{
  \describe{
    \item{\code{size}:}{the number of components (hyper Erlang components).}
    \item{\code{mixrate}:}{a vector of mixed rates (probability for selecting a component).}
    \item{\code{shape}:}{shape parameters for Erlang distributions.}
    \item{\code{rate}:}{rate parameters for Erlang distributions.}
  }
}
\section{Methods}{
  \describe{
    \item{ph.moment}{\code{signature(ph = "herlang")}: ... }
    % \item{emfit.init}{\code{signature(model = "herlang", data = "phdata.wtime")}: ... }
    % \item{emfit.init}{\code{signature(model = "herlang", data = "phdata.group")}: ... }
    % \item{emfit.estep}{\code{signature(model = "herlang", data = "phdata.wtime")}: ... }
    % \item{emfit.estep}{\code{signature(model = "herlang", data = "phdata.group")}: ... }
    % \item{emfit.mstep}{\code{signature(model = "herlang")}: ... }
  }
}
%\references{}
%\author{}
\note{
  This class can be converted to \code{\linkS4class{ph}}.
}
\seealso{
  Classes \code{\linkS4class{ph}} and \code{\linkS4class{cf1}}.

}
\examples{
## create a hyper Erlang consisting of two Erlang
## with shape parameters 2 and 3.
(param1 <- herlang(c(2,3)))

## create a hyper Erlang consisting of two Erlang
## with shape parameters 2 and 3.
(param1 <- herlang(shape=c(2,3)))

## create a hyper Erlang with specific parameters
(param2 <- herlang(shape=c(2,3), mixrate=c(0.3,0.7),
                   rate=c(1.0,10.0)))

## convert to a general PH
as(param2, "ph")

## p.d.f. for 0, 0.1, ..., 1
(dherlang(x=seq(0, 1, 0.1), herlang=param2))

## c.d.f. for 0, 0.1, ..., 1
(pherlang(q=seq(0, 1, 0.1), herlang=param2))

## generate 10 samples
(rherlang(n=10, herlang=param2))

}
\keyword{classes}
