% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_ph.R
\name{ph}
\alias{ph}
\title{Create GPH distribution}
\usage{
ph(size, alpha, Q, xi)
}
\arguments{
\item{size}{An integer for the number of phases}

\item{alpha}{A vector of initial probability}

\item{Q}{An infinitesimal generator}

\item{xi}{An exit rate vector}
}
\value{
An instance of GPH
}
\description{
Create an instance of GPH
}
\note{
This function can omit several patterns of arguments. For example, \code{ph(5)}
omit the arguments \code{alpha}, \code{Q} and \code{xi}. In this case, the default values are
assigned to them.
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

}
