\name{phfit.3mom}
\alias{phfit.3mom}
\title{PH fitting with three moments}

\description{
  estimates PH parameters from three moments.
}
\usage{
phfit.3mom(m1, m2, m3, method = c("Osogami06", "Bobbio05"),
	max.phase = 50, epsilon = sqrt(.Machine$double.eps))
}
\arguments{
  \item{m1}{a value of the first moment.}
  \item{m2}{a value of the second moment.}
  \item{m3}{a value of the third moment.}
  \item{method}{the name of moment matching method.}
  \item{max.phase}{an integer for the maximum number of phases in the method "Osogami06".}
  \item{epsilon}{a value of precision in the method "Osogami06".}
}
\value{
  returns an object of S4 class of general PH \code{\linkS4class{ph}}.
}
\note{
  The method "Osogami06" checks the first three moements on whether there exists an PH whose three moements match to them. In such case, the method "Bobbio05" often returns an error.
}
\seealso{
  \code{\linkS4class{ph}}, \code{\link{ph.moment}}
}
\references{
Osogami, T. and Harchol-Balter, M. (2006)
Closed Form Solutions for Mapping General Distributions to Minimal PH Distributions.
\emph{Performance Evaluation},
\bold{63}(6), 524--552.

Bobbio, A., Horvath, A. and Telek, M. (2005)
Matching Three Moments with Minimal Acyclic Phase Type Distributions.
\emph{Stochastic Models},
\bold{21}(2-3), 303--326.
}
\examples{
## Three moment matching
## Moments of Weibull(shpape=2, scale=1); (0.886227, 1.0, 1.32934)
(result1 <- phfit.3mom(0.886227, 1.0, 1.32934))

## Three moment matching
## Moments of Weibull(shpape=2, scale=1); (0.886227, 1.0, 1.32934)
(result2 <- phfit.3mom(0.886227, 1.0, 1.32934, method="Bobbio05"))

## mean
ph.mean(result1)
ph.mean(result2)

## variance
ph.var(result1)
ph.var(result2)

## up to 5 moments 
ph.moment(5, result1)
ph.moment(5, result2)

}
