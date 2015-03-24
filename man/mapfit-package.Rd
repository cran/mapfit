\name{mapfit-package}
\alias{mapfit-package}
\alias{mapfit}
\alias{phfit}
\docType{package}
\title{PH/MAP parameter estimation tool}
\description{
Estimation methods for phase-type
  distribution (PH) and Markovian arrival process (MAP) from
  empirical data (point and grouped data) and density function.
}
\details{
\tabular{ll}{
Package: \tab mapfit\cr
Type: \tab Package\cr
Version: \tab 0.9.7\cr
Date: \tab 2015-03-24\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Hiroyuki Okamura

Maintainer: Hiroyuki Okamura <okamu@rel.hiroshima-u.ac.jp>
}
% \references{
% ~~ Literature or other references for background information ~~
% }
\keyword{ package }
\seealso{
\code{\link{phfit.point}}
\code{\link{phfit.group}}
\code{\link{phfit.density}}
\code{\link{phfit.3mom}}
\code{\link{mapfit.point}}
\code{\link{mapfit.group}}
}
\examples{

### PH fitting with grouped data
## make sample
wgroup <- hist(x=rweibull(n=100, shape=2, scale=1),
	           breaks="fd", plot=FALSE)

## PH fitting for CF1
phfit.group(ph=cf1(2), counts=wgroup$counts,
            breaks=wgroup$breaks)

## PH fitting for hyper Erlang
phfit.group(ph=herlang(3), counts=wgroup$counts,
            breaks=wgroup$breaks)

### MAP fitting with point data
data(BCpAug89)
BCpAug89s <- head(BCpAug89, 50)

## MAP fitting for ER-HMM (fast estimation algorithm)
mapfit.point(map=erhmm(3), x=cumsum(BCpAug89s))


### MAP fitting with grouped data
## make grouped data
BCpAug89.group <- hist(cumsum(BCpAug89s),
                       breaks=seq(0, 0.15, 0.005),
                       plot=FALSE)

## MAP fitting with approximate MMPP
mapfit.group(map=gmmpp(2),
             counts=BCpAug89.group$counts,
             breaks=BCpAug89.group$breaks)

}
