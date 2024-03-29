% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_map.R
\name{map.acf}
\alias{map.acf}
\title{k-lag correlation of MAP}
\usage{
map.acf(map, ...)
}
\arguments{
\item{map}{An instance of MAP}

\item{...}{Others}
}
\value{
A vector of k-lag correlation
}
\description{
Compute k-lag correlation
}
\examples{
## create an MAP with specific parameters
(param1 <- map(alpha=c(1,0,0),
               D0=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-4)),
               D1=rbind(c(1,1,0),c(1,0,1),c(2,0,1))))

## create an ER-HMM with specific parameters
(param2 <- erhmm(shape=c(2,3), alpha=c(0.3,0.7),
                 rate=c(1.0,10.0),
                 P=rbind(c(0.3, 0.7), c(0.1, 0.9))))

map.acf(map=param1)
map.acf(map=param2)

}
