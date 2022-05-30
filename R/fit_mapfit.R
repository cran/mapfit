#' MAP fitting with time point data
#' 
#' Estimates MAP parameters from time point data.
#' 
#' @param map An object of S4 class for MAP. The estimation algorithm is selected depending on thie class.
#' @param x A vector for time sequence of arrivals. This is equivalent to \code{cumsum(intervals)}. Either time or difftime should be given.
#' @param intervals A vector for the data for intrarrival time. This is equivalent to \code{diff(c(0,x)}). Either time or difftime should be given.
#' @param stationary A logical value that determine whether initial probability is given by a stationary vector of underlying Markov process or not.
#' @param method The name of estimation method for ER-HMM (\code{\linkS4class{erhmm}}).
#' @param lbound A value for lower limit for the number of states in ER-HMM (\code{\linkS4class{erhmm}}).
#' @param ubound A value for upper limit for the number of states in ER-HMM (\code{\linkS4class{erhmm}}).
#' @param control A list of parameters for controlling the fitting process.
#' @param verbose A list of parameters for displaying the fitting process.
#' @param ... Further arguments for methods.
#' @return
#' Returns a list with components, which is an object of S3 class \code{mapfit.result};
#' \item{model}{an object for estimated MAP class (\code{\linkS4class{map}}, \code{\linkS4class{erhmm}}).}
#' \item{llf}{a value of the maximum log-likelihood.}
#' \item{df}{a value of degrees of freedom of the model.}
#' \item{aic}{a value of Akaike information criterion.}
#' \item{iter}{the number of iterations.}
#' \item{convergence}{a logical value for the convergence of estimation algorithm.}
#' \item{ctime}{computation time (user time).}
#' \item{stationary}{a logical value for the argument \code{stationary}.}
#' \item{data}{an object for MAP data class}
#' \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
#' \item{rerror}{a value of relative error for llf at the last step of algorithm.}
#' \item{control}{a list of the argument of \code{control}.}
#' \item{verbose}{a list of the argument of \code{verbose}.}
#' \item{call}{the matched call.}
#' 
#' @seealso \code{\link{mapfit.group}}, \code{\linkS4class{map}} and \code{\linkS4class{erhmm}}
#' 
#' @examples 
#' ## load trace data
#' data(BCpAug89)
#' BCpAug89s <- head(BCpAug89, 50)
#' 
#' ## MAP fitting for general MAP
#' (result1 <- mapfit.point(map=map(2), x=cumsum(BCpAug89s)))
#'
#' ## MAP fitting for MMPP
#' (result2 <- mapfit.point(map=mmpp(2), x=cumsum(BCpAug89s)))
#' 
#' ## MAP fitting for ER-HMM
#' (result3 <- mapfit.point(map=erhmm(3), x=cumsum(BCpAug89s)))
#' 
#' ## marginal moments for estimated MAP
#' map.mmoment(k=3, map=result1$model)
#' map.mmoment(k=3, map=result2$model)
#' map.mmoment(k=3, map=as(result3$model, "map"))
#' 
#' ## joint moments for estimated MAP
#' map.jmoment(lag=1, map=result1$model)
#' map.jmoment(lag=1, map=result2$model)
#' map.jmoment(lag=1, map=as(result3$model, "map"))
#' 
#' ## lag-k correlation
#' map.acf(map=result1$model)
#' map.acf(map=result2$model)
#' map.acf(map=as(result3$model, "map"))
#' 
#' @export

mapfit.point <- function(map, x, intervals, stationary = TRUE,
  method = c("all", "increment"), lbound = 1, ubound = NULL,
  control = list(), verbose = list(), ...) {
  data <- mapfit.time.data.frame(x, intervals)
  switch(class(map),
    "map"=mapfit.gen(map=map, data=data, stationary=stationary, control=control, verbose=verbose, ...),
    "erhmm"={
      phsize <- sum(map@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      mapfit.erhmm(phsize=phsize, data=data, method=method,
        lbound=lbound, ubound=ubound, stationary=stationary, control=control, verbose=verbose, ...)
    })
}

#' MAP fitting with grouped data
#' 
#' Estimates MAP parameters from grouped data.
#' 
#' @param map S4 class for MAP. The estimation algorithm is selected depending on thie class.
#' @param counts A vector for the number of arrivals in time interval.
#' @param breaks A vector for time sequence to determine time interval. This is equivalent to \code{c(0,cumsum(intervals))}. If this is missing, it is assigned to \code{0:length(counts)}.
#' @param intervals A vector for a sequence of time length for intervals. This is equivalent to \code{diff(breaks)}). If this is missing, it is assigned to \code{rep(1,length(counts))}.
#' @param instant A vector of integer to indicate whether an arrival occurs at the last time of interval. If instant is 1, an arrival occurs at the last time of interval. If instant is 0, no arrival occurs at the last time of interval. By using instant, time point data can be expressed by grouped data class. If instant is missing, it is given by \code{rep(0,length(counts))}, i.e., there are no arrivals at the end of interval.
#' @param stationary A logical value that determine whether initial probability is given by a stationary vector of underlying Markov process or not.
#' @param control A list of parameters for controlling the fitting process.
#' @param verbose A list of parameters for displaying the fitting process.
#' @param ... Further arguments for methods.
#' @return
#' Returns a list with components, which is an object of S3 class \code{mapfit.result};
#' \item{model}{an object for estimated MAP class (\code{\linkS4class{map}}, \code{\linkS4class{erhmm}}).}
#' \item{llf}{a value of the maximum log-likelihood.}
#' \item{df}{a value of degrees of freedom of the model.}
#' \item{aic}{a value of Akaike information criterion.}
#' \item{iter}{the number of iterations.}
#' \item{convergence}{a logical value for the convergence of estimation algorithm.}
#' \item{ctime}{computation time (user time).}
#' \item{stationary}{a logical value for the argument \code{stationary}.}
#' \item{data}{an object for MAP data class}
#' \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
#' \item{rerror}{a value of relative error for llf at the last step of algorithm.}
#' \item{control}{a list of the argument of \code{control}.}
#' \item{verbose}{a list of the argument of \code{verbose}.}
#' \item{call}{the matched call.}
#' 
#' @seealso \code{\link{mapfit.point}}, \code{\linkS4class{map}} and \code{\linkS4class{gmmpp}}
#' 
#' @examples 
#' 
#' ## load trace data
#' data(BCpAug89)
#' BCpAug89s <- head(BCpAug89, 50)
#' 
#' ## make grouped data
#' BCpAug89.group <- hist(cumsum(BCpAug89s),
#'                          breaks=seq(0, 0.15, 0.005),
#'                          plot=FALSE)
#'                          
#' ## MAP fitting for general MAP
#' (result1 <- mapfit.group(map=map(2),
#'                         counts=BCpAug89.group$counts,
#'                         breaks=BCpAug89.group$breaks))
#' ## MAP fitting for MMPP
#' (result2 <- mapfit.group(map=mmpp(2),
#'                          counts=BCpAug89.group$counts,
#'                          breaks=BCpAug89.group$breaks))
#'                          
#' ## MAP fitting with approximate MMPP
#' (result3 <- mapfit.group(map=gmmpp(2),
#'                          counts=BCpAug89.group$counts,
#'                          breaks=BCpAug89.group$breaks))
#'
#' ## marginal moments for estimated MAP
#' map.mmoment(k=3, map=result1$model)
#' map.mmoment(k=3, map=result2$model)
#' map.mmoment(k=3, map=result3$model)
#' 
#' ## joint moments for estimated MAP
#' map.jmoment(lag=1, map=result1$model)
#' map.jmoment(lag=1, map=result2$model)
#' map.jmoment(lag=1, map=result3$model)
#' 
#' ## lag-k correlation
#' map.acf(map=result1$model)
#' map.acf(map=result2$model)
#' map.acf(map=result3$model)
#' 
#' @export

mapfit.group <- function(map, counts, breaks, intervals, instant, stationary = TRUE,
  control = list(), verbose = list(), ...) {
  data <- mapfit.group.data.frame(counts, breaks, intervals, instant)
  switch(class(map),
    "map"=mapfit.gen(map=map, data=data, stationary=stationary, control=control, verbose=verbose, ...),
    "gmmpp"=mapfit.gmmpp(map=map, data=data, stationary=stationary, control=control, verbose=verbose, ...)
  )
}

#' @aliases mapfit.point mapfit.group
#' @export

print.mapfit.result <- function (x, ...) {
  cat("\n")
  cat(sprintf("Maximum LLF: %f\n", x$llf))
  cat(sprintf("AIC: %f\n", x$aic))
  cat(sprintf("Iteration:  %d / %d\n", x$iter, x$control$maxiter))
  cat(sprintf("Computation time (user): %f\n", x$ctime))
  cat(sprintf("Convergence: %s\n", x$convergence))
  cat(sprintf("Error (abs): %e (tolerance %e)\n", x$aerror, x$control$abstol))
  cat(sprintf("Error (rel): %e (tolerance %e)\n", x$rerror, x$control$reltol))
  cat("\n")
  emfit.print(x$model)
  cat("\n\n")
  invisible(x)
}
