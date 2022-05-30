#' PH fitting with point data
#' 
#' Estimates PH parameters from point data.
#' 
#' @param ph An object of S4 class for MAP. The estimation algorithm is selected depending on thie class.
#' @param x A vector for point data.
#' @param weights A vector of weights for points.
#' @param method The name of estimation method for hyper Erlang (\code{\linkS4class{herlang}}).
#' @param lbound A value for lower limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).
#' @param ubound A value for upper limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).
#' @param control A list of parameters for controlling the fitting process.
#' @param verbose A list of parameters for displaying the fitting process.
#' @param ... Further arguments for methods.
#' @return
#' Returns a list with components, which is an object of S3 class \code{phfit.result};
#' \item{model}{an object for estimated PH class (\code{\linkS4class{ph}}, \code{\linkS4class{cf1}}, \code{\linkS4class{herlang}}).}
#' \item{llf}{a value of the maximum log-likelihood.}
#' \item{df}{a value of degrees of freedom of the model.}
#' \item{aic}{a value of Akaike information criterion.}
#' \item{iter}{the number of iterations.}
#' \item{convergence}{a logical value for the convergence of estimation algorithm.}
#' \item{ctime}{computation time (user time).}
#' \item{data}{an object for MAP data class}
#' \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
#' \item{rerror}{a value of relative error for llf at the last step of algorithm.}
#' \item{control}{a list of the argument of \code{control}.}
#' \item{verbose}{a list of the argument of \code{verbose}.}
#' \item{call}{the matched call.}
#' 
#' @seealso \code{\link{phfit.group}}, \code{\link{phfit.density}}, \code{\linkS4class{ph}}, \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}
#' @examples 
#' ## make sample
#' wsample <- rweibull(n=100, shape=2, scale=1)
#' 
#' ## PH fitting for general PH
#' (result1 <- phfit.point(ph=ph(2), x=wsample))
#' 
#' ## PH fitting for CF1
#' (result2 <- phfit.point(ph=cf1(2), x=wsample))
#' 
#' ## PH fitting for hyper Erlang
#' (result3 <- phfit.point(ph=herlang(3), x=wsample))
#' 
#' ## mean
#' ph.mean(result1$model)
#' ph.mean(result2$model)
#' ph.mean(result3$model)
#' 
#' ## variance
#' ph.var(result1$model)
#' ph.var(result2$model)
#' ph.var(result3$model)
#' 
#' ## up to 5 moments 
#' ph.moment(5, result1$model)
#' ph.moment(5, result2$model)
#' ph.moment(5, result3$model)
#' 
#' @export

phfit.point <- function(ph, x, weights, method = c("all", "increment"),
  lbound = 1, ubound = NULL, control = list(), verbose = list(), ...) {
  data <- phfit.time.data.frame(time=x, weights=weights)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
      phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
        control=control, verbose=verbose, ...)
    })
}

#' PH fitting with grouped data
#' 
#' Estimates PH parameters from grouped data.
#' 
#' @param ph An object of S4 class for MAP. The estimation algorithm is selected depending on thie class.
#' @param counts A vector of the number of points in intervals.
#' @param breaks A vector for a sequence of points of boundaries of intervals. This is equivalent to \code{c(0,cumsum(intervals))}. If this is missing, it is assigned to \code{0:length(counts)}.
#' @param intervals A vector of time lengths for intervals. This is equivalent to \code{diff(breaks)}). If this is missing, it is assigned to \code{rep(1,length(counts))}.
#' @param instant A vector of integers to indicate whether sample is drawn at the last of interval. If instant is 1, a sample is drawn at the last of interval. If instant is 0, no sample is drawn at the last of interval. By using instant, point data can be expressed by grouped data. If instant is missing, it is given by \code{rep(0L,length(counts))}, i.e., there are no sampels at the last of interval.
#' @param method The name of estimation method for hyper Erlang (\code{\linkS4class{herlang}}).
#' @param lbound A value for lower limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).
#' @param ubound A value for upper limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).
#' @param control A list of parameters for controlling the fitting process.
#' @param verbose A list of parameters for displaying the fitting process.
#' @param ... Rurther arguments for methods.
#' @return
#' Returns a list with components, which is an object of S3 class \code{phfit.result};
#' \item{model}{an object for estimated PH class (\code{\linkS4class{ph}}, \code{\linkS4class{cf1}}, \code{\linkS4class{herlang}}).}
#' \item{llf}{a value of the maximum log-likelihood.}
#' \item{df}{a value of degrees of freedom of the model.}
#' \item{aic}{a value of Akaike information criterion.}
#' \item{iter}{the number of iterations.}
#' \item{convergence}{a logical value for the convergence of estimation algorithm.}
#' \item{ctime}{computation time (user time).}
#' \item{data}{an object for MAP data class}
#' \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
#' \item{rerror}{a value of relative error for llf at the last step of algorithm.}
#' \item{control}{a list of the argument of \code{control}.}
#' \item{verbose}{a list of the argument of \code{verbose}.}
#' \item{call}{the matched call.}
#' 
#' @note
#' In this method, we can handle truncated data using \code{NA} and \code{Inf};
#' \code{phfit.group(ph=cf1(5), counts=c(countsdata, NA), breaks=c(breakdata, +Inf))}
#' \code{NA} means missing of count data at the conrredponding interval, and \code{Inf} ia allowed to put 
#' the last of breaks or intervals which represents a special interval [the last break point,infinity).
#' 
#' @seealso \code{\link{phfit.point}}, \code{\link{phfit.density}}, \code{\linkS4class{ph}}, \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}
#' @examples
#' ## make sample
#' wsample <- rweibull(n=100, shape=2, scale=1)
#' wgroup <- hist(x=wsample, breaks="fd", plot=FALSE)
#' 
#' ## PH fitting for general PH
#' (result1 <- phfit.group(ph=ph(2), counts=wgroup$counts, breaks=wgroup$breaks))
#' 
#' ## PH fitting for CF1
#' (result2 <- phfit.group(ph=cf1(2), counts=wgroup$counts, breaks=wgroup$breaks))
#' 
#' ## PH fitting for hyper Erlang
#' (result3 <- phfit.group(ph=herlang(3), counts=wgroup$counts, breaks=wgroup$breaks))
#' 
#' ## mean
#' ph.mean(result1$model)
#' ph.mean(result2$model)
#' ph.mean(result3$model)
#' 
#' ## variance
#' ph.var(result1$model)
#' ph.var(result2$model)
#' ph.var(result3$model)
#' 
#' ## up to 5 moments 
#' ph.moment(5, result1$model)
#' ph.moment(5, result2$model)
#' ph.moment(5, result3$model)
#' 
#' @export

phfit.group <- function(ph, counts, breaks, intervals, instant,
 method = c("all", "increment"), lbound = 1, ubound = NULL, control = list(), verbose = list(), ...) {
  data <- phfit.group.data.frame(counts=counts, breaks=breaks, difftime=intervals, instant=instant)
  switch(class(ph),
    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
    "herlang"={
      phsize <- sum(ph@shape)
      if (is.null(ubound)) {
        ubound <- phsize
      }
    	phfit.herlang(phsize=phsize, data=data, method=method, lbound=lbound, ubound=ubound,
    		control=control, verbose=verbose, ...)
    })
}

#' PH fitting with density function
#' 
#' Estimates PH parameters from density function.
#' 
#' @param ph An object of S4 class for MAP. The estimation algorithm is selected depending on thie class.
#' @param f A faunction object for a density function.
#' @param method The name of estimation method for hyper Erlang (\code{\linkS4class{herlang}}).
#' @param lbound A value for lower limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).
#' @param ubound A value for upper limit for the number of states in hyper Erlang (\code{\linkS4class{herlang}}).
#' @param deformula An object for formulas of numerical integration. It is not necessary to change it when the density function is defined on the positive domain [0,infinity).
#' @param weight.zero A absolute value which is regarded as zero in numerical integration.
#' @param weight.reltol A value for precision of numerical integration.
#' @param start.divisions A value for starting value of divitions in deformula.
#' @param max.iter A value for the maximum number of iterations to increase divisions in deformula.
#' @param control A list of parameters for controlling the fitting process.
#' @param verbose A list of parameters for displaying the fitting process.
#' @param ... Further arguments for methods, which are also used to send the arguments to density function.
#' 
#' @return
#' Returns a list with components, which is an object of S3 class \code{phfit.result};
#' \item{model}{an object for estimated PH class (\code{\linkS4class{ph}}, \code{\linkS4class{cf1}}, \code{\linkS4class{herlang}}).}
#' \item{llf}{a value of the maximum log-likelihood (a netative value of the cross entropy).}
#' \item{df}{a value of degrees of freedom of the model.}
#' \item{aic}{a value of Akaike information criterion (this is not meaningless in this case).}
#' \item{KL}{a value of Kullback-Leibler divergence.}
#' \item{iter}{the number of iterations.}
#' \item{convergence}{a logical value for the convergence of estimation algorithm.}
#' \item{ctime}{computation time (user time).}
#' \item{data}{an object for MAP data class}
#' \item{aerror}{a value of absolute error for llf at the last step of algorithm.}
#' \item{rerror}{a value of relative error for llf at the last step of algorithm.}
#' \item{control}{a list of the argument of \code{control}.}
#' \item{verbose}{a list of the argument of \code{verbose}.}
#' \item{call}{the matched call.}
#' 
#' @note
#' Any of density function can be applied to the argument \code{f}, where \code{f} should be defined \code{f <- function(x, ...)}. The first argument of \code{f} should be an integral parameter. The other parameters are set in the argument \code{...} of \code{phfit.density}. The truncated density function can also be used directly.
#' 
#' @seealso \code{\link{phfit.point}}, \code{\link{phfit.group}}, \code{\linkS4class{ph}}, \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}
#' @examples
#' ####################
#' ##### truncated density
#' ####################
#' 
#' ## PH fitting for general PH
#' (result1 <- phfit.density(ph=ph(2), f=dnorm, mean=3, sd=1))
#' 
#' ## PH fitting for CF1
#' (result2 <- phfit.density(ph=cf1(2), f=dnorm, mean=3, sd=1))
#' 
#' ## PH fitting for hyper Erlang
#' (result3 <- phfit.density(ph=herlang(3), f=dnorm, mean=3, sd=1))
#' 
#' ## mean
#' ph.mean(result1$model)
#' ph.mean(result2$model)
#' ph.mean(result3$model)
#' 
#' ## variance
#' ph.var(result1$model)
#' ph.var(result2$model)
#' ph.var(result3$model)
#' 
#' ## up to 5 moments 
#' ph.moment(5, result1$model)
#' ph.moment(5, result2$model)
#' ph.moment(5, result3$model)
#' 
#' @export

phfit.density <- function(ph, f,
    method = c("all", "increment"), lbound = 1, ubound = NULL, 
    deformula = deformula.zeroinf, weight.zero = 1.0e-12,
    weight.reltol = 1.0e-8, start.divisions = 8, max.iter = 12,
    control = list(), verbose = list(), ...) {
      x <- deformula(f, ..., zero.eps = weight.zero,
                     rel.tol = weight.reltol,
                     start.divisions = start.divisions, max.iter = max.iter)
      ll <- x$h * sum(x$w * log(f(x$x, ...)))
      data <- phfit.time.data.frame(time=x$x, weights=x$w)
      res <- switch(class(ph),
                    "ph"=phfit.gen(ph=ph, data=data, control=control, verbose=verbose, ...),
                    "cf1"=phfit.cf1(ph=ph, data=data, control=control, verbose=verbose, ...),
                    "herlang"={
                      phsize <- sum(ph@shape)
                      if (is.null(ubound)) {
                        ubound <- phsize
                      }
                      phfit.herlang(phsize=phsize, data=data, method=method,
                                    lbound=lbound, ubound=ubound, control=control, verbose=verbose, ...)
                      })
      res <- c(res, list(KL=ll - res$llf * x$h))
      class(res) <- "phfit.result"
      res
}

#' @aliases phfit.point phfit.group phfit.density
#' @export

print.phfit.result <- function (x, ...) {
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
