#' Class of canonical form 1 (PH distribution)
#' 
#' Parameters for a canonical form 1 which is a subclass of PH.
#' This is extended from \code{\linkS4class{ph}}.
#' 
#' @docType class
#' @name cf1-class
#' @slot rate Transition rates to the next phase.
#' @slot size The number of phases (transient states).
#' @slot alpha A probability (row) vector to decide an initial phase.
#' @slot Q A square matrix that means transition rates between phases.
#' @slot xi A column vector for exiting rates from phases to an absorbing state.
#' @slot df The number of free parameters.
#' 
#' @note 
#' Objects are usually created by a \link{cf1}.
#' 
#' Methods:
#' \describe{
#' \item{\code{\link{ph.moment}}}{\code{signature(ph = "cf1")}: ... }
#' \item{emfit.init}{\code{signature(model = "cf1")}: ... }
#' \item{emfit.mstep}{\code{signature(model = "cf1")}: ... }
#' }
#' 
#' @seealso Classes \code{\linkS4class{ph}} and \code{\linkS4class{herlang}}.
#' 
#' @examples 
#' ## create a CF1 with 5 phases
#' (param1 <- cf1(5))
#' 
#' ## create a CF1 with 5 phases
#' (param1 <- cf1(size=5))
#' 
#' ## create a CF1 with specific parameters
#' (param2 <- cf1(alpha=c(1,0,0), rate=c(1.0,2.0,3.0)))
#' 
#' ## p.d.f. for 0, 0.1, ..., 1
#' (dph(x=seq(0, 1, 0.1), ph=param2))
#' 
#' ## c.d.f. for 0, 0.1, ..., 1
#' (pph(q=seq(0, 1, 0.1), ph=param2))
#' 
#' ## generate 10 samples (this is quiker than rph with general ph)
#' (rph(n=10, ph=param2))
#' 
#' @keywords classes
#' @exportClass cf1
NULL

#' Canonical Form 1 for Phase-Type (PH) Distribution
#' 
#' A function to generate an object of \code{\linkS4class{cf1}}.
#' 
#' - The PH distribution with parameters \eqn{\alpha}, \eqn{Q} and \eqn{\xi = - Q 1}:
#' - Cumulative probability function; \deqn{F(q) = 1 - \alpha \exp( Q q ) 1}
#' - Probability density function; \deqn{f(x) = \alpha \exp( Q x ) \xi,}
#' where \eqn{Q} is a bidiagonal matrix whose entries are sorted.
#'
#' @param size A value for the number of phases.
#' @param alpha A vector for the initial probabilities of PH distribution.
#' @param rate A vector for transition rates to next phase (diagonal elements of Q).
#' @param class Name of Matrix class for \code{Q}.
#' @return \code{cf1} gives an object of canonical form 1 that is a subclass of PH distribution.
#' @seealso \code{\link{ph}}, \code{\link{herlang}}
#' @examples
#' ## create a CF1 with 5 phases
#' (param1 <- cf1(5))
#' 
#' ## create a CF1 with 5 phases
#' (param1 <- cf1(size=5))
#' 
#' ## create a CF1 with specific parameters
#' (param2 <- cf1(alpha=c(1,0,0), rate=c(1.0,2.0,3.0)))
#' 
#' ## p.d.f. for 0, 0.1, ..., 1
#' (dph(x=seq(0, 1, 0.1), ph=param2))
#' 
#' ## c.d.f. for 0, 0.1, ..., 1
#' (pph(q=seq(0, 1, 0.1), ph=param2))
#' 
#' ## generate 10 samples (this is quiker than rph with general ph)
#' (rph(n=10, ph=param2))
#' 
#' @keywords distribution
#' @export

cf1 <- function(size, alpha, rate, class="CsparseMatrix") {
  if (missing(size)) {
    if (missing(alpha) || missing(rate)) {
      stop("alpha and rate are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(rate)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (!missing(alpha) || !missing(rate)) {
        warning("alpha and rate are ignored.")
      }
      alpha <- rep(1.0/size, size)
      rate <- rep(1.0, size)      
    }
  }

  xi <- numeric(size)
  xi[size] <- rate[size]
  if (size >= 2) {
    i <- c(1:size, 1:(size-1))
    j <- c(1:size, 2:size)
    x <- c(-rate, rate[1:(size-1)])
    Q <- sparseMatrix(i=i, j=j, x=x)
  } else {
    Q <- matrix(-rate[1],1,1)
  }
  if (!is(Q, class)) {
    Q <- as(Q, class)
  }
  zero <- 1.0e-8
  df <- sum(abs(alpha) > zero) - 1 + sum(abs(rate) > zero)
  new("cf1", size=size, alpha=alpha, Q=Q, xi=xi, rate=rate, df=df)
}

cf1.sample <- function(n, ph) {
  res <- rexp(n, ph@rate[ph@size])
  if (ph@size == 1) {
    return(res)
  }
  x <- cumsum(rmultinom(n=1, size=n, prob=ph@alpha))
  for (l in (ph@size-1):1) {
    y <- x[l]
    if (y == 0) break
    res[1:y] <- res[1:y] + rexp(y, ph@rate[l])
  }
  sample(res)
}

cf1.param <- function(size, data,
  diff.init = c(1, 4, 16, 64, 256, 1024),
  scale.init = c(0.5, 1.0, 2.0), maxiter.init = 5, verbose, class) {
  if (verbose$cf1init) cat("Initializing CF1 ...\n")
  m <- mapfit.mean(data)
  maxllf <- -Inf
  for (x in scale.init) {
    for (s in diff.init) {
      ph <- cf1.param.linear(size, m * x, s, class)
      phres <- try(emfit(model=ph, data=data,
        initialize = FALSE, control = list(maxiter=maxiter.init)), silent = TRUE)
      if (!is(phres, "try-error")) {
        if (is.finite(phres$llf)) {
          if (maxllf < phres$llf) {
            maxllf <- phres$llf
            maxph <- ph
            if (verbose$cf1init) cat("o")
          }
          else {
            if (verbose$cf1init) cat("x")
          }
        }
        else {
          if (verbose$cf1init) cat("-")
        }
      }
      else {
        if (verbose$cf1init) cat("-")
      }
    }
    if (verbose$cf1init) cat("\n")
    for (s in diff.init) {
      ph <- cf1.param.power(size, m * x, s, class)
      phres <- try(emfit(model=ph, data=data,
        initialize = FALSE, control = list(maxiter=maxiter.init)), silent = TRUE)
      if (!is(phres, "try-error")) {
        if (is.finite(phres$llf)) {
          if (maxllf < phres$llf) {
            maxllf <- phres$llf
            maxph <- ph
           if (verbose$cf1init) cat("o")
          }
          else {
            if (verbose$cf1init) cat("x")
          }
        }
        else {
          if (verbose$cf1init) cat("-")
        }
      }
      else {
        if (verbose$cf1init) cat("-")
      }
    }
    if (verbose$cf1init) cat("\n")
  }
##  if (verbose$cf1init) cat("done\n")
  maxph
}

cf1.param.power <- function(size, mean, s, class) {
  alpha <- rep(1.0/size, size)
  rate <- numeric(size)

  p <- exp(1.0/(size-1) * log(s))
  total <- 1.0
  tmp <- 1.0
  for (i in 2:size) {
    tmp <- tmp * i / ((i-1) * p)
    total <- total + tmp
  }
  base <- total / (size * mean)
  tmp <- base
  for (i in 1:size) {
    rate[i] <- tmp
    tmp <- tmp * p
  }
  cf1(alpha=alpha, rate=rate, class=class)
}

cf1.param.linear <- function(size, mean, s, class) {
  alpha <- rep(1.0/size, size)
  rate <- numeric(size)

  al <- (s-1)/(size-1)
  total <- 1.0
  for (i in 2:size) {
    total <- total + i / (al * (i-1) + 1)
  }
  base <- total / (size * mean)
  for (i in 1:size) {
    tmp <- base * (al * (i-1) + 1)
    rate[i] <- tmp
  }
  cf1(alpha=alpha, rate=rate, class=class)
}

setMethod("emfit.print", signature(model = "cf1"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Initial : ", model@alpha, "\n")
    cat("Rate    : ", model@rate, "\n")
  }
)

setMethod("emfit.init", signature(model = "cf1"),
  function(model, data, verbose=list(), ...) {
    cf1.param(size=model@size, data=data, verbose=verbose, class=class(model@Q), ...)
  }
)

#### mstep

setMethod("emfit.mstep", signature(model = "cf1"),
  function(model, eres, data, ...) {
    res <- .Call('phfit_mstep_cf1', PACKAGE='mapfit', model, eres, data)
    model@alpha <- res[[1]]
    model@xi <- res[[2]]
    model@Q@x <- res[[3]]
    model@rate <- -diag(model@Q)
    model
  })

