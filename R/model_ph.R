#' Class of general PH distributions
#' 
#' Parameters for a general PH distribution.
#'
#' @docType class
#' @name ph-class
#' @slot size The number of phases (transient states).
#' @slot alpha A probability (row) vector to decide an initial phase.
#' @slot Q A square matrix that means transition rates between phases.
#' @slot xi A column vector for exiting rates from phases to an absorbing state.
#' @slot df The number of free parameters.
#' 
#' @note 
#' Objects are usually created by a \link{ph}.
#' 
#' The methods of this class:
#' \describe{
#' \item{ph.moment}{\code{signature(ph = "ph")}: ... }
#' \item{emfit.init}{\code{signature(model = "ph")}: ... }
#' \item{emfit.estep}{\code{signature(model = "ph", data = "phdata.wtime")}: ... }
#' \item{emfit.estep}{\code{signature(model = "ph", data = "phdata.group")}: ... }
#' \item{emfit.mstep}{\code{signature(model = "ph")}: ... }
#' }
#' 
#' @seealso Classes \code{\linkS4class{cf1}} and \code{\linkS4class{herlang}}.
#' 
#' @examples
#' ## create a PH (full matrix) with 5 phases
#' (param1 <- ph(5))
#' 
#' ## create a PH (full matrix) with 5 phases
#' (param1 <- ph(size=5))
#' 
#' ## create a PH with specific parameters
#' (param2 <- ph(alpha=c(1,0,0),
#'               Q=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-1)),
#'               xi=c(2,2,0)))
#'
#' ## p.d.f. for 0, 0.1, ..., 1
#' (dph(x=seq(0, 1, 0.1), ph=param2))
#' 
#' ## c.d.f. for 0, 0.1, ..., 1
#' (pph(q=seq(0, 1, 0.1), ph=param2))
#' 
#' ## generate 10 samples
#' (rph(n=10, ph=param2))
#' 
#' @keywords classes
#' @exportClass ph
NULL

#' Phase-Type (PH) Distribution
#' 
#' Density function, distribution function and
#' random generation for the PH distribution, and
#' a function to generate an object of \code{\linkS4class{ph}}.
#' 
#' @aliases dph pph rph
#'
#' 
#' @param size A value for the number of phases.
#' @param alpha A vector for the initial probabilities of PH distribution.
#' @param Q An object of Matrix class for the initesmal generator of PH distribution.
#' @param xi A vector for the exit rates of PH distribution.
#' @param class Name of Matrix class for \code{Q}.
#' @param x Vectors of quantiles.
#' @param q Vectors of quantiles.
#' @param p A vector of probabilities.
#' @param n Number of observations.
#' @param ph An object of S4 class of PH (\code{\linkS4class{ph}}).
#' @param log Logical; if \code{TRUE}, the log density is returned.
#' @param lower.tail Logical; if \code{TRUE}, probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log.p Logical; if \code{TRUE}, the log probability is returned.
#' @return
#' \code{ph} gives an object of general PH distribution.
#' \code{dph} gives the density function, \code{pph} gives the distribution function,
#' and \code{rph} generates random samples.
#' 
#' @details 
#' The PH distribution with parameters \eqn{\alpha}, \eqn{Q} and \eqn{\xi}:
#' Cumulative probability function; \deqn{F(q) = 1 - \alpha \exp( Q q ) 1}
#' Probability density function; \deqn{f(x) = \alpha \exp( Q x ) \xi}
#' 
#' @note
#' \code{ph} requires either \code{size} or (\code{alpha}, \code{Q}, \code{xi}).
#' \code{rph} for \code{\linkS4class{ph}} is too slow. It is recommended to use
#' \code{rph} for \code{\linkS4class{cf1}}.
#' 
#' @seealso \code{\link{cf1}}, \code{\link{herlang}}
#' @examples
#' ## create a PH (full matrix) with 5 phases
#' (param1 <- ph(5))
#' 
#' ## create a PH (full matrix) with 5 phases
#' (param1 <- ph(size=5))
#' 
#' ## create a PH with specific parameters
#' (param2 <- ph(alpha=c(1,0,0), 
#'               Q=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-1)),
#'               xi=c(2,2,0))) 
#'
#' ## p.d.f. for 0, 0.1, ..., 1
#' (dph(x=seq(0, 1, 0.1), ph=param2))
#' 
#' ## c.d.f. for 0, 0.1, ..., 1
#' (pph(q=seq(0, 1, 0.1), ph=param2))
#' 
#' ## generate 10 samples
#' (rph(n=10, ph=param2))
#' @keywords distribution
#' @export

ph <- function(size, alpha, Q, xi, class="CsparseMatrix") {
  if (missing(size)) {
    if (missing(alpha) || missing(Q) || missing(xi)) {
      stop("alpha, Q and xi are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(Q) && !missing(xi)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (!missing(alpha) || !missing(Q) || !missing(xi)) {
        warning("alpha, Q and xi are ignored.")
      }
      alpha <- rep(1.0/size, size)
      Q <- matrix(1.0, size, size)
      diag(Q) <- rep(-size, size)
      xi <- rep(1.0, size)
    }
  }
  if (!is(Q, class)) {
    Q <- as(Q, class)
  }
  if (missing(xi)) {
    xi = -apply(Q, 1, sum)
  }
  zero <- 1.0e-8
  df <- sum(abs(alpha) > zero) - 1 + sum(abs(Q) > zero) +
    sum(abs(xi) > zero) - size + sum(abs(Matrix::diag(Q)) < zero)
  new("ph", size=size, alpha=alpha, Q=Q, xi=xi, df=df)
}

ph.bidiag <- function(size, class="CsparseMatrix") {
  if (size <= 1) {
    ph(size, class=class)
  } else {
  alpha <- rep(1/size,size)
  xi <- rep(0, size)
  Q <- matrix(0, size, size)
  for (i in 1:(size-1)) {
    Q[i,i] <- -1
    Q[i,i+1] <- 1
  }
  Q[size,size] <- -1
  xi[size] <- 1
  ph(alpha=alpha, Q=Q, xi=xi, class=class)
  }
}

ph.tridiag <- function(size, class="CsparseMatrix") {
  if (size <= 2) {
    ph(size, class=class)
  } else {
  alpha <- rep(1/size,size)
  xi <- rep(0, size)
  Q <- matrix(0, size, size)
  Q[1,1] <- -1
  Q[1,2] <- 1
  for (i in 2:(size-1)) {
    Q[i,i] <- -2
    Q[i,i-1] <- 1
    Q[i,i+1] <- 1
  }
  Q[size,size-1] <- 1
  Q[size,size] <- -2
  xi[size] <- 1
  ph(alpha=alpha, Q=Q, xi=xi, class=class)
  }
}

#' Moments for Phase-Type (PH) Distribution
#' 
#' Moments for PH distribution.
#' 
#' @name ph.moment
#' @aliases ph.moment-method
#' @aliases ph.mean ph.var
#' @aliases ph.moment,ANY,ph-method
#' @aliases ph.moment,ANY,herlang-method
#' 
#' @param ph An object of S4 class of PH (\code{\linkS4class{ph}}) or Hyper-Erlang (\code{\linkS4class{herlang}}).
#' @param k An integer of dgrees of moments.
#' @param ... Further arguments for methods.
#' @return
#' \code{ph.mean} and \code{ph.var} give mean and variance of PH.
#' \code{ph.moment} gives a vector of up to k moments.
#' 
#' @details 
#' The PH distribution with parameters \eqn{alpha}, \eqn{Q} and \eqn{xi}:
#' k-th moment; \deqn{k! \alpha (-Q)^{-k} 1}
#' 
#' @note 
#' \code{ph.moment} is a generic function for \code{\linkS4class{ph}} and \code{\linkS4class{herlang}}.
#' 
#' @seealso \code{\link{ph}}, \code{\link{cf1}}, \code{\link{herlang}}
#' 
#' @examples 
#' ## create a PH with specific parameters
#' (param1 <- ph(alpha=c(1,0,0), 
#'               Q=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-1)), 
#'               xi=c(2,2,0)))
#'
#' ## create a CF1 with specific parameters
#' (param2 <- cf1(alpha=c(1,0,0), rate=c(1.0,2.0,3.0)))
#' 
#' ## create a hyper Erlang with specific parameters
#' (param3 <- herlang(shape=c(2,3), mixrate=c(0.3,0.7), rate=c(1.0,10.0)))
#' 
#' ## mean
#' ph.mean(param1)
#' ph.mean(param2)
#' ph.mean(param3)
#' 
#' ## variance
#' ph.var(param1)
#' ph.var(param2)
#' ph.var(param3)
#' 
#' ## up to 5 moments 
#' ph.moment(5, param1)
#' ph.moment(5, param2)
#' ph.moment(5, param3)
#' 
#' @docType methods
#' @keywords distribution
#' @export

setMethod("ph.moment", signature(ph = "ph"),
  function(k, ph, ...) {
  tmp <- ph@alpha
  tmp2 <- 1.0
  res <- numeric(0)
  for (i in 1:k) {
    tmp <- msolve(alpha=1.0, A=-as.matrix(ph@Q), x=tmp, transpose=TRUE)
    tmp2 <- tmp2 * i
    res <- c(res, tmp2 * sum(tmp))
  }
  res
}
)

#' @rdname ph.moment
#' @aliases ph.mean
#' @export

ph.mean <- function(ph) ph.moment(1, ph)

#' @rdname ph.moment
#' @aliases ph.var
#' @export

ph.var <- function(ph) {
  res <- ph.moment(2, ph)
  res[2] - res[1]^2
}

#' @rdname ph
#' @aliases dph
#' @export

dph <- function(x, ph = ph(1), log = FALSE) {
  inv <- order(order(x))
  x <- c(0,sort(x))
  res <- mexp(t=x, transpose=TRUE, x=ph@alpha, A=ph@Q)
  y <- as.vector(ph@xi %*% res$result[,2:length(x)])[inv]
  if (log) {
    log(y)
  } else {
    y
  }
}

#' @rdname ph
#' @aliases pph
#' @export

pph <- function(q, ph = ph(1), lower.tail = TRUE, log.p = FALSE) {
  inv <- order(order(q))
  x <- c(0,sort(q))
  if (lower.tail) {
    al <- c(ph@alpha,0)
    A <- diag.padding.zero(rbind2(cbind2(ph@Q,ph@xi),rep(0,ph@size+1)))
    res <- mexp(t=x, transpose=TRUE, x=al, A=A)
    y <- res$result[ph@size+1,2:length(x)][inv]
    if (log.p) {
      log(y)
    } else {
      y
    }
  } else {
    res <- mexp(t=x, transpose=TRUE, x=ph@alpha, A=ph@Q)
    y <- as.vector(rep(1,ph@size) %*% res$result[,2:length(x)])[inv]
    if (log.p) {
      log(y)
    } else {
      y
    }
  }
}

#' @rdname ph
#' @aliases rph
#' @export

rph <- function(n, ph = ph(1)) {
  if (is(ph, "cf1")) {
    cf1.sample(n, ph)
  } else {
    sapply(1:n, function(a) ph.sample(ph))
  }
}

ph.sample <- function(ph) {
  s <- which(as.vector(rmultinom(n=1, size=1, prob=c(ph@alpha,0)))==1)
  t <- 0
  while (s != ph@size+1) {
    x <- c(ph@Q[s,], ph@xi[s])
    r <- -x[s]
    p <- x / r
    p[s] <- p[s] + 1
    t <- t + rexp(n=1, rate=r)
    s <- which(as.vector(rmultinom(n=1, size=1, prob=p))==1)
  }
  t
}

## S4 methods

setMethod("emfit.df", signature(model = "ph"),
  function(model, ...) {
    model@df
  }
)

setMethod("emfit.print", signature(model = "ph"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Initial : ", model@alpha, "\n")
    cat("Exit    : ", model@xi, "\n")
    cat("Infinitesimal generator : \n")
    print(model@Q)
  }
)

setMethod("emfit.init", signature(model = "ph"),
  function(model, data, verbose = list(), ...) {
    ph.param.random(size=model@size, data=data,
      skelpi=model@alpha, skelQ=as.matrix(model@Q),
      skelxi=model@xi, verbose=verbose, class=class(model@Q))
  }
)

ph.param.random <- function(size, data, skelpi, skelQ, skelxi, verbose, class) {
  if (missing(size)) size <- length(skelpi)
  if (missing(skelpi)) skelpi <- rep(1, size)
  if (missing(skelQ)) skelQ <- matrix(1, size, size)
  if (missing(skelxi)) skelxi <- rep(1, size)

  mean <- mapfit.mean(data)

  alpha <- skelpi * runif(size)
  alpha <- alpha / sum(alpha)

  diag(skelQ) <- 0
  Q <- skelQ * matrix(runif(size*size), size, size)
  xi <- skelxi * runif(size)
  diag(Q) <- -(apply(Q, 1, sum) + xi)

  p <- ph(alpha=alpha, Q=Q, xi=xi, class=class)
  m <- ph.mean(p) / mean
  p@Q <- as(as.matrix(p@Q * m), class)
  p@xi <- p@xi * m
  p
}

#### estep

setMethod("emfit.estep", signature(model = "ph", data = "phdata.wtime"),
  function(model, data, ufact = 1.01, eps = sqrt(.Machine$double.eps), ...) {
    res <- .Call('phfit_estep_gen_wtime', PACKAGE='mapfit', model, data, eps, ufact)
    list(eres=list(etotal=res[[1]], eb=res[[2]], ey=res[[3]], ez=res[[4]], en=res[[5]]), llf=res[[6]])
  })

setMethod("emfit.estep", signature(model = "ph", data = "phdata.group"),
  function(model, data, ufact = 1.01, eps = sqrt(.Machine$double.eps), ...) {

  data@data$instant[is.na(data@data$counts)] <- 0
  data@data$counts[is.na(data@data$counts)] <- -1
  l <- data@size
  if (is.infinite(data@data$time[l])) {
    gdatlast <- data@data$counts[l]
    data@data <- data@data[-l,]
    data@size <- data@size - 1
  } else {
    gdatlast <- 0
  }

  ba <- msolve(alpha=1.0, A=-as.matrix(model@Q), x=model@alpha, transpose=TRUE)

  res <- .Call('phfit_estep_gen_group', PACKAGE='mapfit', model, ba, data, gdatlast, eps, ufact)
  list(eres=list(etotal=res[[1]], eb=res[[2]], ey=res[[3]], ez=res[[4]], en=res[[5]]), llf=res[[6]])
  })

#### mstep

setMethod("emfit.mstep", signature(model = "ph"),
  function(model, eres, data, ...) {
    res <- .Call('phfit_mstep_gen', PACKAGE='mapfit', model, eres, data)
    model@alpha <- res[[1]]
    model@xi <- res[[2]]
    model@Q@x <- res[[3]]
    model
  })

