#' Classes of MAP
#' 
#' Parameters for MAP and MMPP.
#' 
#' @docType class
#' @name map-class
#' @aliases gmmpp-class
#' @slot size The number of phases (internal states).
#' @slot alpha A probability (row) vector to decide an initial phase.
#' @slot D0 A square matrix that means transition rates without arrivals.
#' @slot D1 A square matrix that means transition rates with arrivals. In the case of MMPP, D1 should be a diagonal matrix.
#' @slot df The number of free parameters.
#' 
#' @note 
#' Objects are usually created by \link{map}, \link{mmpp} or \link{gmmpp}.
#' 
#' @seealso Classes \code{\linkS4class{erhmm}}.
#' 
#' @examples 
#' ## create an MAP (full matrix) with 5 phases
#' map(5)
#' 
#' ## create an MAP (full matrix) with 5 phases
#' map(size=5)
#' 
#' ## create an MMPP with 5 states
#' mmpp(5)
#' 
#' ## create an MMPP with 5 states for approximate
#' ## estimation
#' gmmpp(5)
#' 
#' ## create an MAP with specific parameters
#' (param <- map(alpha=c(1,0,0),
#'               D0=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-4)),
#'               D1=rbind(c(1,1,0),c(1,0,1),c(2,0,1))))
#'
#' ## marginal moments of MAP
#' map.mmoment(k=3, map=param)
#' 
#' ## joint moments of MAP
#' map.jmoment(lag=1, map=param)
#' 
#' ## k-lag correlation
#' map.acf(map=param)
#' 
#' @keywords classes
#' @exportClass map
NULL

#' Markovian Arrival Process (MAP)
#' 
#' Functions to generate an object of \code{\linkS4class{map}}.
#' 
#' @aliases mmpp gmmpp
#' 
#' @param size An integer for the number of phases.
#' @param alpha A vector of probabilities for determing an initial phase.
#' @param D0 An object of Matrix class for the initesmal generator without arrivals.
#' @param D1 An object of Matrix class for the initesmal generator with arrivals.
#' @param class Name of Matrix class for \code{D0} and \code{D1}.
#' @return 
#' \code{map} gives an object of general MAP.
#' \code{mmpp} gives an object of MMPP with default parameters.
#' \code{gmmpp} gives an object of MMPP which uses an approximate estimation algorithm.
#' 
#' @details 
#' MAP parameters are \eqn{alpha}, \eqn{D_0} and \eqn{D_1}. \eqn{alpha} is the
#' probability vector to determine an initial phase at time 0. \eqn{D_0} is an
#' infinitesimal generator of underlyinc continuous-time Markov chain (CTMC)
#' without arrival. \eqn{D_1} is an infinitesimal generator of CTMC with arrival.
#' The infinitesimal generator of underlying CTMC becomes \eqn{D_0+D_1}.
#' In the stationary case, \eqn{\alpha} is often given by a stationary vector
#' satisfying \eqn{\alpha (D_0+D_1) = \alpha}.
#' 
#' \code{mmpp} generates an object of a specific MAP called MMPP.
#' MMPP (Markov modulated Poisson process) is an MAP whose \eqn{D_1} is given by
#' a diagonal matrix. Unlike to general MAPs, MMPP never changes the phase at
#' which an arrival occurs.
#' 
#' \code{gmmpp} generates an object of \code{\linkS4class{gmmpp}}, which is
#' exactly same as MMPP. In the estimation algorithm, \code{\linkS4class{gmmpp}}
#' class uses an approximate method.
#' 
#' @note 
#' \code{map} and \code{gmmpp} require either \code{size} or (\code{alpha}, \code{D0}, \code{D1}).
#' 
#' @seealso \code{\link{erhmm}}, \code{\link{map.mmoment}}, \code{\link{map.jmoment}}, \code{\link{map.acf}}
#' 
#' @examples
#' ## create an MAP (full matrix) with 5 phases
#' map(5)
#' 
#' ## create an MAP (full matrix) with 5 phases
#' map(size=5)
#' 
#' ## create an MMPP with 5 states
#' mmpp(5)
#' 
#' ## create an MMPP with 5 states for approximate
#' ## estimation
#' gmmpp(5)
#' 
#' ## create an MAP with specific parameters
#' (param <- map(alpha=c(1,0,0),
#'               D0=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-4)),
#'               D1=rbind(c(1,1,0),c(1,0,1),c(2,0,1))))
#'
#' ## marginal moments of MAP
#' map.mmoment(k=3, map=param)
#' 
#' ## joint moments of MAP
#' map.jmoment(lag=1, map=param)
#' 
#' ## k-lag correlation
#' map.acf(map=param)
#' 
#' @export
  
map <- function(size, alpha, D0, D1, class="CsparseMatrix") {
  if (missing(size)) {
    if (missing(alpha) || missing(D0) || missing(D1)) {
      stop("alpha, D0 and D1 are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(D0) && !missing(D1)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (missing(alpha)) {
        # warning("alpha is set by default")
        alpha <- rep(1.0/size, size)
      }
      if (missing(D1)) {
        # warning("D1 is set by default")
        D1 <- matrix(1.0, size, size)
      }
      if (missing(D0)) {
        # warning("D0 is set by default")
        D0 <- matrix(1.0, size, size)
        diag(D0) <- -(rep(size-1, size) + apply(D1,1,sum))
      }
    }
  }
  if (!is(D0, class)) {
    D0 <- as(D0, class)
  }
  if (!is(D1, class)) {
    D1 <- as(D1, class)
  }
  zero <- 1.0e-8
  df <- sum(abs(alpha) > zero) - 1 + sum(abs(D0) > zero) +
    sum(abs(D1) > zero) - size + sum(abs(diag(D0)) < zero)
  new("map", size=size, alpha=alpha, D0=D0, D1=D1, df=df)
}

#' @rdname map
#' @aliases mmpp
#' @export

mmpp <- function(size, class="CsparseMatrix") {
  alpha <- rep(1.0/size, size)
  D1 <- diag(1.0, size)
  D0 <- matrix(1.0, size, size)
  diag(D0) <- -(rep(size-1, size) + apply(D1,1,sum))
  map(alpha=alpha, D0=D0, D1=D1, class=class)
}

map.bidiag <- function(size, class="CsparseMatrix") {
  if (size <= 1) {
    map(size, class=class)
  } else {
    alpha <- rep(1.0/size, size)
    D0 <- matrix(0, size, size)
    D1 <- matrix(1, size, size)
    for (i in 1:(size-1)) {
      D0[i,i] <- -1-size
      D0[i,i+1] <- 1
    }
    D0[size,size] <- -size
    map(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

map.tridiag <- function(size, class="CsparseMatrix") {
  if (size <= 2) {
    map(size, class=class)
  } else {
    alpha <- rep(1.0/size, size)
    D0 <- matrix(0, size, size)
    D1 <- matrix(1, size, size)
    D0[1,1] <- -1-size
    D0[1,2] <- 1
    for (i in 2:(size-1)) {
      D0[i,i] <- -2-size
      D0[i,i-1] <- 1
      D0[i,i+1] <- 1
    }
    D0[size,size] <- -1-size
    D0[size,size-1] <- 1
    map(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

mmpp.tridiag <- function(size, class="CsparseMatrix") {
  if (size <= 2) {
    map(size, class=class)
  } else {
    alpha <- rep(1.0/size, size)
    D0 <- matrix(0, size, size)
    D1 <- diag(1, size)
    D0[1,1] <- -2
    D0[1,2] <- 1
    for (i in 2:(size-1)) {
      D0[i,i] <- -3
      D0[i,i-1] <- 1
      D0[i,i+1] <- 1
    }
    D0[size,size] <- -2
    D0[size,size-1] <- 1
    map(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

#' Moments for Markovian arrival pcess (MAP)
#' 
#' Moments for MAP.
#' 
#' @aliases map.mmoment map.jmoment map.acf
#' 
#' @param map An object of S4 class of MAP (\code{\linkS4class{map}}, \code{\linkS4class{gmmpp}}).
#' @param k An integer of dgrees of moments.
#' @param lag An integer of time lag for corrleation.
#' @return 
#' \code{map.mmoment} gives a vector of up to k moments.
#' \code{map.jmoment} gives a matrix of \eqn{s_{ij}(lag), i=1,..,n, j=1,..,n} where n is the size of phases.
#' \code{map.acf} gives a vector of up to n-lag correlation, where n is the size of phases.
#' 
#' @details 
#' MAP parameters are \eqn{\alpha}, \eqn{D_0} and \eqn{D_1};
#' \deqn{P = (-D_0)^{-1} D_1} and \deqn{s P = s.}
#' 
#' Then the moments for MAP are marginal moment; \deqn{m_k = k! s (-D_0)^{-k} 1,}
#' joint moment; \deqn{s_{ij}(lag) = i! j! s (-D_0)^{-i} P^{lag} (-D_0)^{-j} 1,}
#' k-lag correlation (autocorrelation); \deqn{rho(lag) = (s_{11}(lag) - m_1^2)/(m_2 - m_1^2)}
#' 
#' @note 
#' \code{map.mmoment} is a generic function for \code{\linkS4class{ph}} and \code{\linkS4class{herlang}}.
#' 
#' @seealso \code{\link{map}}, \code{\link{gmmpp}}, \code{\link{erhmm}}
#' 
#' @examples 
#' ## create an MAP with specific parameters
#' (param1 <- map(alpha=c(1,0,0),
#'                D0=rbind(c(-4,2,0),c(2,-5,1),c(1,0,-4)),
#'                D1=rbind(c(1,1,0),c(1,0,1),c(2,0,1))))
#'
#' ## create an ER-HMM with specific parameters
#' (param2 <- erhmm(shape=c(2,3), alpha=c(0.3,0.7),
#'                  rate=c(1.0,10.0),
#'                  P=rbind(c(0.3, 0.7), c(0.1, 0.9))))
#'
#' ## marginal moments of MAP
#' map.mmoment(k=3, map=param1)
#' map.mmoment(k=3, map=as(param2, "map"))
#' 
#' ## joint moments of MAP
#' map.jmoment(lag=1, map=param1)
#' map.jmoment(lag=1, map=as(param2, "map"))
#' 
#' ## k-lag correlation
#' map.acf(map=param1)
#' map.acf(map=as(param2, "map"))
#' 
#' @keywords distribution
#' @export

map.mmoment <- function(k, map) {
  D0 <- as.matrix(map@D0)
  D1 <- as.matrix(map@D1)
  piv <- as.vector(ctmc.st(D0+D1) %*% D1)
  piv <- piv / sum(piv)
  tmp <- piv
  tmp2 <- 1.0
  res <- numeric(0)
  for (i in 1:k) {
    tmp <- msolve(alpha=1.0, A=-D0, x=tmp, transpose=TRUE)
    tmp2 <- tmp2 * i
    res <- c(res, tmp2 * sum(tmp))
  }
  res
}

#' @rdname map.mmoment
#' @aliases map.jmoment
#' @export

map.jmoment <- function(lag, map) {
  D0 <- as.matrix(map@D0)
  D1 <- as.matrix(map@D1)
  piv <- as.vector(ctmc.st(D0+D1) %*% D1)
  piv <- piv / sum(piv)
  P <- mpow(mpow(A=-D0, m=-1) %*% D1, m=lag)
  vone <- rep(1, map@size)

  fmat <- matrix(0, map@size, map@size)
  fmat[,1] <- piv
  for (i in 2:map@size) {
    fmat[,i] <- msolve(alpha=1.0, A=-D0, x=fmat[,i-1], transpose=TRUE)
  }

  res <- matrix(0, map@size, map@size)
  tmp <- vone
  for (j in 1:map@size) {
    for (i in 1:map@size) {
      res[i,j] <- gamma(i) * gamma(j) * as.vector(fmat[,i] %*% P %*% tmp)
    }
    tmp <- msolve(alpha=1.0, A=-D0, x=tmp)
  }
  res
}

#' @rdname map.mmoment
#' @aliases map.acf
#' @export

map.acf <- function(map) {
  D0 <- as.matrix(map@D0)
  D1 <- as.matrix(map@D1)
  piv <- as.vector(ctmc.st(D0+D1) %*% D1)
  piv <- piv / sum(piv)
  P <- mpow(A=-D0, m=-1) %*% D1
  vone <- rep(1, map@size)

  piv <- msolve(alpha=1.0, A=-D0, x=piv, transpose=TRUE)
  vone <- msolve(alpha=1.0, A=-D0, x=vone)

  mres <- map.mmoment(2, map)
  (sapply(1:map@size, function(k) as.vector(piv %*% mpow(P,k) %*% vone)) - mres[1]^2) / (mres[2] - mres[1]^2)
}

## S4 methods

setMethod("emfit.print", signature(model = "map"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Initial : ", model@alpha, "\n")
    cat("Infinitesimal generator : \n")
    print(model@D0)
    print(model@D1)
  }
)

setMethod("emfit.df", signature(model = "map"),
  function(model, stationary = FALSE, ...) {
    if (stationary)
      model@df - model@size + 1
    else
      model@df
  }
)

## init

setMethod("emfit.init", signature(model = "map", data = "mapdata"),
  function(model, data, verbose = list(), ...) {
    map.param.kmeans(size=model@size, data=data@data, 
      skelal=model@alpha, skelD0=as.matrix(model@D0),
      skelD1=as.matrix(model@D1), verbose=verbose, class=class(model@D0), ...)
  }
)

map.param.kmeans <- function(size, data, skelal, skelD0, skelD1, verbose, class, ...) {
  if (missing(skelal)) skelal <- rep(1, size)
  if (missing(skelD0)) skelD0 <- matrix(1, size, size)
  if (missing(skelD1)) skelD1 <- matrix(1, size, size)

  maxt <- max(data$time)
  maxg <- max(data$counts + data$instant)
  x <- cbind(data$time, (data$counts + data$instant))
  v <- cbind(data$time / maxt, (data$counts + data$instant) / maxg)
  result <- kmeans(v, size)
  
  diagelem <- sapply(1:size, function(k) {
    (sum(x[result$cluster == k,2]) + 1) /
    sum(x[result$cluster == k,1])})

  alpha <- skelal * runif(size)
  alpha <- alpha / sum(alpha)

  diag(skelD0) <- 0
  D0 <- skelD0 * matrix(runif(size*size), size, size)
  D1 <- skelD1 * matrix(runif(size*size), size, size)
  d <- diagelem / (apply(D0, 1, sum) + apply(D1, 1, sum))
  D0 <- D0 * d
  D1 <- D1 * d
  diag(D0) <- -diagelem
  map(alpha=alpha, D0=D0, D1=D1, class=class)
}

#### estep

setMethod("emfit.estep", signature(model = "map", data = "mapdata"),
  function(model, data, ufact = 1.01, eps = 1.0e-8, ...) {
    data@data$instant[is.na(data@data$counts)] <- 0
    data@data$counts[is.na(data@data$counts)] <- -1
    res <- .Call('mapfit_estep_gen_group', PACKAGE='mapfit', model, data, eps, ufact)
    list(eres=list(eb=res[[1]], ez=res[[2]], en0=res[[3]], en1=res[[4]]), llf=res[[5]])
  })

#### mstep

setMethod("emfit.mstep", signature(model = "map"),
  function(model, eres, data, stationary = TRUE, ...) {
    res <- .Call('mapfit_mstep_gen', PACKAGE='mapfit', model, eres, data)
    model@D0@x <- res[[2]]
    model@D1@x <- res[[3]]
    if (stationary)
      model@alpha <- ctmc.st(as.matrix(model@D0 + model@D1))
    else
      model@alpha <- res[[1]]
    model
  })

