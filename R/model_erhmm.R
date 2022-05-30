#' Class of ER-HMM
#' 
#' Parameters for an ER-HMM (Hidden Markov Model with Erlang outputs).
#' 
#' @docType class
#' @name erhmm-class
#' @slot size The number of HMM states.
#' @slot alpha A vector of initial probabilities for HMM states.
#' @slot shape Shape parameters for Erlang distributions. The sum of shape parameters is the number of phases of MAP.
#' @slot rate Rate parameters for Erlang distributions.
#' @slot P An object of Matrix class for a transition probability matrix of HMM.
#' 
#' @note 
#' Objects are usually created by an \link{erhmm}.
#' 
#' This class can be converted to \code{\linkS4class{map}}.
#' 
#' Methods:
#' \describe{
#' \item{ph.moment}{\code{signature(ph = "herlang")}: ... }
#' \item{emfit.init}{\code{signature(model = "herlang", data = "phdata.wtime")}: ... }
#' \item{emfit.init}{\code{signature(model = "herlang", data = "phdata.group")}: ... }
#' \item{emfit.estep}{\code{signature(model = "herlang", data = "phdata.wtime")}: ... }
#' \item{emfit.estep}{\code{signature(model = "herlang", data = "phdata.group")}: ... }
#' \item{emfit.mstep}{\code{signature(model = "herlang")}: ... }
#' }
#' 
#' @seealso Classes \code{\linkS4class{map}} and \code{\linkS4class{gmmpp}}.
#' 
#' @examples 
#' ## create an ER-HMM consisting of two Erlang components with
#' ## shape parameters 2 and 3.
#' erhmm(c(2,3))
#' 
#' ## create an ER-HMM consisting of two Erlang components with
#' ## shape parameters 2 and 3.
#' erhmm(shape=c(2,3))
#' 
#' ## create an ER-HMM with specific parameters
#' (param <- erhmm(shape=c(2,3), alpha=c(0.3,0.7),
#'                 rate=c(1.0,10.0),
#'                 P=rbind(c(0.3, 0.7), c(0.1, 0.9))))
#'
#' ## convert to a general MAP
#' as(param, "map")
#' 
#' ## marginal moments of MAP
#' map.mmoment(k=3, map=as(param, "map"))
#' 
#' ## joint moments of MAP
#' map.jmoment(lag=1, map=as(param, "map"))
#' 
#' ## k-lag correlation
#' map.acf(map=as(param, "map"))
#' 
#' @keywords classes
#' @exportClass erhmm
NULL

#' ER-HMM (HMM with Erlang outputs)
#' 
#' A function to generate an object of \code{\linkS4class{erhmm}}.
#' 
#' @param shape An integer vector of shape parameters of Erlang outputs.
#' @param alpha A vector for initial probabilities of HMM states.
#' @param rate A vector of rate parameters of Erlang outputs.
#' @param P An object of Matrix class for a transition probability matrix of HMM.
#' @param class Name of Matrix class for \code{P}.
#' @return 
#' \code{erhmm} gives an object of ER-HMM.
#' 
#' @details 
#' ER-HMM has parameters \eqn{\alpha}, \eqn{shape}, \eqn{rate} and \eqn{P}.
#' HMM state chages according to a discrete-time Markov chain with transition
#' matrix \eqn{P}. At each HMM state, there is an inherent Erlang distriution
#' as an output. This model can be converted to a MAP.
#' 
#' @note
#' \code{erhmm} requires shape parameters. Other parameters have default values.
#' 
#' @seealso
#' \code{\link{map}}, \code{\link{gmmpp}}, \code{\link{map.mmoment}},
#' \code{\link{map.jmoment}}, \code{\link{map.acf}}
#' 
#' @examples 
#' ## create an ER-HMM consisting of two Erlang components with
#' ## shape parameters 2 and 3.
#' erhmm(c(2,3))
#' 
#' ## create an ER-HMM consisting of two Erlang components with
#' ## shape parameters 2 and 3.
#' erhmm(shape=c(2,3))
#' 
#' ## create an ER-HMM with specific parameters
#' (param <- erhmm(shape=c(2,3), alpha=c(0.3,0.7),
#'                 rate=c(1.0,10.0),
#'                 P=rbind(c(0.3, 0.7), c(0.1, 0.9))))
#'
#' ## convert to a general MAP
#' as(param, "map")
#' 
#' ## marginal moments of MAP
#' map.mmoment(k=3, map=as(param, "map"))
#' 
#' ## joint moments of MAP
#' map.jmoment(lag=1, map=as(param, "map"))
#' 
#' ## k-lag correlation
#' map.acf(map=as(param, "map"))
#' 
#' @export

erhmm <- function(shape, alpha, rate, P, class="CsparseMatrix") {
    size <- length(shape)
    if (missing(alpha)) {
      alpha <- rep(1/size, size)
    }
    if (missing(rate)) {
      rate <- rep(1, size)
    }
    if (missing(P)) {
      P <- matrix(1/size, size, size)
    }
    if (!is(P, class)) {
      P <- as(P, class)
    }
    new("erhmm", size=size, alpha=alpha, shape=shape, rate=rate, P=P)
}

setAs("erhmm", "map", function(from, to) {
  phsize <- sum(from@shape)
  index <- cumsum(from@shape)
  sindex <- c(1, index + 1)[1:from@size]
  eindex <- index
  alpha <- numeric(phsize)
  alpha[sindex] <- from@alpha
##  xi[eindex] <- from@rate

  v <- numeric(0)
  i <- numeric(0)
  j <- numeric(0)
  for (k in 1:from@size) {
    i <- c(i, sindex[k]:eindex[k])
    j <- c(j, sindex[k]:eindex[k])
    v <- c(v, rep(-from@rate[k], from@shape[k]))
  }
  for (k in 1:from@size) {
    if (from@shape[k] != 1) {
      i <- c(i, sindex[k]:(eindex[k]-1))
      j <- c(j, (sindex[k]+1):eindex[k])
      v <- c(v, rep(from@rate[k], from@shape[k]-1))
    }
  }
  D0 <- sparseMatrix(dims=c(phsize,phsize), i=i, j=j, x=v)

  v <- numeric(0)
  i <- numeric(0)
  j <- numeric(0)
  for (k1 in 1:from@size) {
    for (k2 in 1:from@size) {
      i <- c(i, eindex[k1])
      j <- c(j, sindex[k2])
      v <- c(v, from@rate[k1] * from@P[k1,k2])
    }
  }
  D1 <- sparseMatrix(dims=c(phsize,phsize), i=i, j=j, x=v)

  map(alpha=alpha, D0=D0, D1=D1)
  })

## S4 methods

setMethod("emfit.print", signature(model = "erhmm"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Shape   : ", model@shape, "\n")
    cat("Initial : ", model@alpha, "\n")
    cat("Rate    : ", model@rate, "\n")
    cat("Transition probability : \n")
    print(model@P)
  }
)

setMethod("emfit.df", signature(model = "erhmm"),
  function(model, stationary = FALSE, ...) {
    if (stationary) {
      d <- model@size * (model@size - 1) + 2 * model@size - 1
    } else {
      d <- model@size * (model@size - 1) + 3 * model@size - 2
    }
    return(d)
  }
)

# init

setMethod("emfit.init", signature(model = "erhmm", data = "mapdata.time"),
  function(model, data, verbose = list(), ...) {
    erhmm.param.kmeans(shape=model@shape, data=data@data, 
      skelal=model@alpha, skelP=as.matrix(model@P), verbose=verbose, ...)
  }
)

erhmm.param.kmeans <- function(shape, data, skelal, skelP, verbose, ...) {
  if (missing(skelal)) skelal <- rep(1, size)
  if (missing(skelP)) skelP <- matrix(1, size, size)
  size <- length(shape)
  tmp <- numeric(size)
  rate <- numeric(size)

  if (size >= 2) {
    result <- kmeans(data$time, size)
    for (k in 1:size) {
      m <- base::mean(data$time[result$cluster == k])
      s2 <- var(data$time[result$cluster == k])
      tmp[k] <- round(m^2 / s2)
      rate[k] <- 1.0 / m
    }
    rate <- rate[rank(tmp)] * shape
  } else {
    m <- base::mean(data$time)
    rate[1] <- shape[1] / m
  }
  alpha <- skelal * runif(size)
  alpha <- alpha / sum(alpha)
  P <- skelP * matrix(runif(size*size), size, size)
  P <- P / apply(P, 1, sum)
  erhmm(shape=shape, alpha=alpha, rate=rate, P=P, ...)
}

### estep

setMethod("emfit.estep", signature(model = "erhmm", data = "mapdata.time"),
  function(model, data, ...) {
    res <- .Call('mapfit_hmm_erlang_estep', PACKAGE='mapfit', model, data)
    list(eres=list(eb=res[[1]], en=res[[2]], ew0=res[[3]], ew1=res[[4]]), llf=res[[5]])
  })

#### mstep

setMethod("emfit.mstep", signature(model = "erhmm"),
  function(model, eres, data, stationary = TRUE, ...) {
    res <- .Call('mapfit_hmm_erlang_mstep', PACKAGE='mapfit', model, eres, data)
    model@rate <- res[[2]]
    model@P@x <- res[[3]]
    if (stationary) {
      model@alpha <- ctmc.st(as.matrix(model@P))
    } else {
      model@alpha <- res[[1]]      
    }
    model
  })

