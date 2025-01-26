#' @useDynLib DEdesigns, .registration=TRUE
NULL

#' Differential Evolution
#'
#' DE: A function to construct Designs using various criterions
#' algorithm
#' @param n number of rows(observations)
#' @param m number of columns(variables/factors)
#' @param s number of levels. Default = n
#' @param NP initial population size. Default = 100
#' @param itermax maximum number of iterations. Default = 1000
#' @param pMut the probability of mutation
#' @param pCR Probability of CrossOver
#' @param pGBest Probability of using global best
#' @param replicates number of replications. Determines the output
#' @param seed seed to be set for reproducibility
#' @param ncores number f cores to be used for parallelization/multi-threading.
#' @param method The criterion to be optimized. Only UniPro, MaxPro are implemented
#' @return A list with the phi value, the final design and total computational time taken:
#' \item{bestX}{The final UniPro Design. Only returned when replicates = 1}
#' \item{optValues}{The phi value of the final Design. A vector of length replicates}
#' \item{timeTaken}{Total computational time taken}
#' @examples
#' DE(10, 3) # uses UniPro by default
#' DE(30, 3, replicates=10, method = 'MaxiPro')
#' @usage DE(n, m, s = n, NP = 100, itermax = 1000, pMut =0.2, pCR = 0.3,
#'           pGBest = 0.9, replicates = 1, seed = sample(1e7,1), ncores = NULL,
#'           method = c("UniPro", "MaxPro"))
#'
#' @export
#' @name DE
#'
DE <-function(n, m, s = n, NP = 100, itermax = 1000, pMut = 0.2,
              pCR = 0.3, pGBest = 0.9, replicates = 1,
              seed = sample(1e7,1), cores = NULL,
              method = c("UniPro", "MaxPro", "maximinLHD")){
  methods <- c("UniPro", "MaxPro", "maximinLHD")
  if(is.character(method))
    method <- charmatch(method, methods, 0)
  else
    method <- as.integer(method)
  if(method < 1 || method > 3) stop("invalid method")

  if(is.null(cores)) cores <- as.integer(max(1, detectCores() - 2))

  .Call("DE",
        as.integer(n), as.integer(m), as.integer(s), as.integer(NP),
        as.integer(itermax), as.double(pMut), as.double(pCR),
        as.double(pGBest), as.integer(replicates),
        as.integer(seed), cores, method, PACKAGE = 'DEdesigns')
}

#' Uniform Projection Designs
#'
#' UniPro: A function to construct a Uniform Projection Design using Differential Evolution
#' algorithm
#'
#' @usage UniPro(n, m, s = n, NP = 100, itermax = 1000, pMut =0.2, pCR = 0.3,
#'        pGBest = 0.9, replicates = 1, seed = sample(1e7,1), ncores =NULL)
#'
#' @rdname DE
#' @examples
#' UniPro(10, 3)
#' UniPro(30, 3, replicates=10)
#'
#' @export
#'
#'

UniPro <-function(n, m, s = n, NP = 100, itermax = 1000, pMut = 0.2,
                  pCR = 0.3, pGBest = 0.9, replicates = 1,
                  seed = sample(1e7,1), ncores = NULL){
  mf <- match.call()
  mf$method <- 1L
  mf[[1]] <- quote(DE)
  eval(mf)
}

#' Maximum Projection Designs
#'
#' MaxPro: A function to construct a Maximum Projection Designs using Differential Evolution
#' algorithm
#'
#' @usage MaxPro(n, m, s = n, NP = 100, itermax = 1000, pMut =0.2, pCR = 0.3,
#'        pGBest = 0.9, replicates = 1, seed = sample(1e7,1), ncores =NULL)
#'
#' @rdname DE
#' @examples
#' MaxPro(10, 3)
#' MaxPro(30, 3, replicates=10)
#'
#' @export
#'
#'
MaxPro <-function(n, m, s = n, NP = 100, itermax = 1000, pMut = 0.2,
                  pCR = 0.3, pGBest = 0.9, replicates = 1,
                  seed = sample(1e7,1), ncores = NULL){
  mf <- match.call()
  mf$method <- 2L
  mf[[1]] <- quote(DE)
  eval(mf)
}

#' Maximin-Distance Latin Hypercube Designs
#'
#' maximinLHD: A function to construct  Maximin-Distance Latin Hypercube Designs using Differential Evolution
#' algorithm
#'
#' @usage maximinLHD(n, m, s = n, NP = 100, itermax = 1000, pMut =0.2, pCR = 0.3,
#'        pGBest = 0.9, replicates = 1, seed = sample(1e7,1), ncores =NULL)
#'
#' @rdname DE
#' @examples
#' maximinLHD(10, 3)
#' maximinLHD(30, 3, replicates=10)
#'
#' @export
maximinLHD <-function(n, m, s = n, NP = 100, itermax = 1000, pMut = 0.2,
                  pCR = 0.3, pGBest = 0.9, replicates = 1,
                  seed = sample(1e7,1), ncores = NULL){
  mf <- match.call()
  mf$method <- 3L
  mf[[1]] <- quote(DE)
  eval(mf)
}


printFun <- function(x){
  len <- length(x)
  y <- c(sprintf("%5.3f",head(x, 3)),
         if(len > 10) "...",
         sprintf("%5.3f",tail(tail(x, -3), 3)))|>
    toString()
  cat(" phi Values", if(len>1)c(" [1:",len,"]"),": [", y, "]\n", sep="")
}

#' @export
print.DE <- function(x){
  cat(sprintf("\n Total Time Taken: %8.4f Secs\n", x$timeTaken))
  printFun(x$optValues)
  invisible(x)
}


#' phi
#'
#'
#' A function that computes the UniPro criterion value for a balanced design.
#' Equation 5.1 and 5.2 of the Uniform Projection Designs paper by Sun, Wang and Xu.
#'
#' @usage phi(X, s = nrow(X))
#' @param X Design
#' @param s number of levels for the design X
#' @return the phi value.
#' @examples
#' phi(replicate(3, sample(10)))
#'
#' @export
#'
phi <- function(x, s = nrow(x)){
  .phi(x, s, 0L)
}

#' psi
#'
#'
#' A function that computes the MaxPro criterion value for a given design.
#' Equation 5 of the  Maximum projection designs for computer experiments
#' paper BY V. ROSHAN JOSEPH,EVRENGUL
#'
#' @usage psi(X, s = nrow(X))
#' @param X Design
#' @param s number of levels for the design X
#' @return the phi value.
#' @examples
#' phi(replicate(3, sample(10)))
#'
#' @export
#'
psi <- function(x){
  .phi(x, nrow(x), 1L)
}




.phi <- function(X, s = nrow(X), method = 0L){
  .Call("phi", `mode<-`(X, "integer"), as.integer(s),
        as.matrix(method),  PACKAGE = 'DEdesigns')
}

detectCores <- function(){
  .Call("detectCores", PACKAGE = 'DEdesigns')
}
