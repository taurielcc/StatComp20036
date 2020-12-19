#' @title Non-parametric Bootstrapping
#' @description Estimate the standard error and the bias of an estimator using R
#' @param data he data as a vector, matrix or data frame
#' @param func function to be bootstrapped
#' @param B the number of replicates
#' @return a list of standard error and bias
#' @examples
#' \dontrun{
#' data <- 20 * rbeta(1000,2,3)
#' boot(data = data, func = mean, B = 2000)
#' }
#' @importFrom stats sd
#' @export
boot <- function(data,func=NULL, B){
  theta.hat <- func(data)
  #set up the bootstrap
  n <- length(data)      #sample size
  theta.b <- numeric(B)     #storage for replicates
  for (b in 1:B) {
    #randomly select the indices
    i <- sample(1:n, size = n, replace = TRUE)
    dat <- data[i]       #i is a vector of indices
    theta.b[b] <- func(dat)
  }
  #bootstrap estimate of standard error of R
  bias.theta <- mean(theta.b - theta.hat)
  se <- sd(theta.b)
  return(list(bias.b = bias.theta,se.b = se))
}

#' @title Non-parametric Bootstrapping
#' @description Compute the jackknife estimate of standard error using R
#' @param data he data as a vector
#' @param func function to be bootstrapped
#' @return the standard error to be estimated
#' @examples
#' \dontrun{
#' data <- 20 * rbeta(1000,2,3)
#' jack(data = data, func = mean)
#' }
#' @export
jack <- function(data,func=NULL){
  theta.hat <- func(data)
  #set up the bootstrap
  #B is the number of replicates
  n <- length(data)      #sample size
  M <- numeric(n)
  for (i in 1:n) { #leave one out
    y <- data[-i]
    M[i] <- func(y)
  }
  Mbar <- mean(M)
  se.jack <- sqrt(((n - 1)/n) * sum((M - Mbar)^2))
  return(se.jack)
}

#' @title Non-parametric Bootstrapping
#' @description Computes an estimate for each leave-one-out sample using R
#' @param data the data as a vector
#' @param func function to be bootstrapped
#' @param B the number of replicates
#' @return the standard error of both bootstrap method and jackknife-after-bootstrap method
#' @examples
#' \dontrun{
#' data <- 20 * rbeta(1000,2,3)
#' jackafterboot(data = data, func=mean, B = 2000)
#' }
#' @importFrom stats sd
#' @export
jackafterboot <- function(data, func=NULL, B){
  n <- length(data)
  theta.b <- numeric(B)
  # set up storage for the sampled indices
  indices <- matrix(0, nrow = B, ncol = n)
  # jackknife-after-bootstrap step 1: run the bootstrap
  for (b in 1:B) {
    i <- sample(1:n, size = n, replace = TRUE)
    y <- data[i]
    theta.b[b] <- func(y)
    #save the indices for the jackknife
    indices[b, ] <- i
  }
  #jackknife-after-bootstrap to est. se(se)
  se.jack <- numeric(n)
  for (i in 1:n) {
    #in i-th replicate omit all samples with x[i]
    keep <- (1:B)[apply(indices, MARGIN = 1,
                        FUN = function(k) {!any(k == i)})]
    se.jack[i] <- sd(theta.b[keep])
  }
  se.boot <- sd(theta.b)
  se.jackafterboot <- sqrt((n-1) * mean((se.jack - mean(se.jack))^2))
  return(list(se.boot = se.boot, se.jackafterboot=se.jackafterboot))
}

