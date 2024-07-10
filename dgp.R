require(stabledist)
require(sn)

#' @description Generate tensor-valued time series with factor structures
#' See Section 5 for more information
#' @param n sample size
#' @param case if \code{case = 1}, \code{p = (10, 10, 10)}; 
#' else if \code{case = 2}, \code{p = (100, 10, 10)};
#' else if \code{case = 3}, \code{p = (20, 30, 40)}. In all cases, \code{r = (3, 3, 3)}
#' @param dist 'gauss' or 't'
#' @param perc percentage of outliers (\code{0 <= perc <= 100})
#' @return a list containing:
#' \item{x}{ \code{p x n}-array containing the simulated data }
#' \item{r}{ factor number }
#' \item{chi}{ \code{p x n}-array containing the common component }
#' \item{loading}{ List containing the loading matrices }
#' @references Matteo Barigozzi, Haeran Cho and Hyeyoung Maeng (2024) 
#' Tail-robust factor modelling of vector and tensor time series in high dimensions
#' @export
tensor_dgp <- function(n, case = 1, dist = c('gauss', 't'), perc = 0){
  
  K <- 3
  phi <- psi <- .3
  nu <- 3
  burnin <- n
  r <- rep(3, 3)
  
  stopifnot(case %in% c(1, 2, 3) && length(case) == 1)  
  dist <- match.arg(dist)
  stopifnot(perc >= 0 && perc <= 100)
  
  if(case == 1){
    p <- 10 * c(1, 1, 1)
  } else if(case == 2){
    p <- 10 * c(10, 1, 1)
  } else if(case == 3){
    p <- 10 * c(2, 3, 4)
  }
  
  sigma.list <- list()
  for(kk in 1:3){
    sigma <- matrix(1/p[kk], p[kk], p[kk])
    diag(sigma) <- 1
    sigma.list <- c(sigma.list, list(t(chol(sigma))))
  }
  
  if(dist == 'gauss'){
    xi <- array(rnorm((n + burnin) * prod(p)), dim = c(p, n + burnin))
    vecf <- matrix(rnorm((n + burnin) * prod(r)), nrow = prod(r))
  } else if(dist == 't'){
    xi <- array(rt((n + burnin) * prod(p), nu) * sqrt((nu - 2) / nu), dim = c(p, n + burnin))
    vecf <- matrix(rt((n + burnin) * prod(r), nu) * sqrt((nu - 2) / nu), nrow = prod(r))
  }
  
  xi <- as.array(rTensor::ttl(rTensor::as.tensor(xi), list_mat = sigma.list, ms = 1:K)@data)
  vecxi <- t(rTensor::k_unfold(as.tensor(xi), K + 1)@data)
  
  vecxi[, 1] <- sqrt(1 - psi^2) * vecxi[, 1]
  for(tt in 2:dim(vecf)[2]){
    vecf[, tt] <- phi * vecf[, tt - 1] + sqrt(1 - phi^2) * vecf[, tt]
    vecxi[, tt] <- psi * vecxi[, tt - 1] + sqrt(1 - psi^2) * vecxi[, tt]
  }
  vecf <- vecf[, -(1:burnin), drop = FALSE]
  vecxi <- vecxi[, -(1:burnin), drop = FALSE]
  
  LL.list <- list()
  lambda <- 1
  for(kk in K:1){
    LL <- matrix(runif(p[kk] * r[kk], -1, 1), nrow = p[kk])
    LL.list <- c(LL.list, list(LL))
    lambda <- kronecker(lambda, LL)
  }
  LL.list <- rev(LL.list)

  vecchi <- lambda %*% vecf
  chi <- array(vecchi, dim = c(p, n))
  xi <- array(as.vector(vecxi), dim = c(p, n))
  x <- chi + xi
  
  if(perc > 0){
    no <- floor(n * prod(p) * perc / 100)
    so <- sample(prod(dim(x)), no)
    qq <- quantile(abs(x), max(1 - 100/prod(dim(x)), .999))
    x[so] <- sample(c(-1, 1), no, replace = TRUE) * runif(no, qq + 10, qq + 15)
  }
  
  out <- list(x = x, r = r, chi = chi, loading = LL.list)
  
}

#' @description Generate vector-valued time series with factor structures
#' See Appendix B for more information
#' @param n sample size
#' @param p dimensionality
#' @param r factor number
#' @param dep boolean; whether the data are serailly dependent or not;
#' also determines whether the idiosyncratic component is cross-sectionaly dependent
#' @param case taking values from 1 to 5, specifies the model
#' @param perc percentage of outliers (\code{0 <= perc <= 100})
#' @return a list containing:
#' \item{x}{ \code{p x n}-matrix containing the simulated data }
#' \item{r}{ factor number }
#' \item{chi}{ \code{p x n}-matrix containing the common component }
#' \item{loading}{ loading matrix }
#' @references Matteo Barigozzi, Haeran Cho and Hyeyoung Maeng (2024) 
#' Tail-robust factor modelling of vector and tensor time series in high dimensions
#' @export
vector_dgp <- function(n, p, r, dep = FALSE, case = 1, perc = 0){
  
  burnin <- n
  theta <- 1
  stopifnot(perc >= 0 && perc <= 100)
  stopifnot(case %in% c(1, 2, 3, 4, 5) && length(case) == 1)
  
  if(dep){
    rho <- 0.5; beta <- .2; J <- max(10, p / 20)
  } else rho <- beta <- J <- 0
  
  L <- matrix(rnorm(p * r), ncol = r)
  
  if(case == 1){
    f <- matrix(rnorm((n + burnin) * r), nrow = r)
    v <- matrix(rnorm((n + burnin) * p), nrow = p)
  } else if(case == 2){
    f <- matrix(rt((n + burnin) * r, 3), nrow = r) / sqrt(3)
    v <- matrix(rt((n + burnin) * p, 3), nrow = p) / sqrt(3)
  } else if(case == 3){
    f <- matrix(rnorm((n + burnin) * r), nrow = r)
    v <- matrix(rt((n + burnin) * p, 3), nrow = p) / sqrt(3)
  } else if(case == 4){
    f <- matrix(stabledist::rstable((n + burnin) * r, alpha = 1.9, beta = 0, gamma = 1, delta = 0, pm = 0), nrow = r)
    v <- matrix(stabledist::rstable((n + burnin) * p, alpha = 1.9, beta = 0, gamma = 1, delta = 0, pm = 0), nrow = p) 
  } else if(case == 5){
    f <- t(rmst(n = (n + burnin), xi = rep(0, r), Omega = diag(r), alpha = rep(20, r), nu = 3))
    v <- matrix(stabledist::rstable((n + burnin) * p, alpha = 1.9, beta = 0, gamma = 1, delta = 0, pm = 0), nrow = p) 
  }     

  u <- v
  for(t in 2:(n + burnin)){
    for(i in 1:p) u[i, t] <- rho * u[i, t - 1] + (1 - beta) * v[i, t] + beta * sum(v[max(1, i - J):min(i + J, p), t])
  }
  u <- sqrt((1 - rho^2)/(1 + 2 * J * beta^2)) * u
  
  f <- f[, -(1:burnin), drop = FALSE]
  u <- u[, -(1:burnin), drop = FALSE]

  chi <- L %*% f
  x <- chi + sqrt(theta) * u
  
  if(perc > 0){
    no <- floor(n * prod(p) * perc / 100)
    so <- sample(prod(dim(x)), no)
    qq <- quantile(abs(x), max(1 - 100/prod(dim(x)), .999))
    x[so] <- sample(c(-1, 1), no, replace = TRUE) * runif(no, qq + 10, qq + 15)
  }
  
  out <- list(x = x, r = r, chi = chi, loading = L)
  
  out
  
}
