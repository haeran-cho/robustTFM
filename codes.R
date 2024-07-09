require(rTensor)
require(tensor)

#' @description Robust tensor factor analysis via truncation
#' @param x input tensor array; the last mode is reserved for the temporal index
#' @param r a vector of integers containing numbers of factors; 
#' if \code{r = NULL}, eigenvalue ratio-based estimators are calculated
#' @param standardize method for centering and standardization; 
#' if \code{standardize = 'median'}, median and MAD are used while 
#' if \code{standardize = 'mean'}, mean and standard deviation are used;
#' if \code{standardize = 'none'}, no centering and standardization is performed
#' @param tau truncation parameter, either a vector or a single value; 
#' if \code{tau = NULL}, a sequence is generated as described in Barigozzi et al (2024)
#' @param kappa truncation parameter, see \code{tau}
#' @param t.ind a vector of integers containing the time indices for which factor and
#' common component are estimated
#' @param nfold number of folds for the cross validation procedure to select tau and kappa
#' @return a list containing:
#' \item{loading}{ List containing the final estimators of loadings }
#' \item{factor}{ \code{(K + 1)}-array containing the factors }
#' \item{chi}{ \code{(K + 1)}-array containing the estimated common component }
#' \item{sc}{ \code{K}-array containing the estimated scale of individual time series }
#' \item{center}{ \code{K}-array containing the estimated center of individual time series }
#' \item{eigval}{ List containing the eiganvalues of the final estimator of mode-wise second moment matrices }
#' \item{init.loading}{ List containing the initial estimators of loadings }
#' \item{init.eigval}{ List containing the eiganvalues of the initial estimator of mode-wise second moment matrices }
#' \item{tau, kappa}{ Choise of truncation parameters }
#' \item{r}{ \code{K}-vector containing the factor numbers }
#' \item{tau.cv}{ Array containing output from the cross validation procedure }
#' @references Matteo Barigozzi, Haeran Cho and Hyeyoung Maeng (2024) 
#' Tail-robust factor modelling of vector and tensor time series in high dimensions
#' @export
rob.tfa <- function(x, r = NULL, standardize = c('median', 'mean', 'none'),
                    tau = NULL, kappa = NULL,
                    t.ind = NULL, nfold = 3){
  
  # main function
  # more data checks
  
  standardize <- match.arg(standardize)
  
  pn <- dim(x)
  K <- length(pn) - 1
  p <- pn[1:K]
  n <- pn[K + 1]
  
  tau.seq <- tau; kappa.seq <- kappa
  
  if(!is.null(r) && length(r) != K) stop('some warning about r')
  if(is.null(t.ind)) t.ind <- 1:n
  
  if(standardize == 'median'){
    center <- apply(x, 1:K, median)
    sc <- apply(x, 1:K, function(z){ ifelse(mad(z) > 0, mad(z), sd(z)) })
    sx <- apply(x, 1:K, function(z){ z <- z - median(z); z / ifelse(mad(z) > 0, mad(z), sd(z)) })
    sx <- aperm(sx, c(2:(K + 1), 1))
  } else if(standardize == 'mean'){
    center <- apply(x, 1:K, mean)
    sc <- apply(x, 1:K, sd)
    sx <- apply(x, 1:K, function(z){ z <- z - mean(z); z / sd(z) })
    sx <- aperm(sx, c(2:(K + 1), 1))
  } else{
    center <- array(0, dim = p)
    sc <- array(1, dim = p)
    sx <- x
  }  
  
  if(is.null(tau.seq)){
    tau.seq <- exp(seq(log(max(abs(sx))), log(median(abs(sx))), length.out = 50))
    ind <- unique(c(1, which(- diff(tau.seq) > diff(range(sx)) * .001), 50))
    tau.seq <- tau.seq[ind]
  }
  
  if(length(tau.seq) == 1){
    tau <- tau.seq; tau.cv.err <- NULL
    if(is.null(r)){
      r.max <- sapply(p, function(z){ min(p - 1, floor(p/2), 20) })
      r <- tnsr.r.est(sx, tau, r.max, c = NULL)
    }
    
    flag <- FALSE
  } else{
    if(!is.null(r)){
      tau.cv.err <- tnsr.cv.tau(sx, r, tau.seq, nfold)
      tau <- tau.seq[which.min(apply(tau.cv.err[, 2,,, drop = FALSE], 1, sum))]
      
      flag <- FALSE
    } else{
      r.max <- sapply(p, function(z){ min(z - 1, floor(z/2), 20) })
      flag <- TRUE; count <- 1
      r.mat <- r <- tnsr.r.est(sx, tau.seq[1], r.max, c = NULL)
      while(flag && count <= 10){
        r0 <- r
        tau.cv.err <- tnsr.cv.tau(sx, r, tau.seq, nfold)
        tau <- tau.seq[which.min(apply(tau.cv.err[, 2,,, drop = FALSE], 1, sum))]
        r <- tnsr.r.est(sx, tau, r.max, c = NULL)
        if(identical(r0, r)) flag <- FALSE
        count <- count + 1
        r.mat <- cbind(r.mat, r)
      }
      if(flag) r <- apply(r.mat, 1, max)
    }
  }
  
  out <- rob.tnsr.loading.est(sx, tau, r)
  eigvec <- out$final$eigvec
  
  if(is.null(kappa.seq) || length(kappa.seq) > 1) kappa <- tau
  
  f <- rob.tnsr.factor.est(sx, kappa, eigvec, r, t.ind)
  
  init.loading <- loading <- list()
  for(kk in 1:K){
    init.loading <- c(init.loading, list(out$first$eigvec[[kk]] * sqrt(p[kk])))
    loading <- c(loading, list(eigvec[[kk]] * sqrt(p[kk])))
  }
  
  chi <- rTensor::ttl(f, list_mat = loading, ms = 1:K)@data
  chi <- apply(chi, K + 1, function(z){ z * sc })
  chi <- array(as.vector(chi), c(p, length(t.ind)))
  
  ans <- list(loading = loading, f = f@data, chi = chi, sc = sc, center = center,
              eigval = out$final$eigval, 
              init.loading = init.loading, init.eigval = out$first$eigval,
              tau = tau, kappa = kappa, r = r, # flag = !flag,
              tau.cv = tau.cv.err)
  
  return(ans)
  
}

#' @description Produces the tensor loading estimator
#' @param 
#' @return a list containing:
#' @noRd
rob.tnsr.loading.est <- function(sx, tau, r){
  
  K <- length(r)
  trunc.x <- trunc.tnsr(sx, tau)
  
  first <- init.loading.est(trunc.x, r)
  if(K > 1) final <- final.loading.est(trunc.x, r, first$eigvec) else final <- first
  
  out <- list(first = first, final = final)
  return(out)
  
}

#' @description Produces the initial estimator
#' @param 
#' @return a list containing:
#' @noRd
init.loading.est <- function(x, r){
  
  K <- length(r)
  p <- dim(x)[1:K]
  n <- dim(x)[K + 1]
  
  eigvec <- eigval <- list()
  for(kk in 1:K){
    if(K == 1){
      Gamma_k <- x %*% t(x) / n
    } else Gamma_k <- tensor::tensor(x, x, (1:(K + 1))[-kk], (1:(K + 1))[-kk]) / (n * prod(p)/p[kk])
    sv <- svd(Gamma_k, nu = r[kk], nv = 0)
    eigval <- c(eigval, list(sv$d))
    eigvec <- c(eigvec, list(sv$u[, 1:r[kk], drop = FALSE]))
  }
  
  out <- list(eigval = eigval, eigvec = eigvec)
  return(out)
  
}

#' @description Produces the final estimator
#' @param 
#' @return a list containing:
#' @noRd
final.loading.est <- function(x, r, init.eigvec){
  
  K <- length(r)
  p <- dim(x)[1:K]
  n <- dim(x)[K + 1]
  
  tx <- rTensor::as.tensor(x)
  eigvec <- eigval <- list()
  for(kk in 1:K){
    dx <- as.array(rTensor::ttl(tnsr = tx, list_mat = lapply(init.eigvec[-kk], t), ms = (1:K)[-kk])@data)
    Gamma_k <- tensor::tensor(dx, dx, (1:(K + 1))[-kk], (1:(K + 1))[-kk]) / (n * prod(p)/p[kk])
    sv <- svd(Gamma_k, nu = r[kk], nv = 0)
    eigval <- c(eigval, list(sv$d))
    eigvec <- c(eigvec, list(sv$u[, 1:r[kk], drop = FALSE]))
  }
  
  out <- list(eigval = eigval, eigvec = eigvec)
  return(out)
  
}

#' @description 
#' @param 
#' @return a list containing:
#' @noRd
rob.tnsr.factor.est <- function(sx, kappa, eigvec, r, tt){
  
  np <- dim(sx)
  K <- length(np) - 1
  p <- np[1:K]
  
  sx <- t(rTensor::k_unfold(rTensor::as.tensor(trunc.tnsr(sx, kappa)), K + 1)@data)[, tt, drop = FALSE]
  sx <- array(as.vector(sx), dim = c(p, length(tt)))
  f <- rTensor::as.tensor(trunc.tnsr(sx, kappa))
  for(kk in 1:K) f <- rTensor::ttm(f, t(eigvec[[kk]]), kk)
  # f <- f@data
  f <- f / sqrt(prod(p))
  
  # vectx <- t(rTensor::k_unfold(rTensor::as.tensor(trunc.tnsr(sx, kappa)), K + 1)@data)[, tt, drop = FALSE]
  # D <- 1
  # for(kk in K:1) D <- kronecker(D, eigvec[[kk]])
  # f <- t(D) %*% vectx / sqrt(prod(p))
  
  return(f)
  
}

#' @description Produces the final estimator
#' @noRd
trunc.tnsr <- function(sx, tau){
  
  sx <- sx * (abs(sx) <= tau) + sign(sx) * tau * (abs(sx) > tau)
  sx
  
}

#' @description cross validation for tau
#' @param sx input data array
#' @param r factor numbers, known
#' @param tau.seq a sequence of truncation parameters
#' @param nfold number of folds
#' @return err
#' @noRd
tnsr.cv.tau <- function(sx, r, tau.seq, nfold = 3){
  
  np <- dim(sx)
  K <- length(np) - 1
  n <- np[K + 1]
  p <- np[1:K]
  vecsx <- t(rTensor::k_unfold(rTensor::as.tensor(sx), K + 1)@data)
  
  cv.err <- array(0, dim = c(length(tau.seq), 2, K, nfold))
  dimnames(cv.err)[[1]] <- tau.seq
  dimnames(cv.err)[[2]] <- c('initial', 'final')
  dimnames(cv.err)[[3]] <- 1:K
  dimnames(cv.err)[[4]] <- 1:nfold
  
  for(mm in 1:nfold){
    ind0 <- (floor(n/nfold) * (mm - 1) + 1):(floor(n/nfold) * mm)
    ind1 <- setdiff(1:n, ind0)
    
    sx0 <- array(as.vector(vecsx[, ind0, drop = FALSE]), dim = c(p, length(ind0)))
    u0 <- rob.tnsr.loading.est(sx0, tau = max(tau.seq) * 1.1, r)$final$eigvec
    ls0.list <- c()
    for(kk in 1:K) ls0.list <- c(ls0.list, list(u0[[kk]] %*% t(u0[[kk]])))
    
    sx1 <- array(as.vector(vecsx[, ind1, drop = FALSE]), dim = c(p, length(ind1)))
    
    for(ii in 1:length(tau.seq)){
      out <- rob.tnsr.loading.est(sx1, tau = tau.seq[ii], r)
      for(kk in 1:K){
        cv.err[ii, 1, kk, mm] <- 1 - sum(diag(out$first$eigvec[[kk]] %*% t(out$first$eigvec[[kk]]) %*% ls0.list[[kk]])) / r[kk]
        cv.err[ii, 2, kk, mm] <- 1 - sum(diag(out$final$eigvec[[kk]] %*% t(out$final$eigvec[[kk]]) %*% ls0.list[[kk]])) / r[kk]
      }
    }
  }
  cv.err
  
}

#' @description 
#' @param tx 
#' @param r.max number of folds
#' @return r.est
#' @noRd
tnsr.r.est <- function(sx, tau, r.max = NULL, c = NULL){
  
  np <- dim(sx)
  K <- length(np) - 1
  n <- np[K + 1]
  p <- np[1:K]
  if(is.null(r.max)) r.max <- sapply(p, function(z){ min(floor(z/2), z - 1, 20) })
  
  tx <- trunc.tnsr(sx, tau)
  
  init.eigvec <- init.eigval <- list()
  for(kk in 1:K){
    if(K == 1){
      Gamma_k <- tx %*% t(tx) / n
    } else Gamma_k <- tensor::tensor(tx, tx, (1:(K + 1))[-kk], (1:(K + 1))[-kk]) / (n * prod(p)/p[kk])
    sv <- svd(Gamma_k, nu = r.max[kk], nv = 0)
    init.eigval <- c(init.eigval, list(sv$d))
    init.eigvec <- c(init.eigvec, list(sv$u))
  }
  
  if(FALSE){
    par(mfcol = c(2, K), mar = c(2, 2, 1, 1))
    for(kk in 1:K){
      if(is.null(c)) c0 <- 1 / init.eigval[[kk]][1] else c0 <- c
      plot( eval <- init.eigval[[kk]][1:r.max[kk]] / (init.eigval[[kk]][1:r.max[kk] + 1] + c0 ))
      print(which.max(eval))
      plot( eval <- init.eigval[[kk]][1:r.max[kk]] / (init.eigval[[kk]][1:r.max[kk] + 1] ))
      print(which.max(eval))
    }
  }
  
  if(K >= 2){
    tx <- rTensor::as.tensor(tx)
    eigval <- list()
    for(kk in 1:K){
      dx <- as.array(rTensor::ttl(tnsr = tx, list_mat = lapply(init.eigvec[-kk], t), ms = (1:K)[-kk])@data)
      Gamma_k <- tensor::tensor(dx, dx, (1:(K + 1))[-kk], (1:(K + 1))[-kk]) / (n * prod(p)/p[kk])
      sv <- svd(Gamma_k, nu = 0, nv = 0)
      eigval <- c(eigval, list(sv$d))
    }
  } else eigval <- init.eigval
  
  r.est <- c()
  for(kk in 1:K){
    if(is.null(c)) c0 <- 1 / eigval[[kk]][1] else c0 <- c
    eval <- eigval[[kk]][1:r.max[kk]] / (eigval[[kk]][1:r.max[kk] + 1] + c0)
    r.est <- c(r.est, which.max(eval))
  }
  
  if(FALSE){
    par(mfcol = c(2, K), mar = c(2, 2, 1, 1))
    for(kk in 1:K){
      if(is.null(c)) c0 <- 1 / eigval[[kk]][1] else c0 <- c
      plot( eval <- eigval[[kk]][1:r.max[kk]] / (eigval[[kk]][1:r.max[kk] + 1] + c0 ))
      print(which.max(eval))
      plot( eval <- eigval[[kk]][1:r.max[kk]] / (eigval[[kk]][1:r.max[kk] + 1] + 0 ))
      print(which.max(eval))
    }
  }  
  
  r.est
  
}
