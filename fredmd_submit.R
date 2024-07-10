# run these codes prior to analysis

install.packages("readr")
install.packages("pracma")
devtools::install_github("cykbennie/fbi")
library(fbi)

source('codes.R')

fluctuation_test <- function (loss1, loss2, mu = 0.2, dmv_fullsample = TRUE, lag_truncate = 0, 
                              time_labels = NULL, conf_level = 0.05, ...) {
  # modified ver. of murphydiagram::fluctuation_test
  
  if (length(loss1) != length(loss2)) {
    stop("Vectors of losses must have the same length")
  }
  if (all(abs(seq(from = 0.1, to = 0.9, by = 0.1) - mu) > 1e-12)) {
    stop("mu must be in {0.1, 0.2, ..., 0.9}")
  }
  if (!lag_truncate %in% 0:5) {
    stop("lag_truncate must be in {0, 1, .., 5}")
  }
  if (!is.null(time_labels) & length(time_labels) != length(loss1)) {
    warning("Specified time labels are inconsistent - simple integers used instead")
    time_labels <- NULL
  }
  if (!conf_level %in% c(0.05, 0.1)) {
    stop("significance_level must be either 0.05 or 0.1")
  }
  CV <- cbind(seq(from = 0.1, to = 0.9, by = 0.1), c(3.393, 
                                                     3.179, 3.012, 2.89, 2.779, 2.634, 2.56, 2.433, 2.248), 
              c(3.17, 2.948, 2.766, 2.626, 2.5, 2.356, 2.252, 2.13, 
                1.95))
  vHAC2 <- function(ld, lag_truncate) {
    murphydiagram:::vHAC(ld, k = lag_truncate, meth = "Bartlett")$hac
  }
  cex_gen <- 1
  P <- length(loss1)
  m <- round(mu * P)
  ld <- loss1 - loss2
  dm_num <- dm_den <- rep(0, P - m + 1)
  for (jj in m:P) {
    ind <- which(m:P == jj)
    ld_tmp <- ld[(jj - m + 1):jj]
    dm_num[ind] <- mean(ld_tmp)
    dm_den[ind] <- sqrt(vHAC2(ld_tmp, lag_truncate)/m)
  }
  s2hat <- c(vHAC2(ld, lag_truncate))
  dm1 <- sqrt(m) * dm_num/sqrt(s2hat)
  dm2 <- dm_num/dm_den
  if (dmv_fullsample) {
    dm_final <- dm1
  }
  else {
    dm_final <- dm2
  }
  if (conf_level == 0.05) {
    CVs <- CV[abs(CV[, 1] - mu) < 1e-12, 2] * c(-1, 1)
  }
  else if (conf_level == 0.1) {
    CVs <- CV[abs(CV[, 1] - mu) < 1e-12, 3] * c(-1, 1)
  }
  plot(x = m:length(loss1), y = dm_final, ylim = max(abs(dm_final)) * c(-1, 1),
       bty = "n", ylab = "", xlab = "Time (End of Rolling Window)", 
       type = "l", col = "cornflowerblue", lwd = 2.5, axes = FALSE, 
       cex.lab = cex_gen, ...)
  axis(2, cex.axis = cex_gen)
  abline(h = CVs, lwd = 3.5)
  abline(h = 0, lwd = 1.8, lty = 2)
  if (is.null(time_labels)) {
    time_labels <- 1:length(loss1)
  }
  inds <- floor(seq(from = m, to = length(time_labels), length.out = 5))
  axis(1, at = inds, labels = time_labels[inds], cex.axis = cex_gen)
  list(df = data.frame(time = time_labels[m:length(loss1)], 
                       dmstat = dm_final), CV = CVs)
}

# data preparation

data <- fredmd("fredmd.csv")
data <- as.data.frame(data)
data$date <- format(as.Date(data$date, "%Y-%m-%d"), "%Y-%m")
data <- data[data$date > "1959-12" & data$date < "2024-01", ]
data <- data[, !is.na(data[1, ])]
data <- data[, apply(data, 2, function(z){ sum(is.na(z)) }) == 0]

x <- t(as.matrix(data[, - 1]))
x <- x[!(rownames(x) %in% c("NONBORRES", "M1SL")), ]
dim(x) # n = 768, p = 109

## rolling window-based forecasting exercise

library(HDRFA)

r <- c(5, 6)[1]
m <- 12 * 10

p <- dim(x)[1]; n <- dim(x)[2]

err <- array(0, dim = c(2, n, 2, 3, 24))
dimnames(err)[[3]] <- c('in-sample', 'forecast')
dimnames(err)[[4]] <- c('trunc', 'no-trunc', 'HDRFA')
dimnames(err)[[5]] <- 1:24
trunc.param.seq <- rep(0, n)

ind <- which(rownames(x) %in% c("INDPRO", "CPIAUCSL"))

for(tt in m:(n - 24)){

  int <- (tt - m + 1):tt
  
  for(ll in 1:dim(err)[[4]]){
    
    if(ll == 1){
      out <- rob.tfa(x[, int], r = r, standardize = 'mean', 
                     tau = NULL, kappa = NULL, nfold = 3)
      max.tau <- as.numeric(dimnames(out$tau.cv)[[1]])[1]
      trunc.param.seq[tt] <- out$tau
      
      loading <- out$loading[[1]]
      eigval <- out$eigval[[1]][1:r]
      ff <- out$f[, dim(out$f)[2]]
      center <- out$center
      sc <- out$sc
      
      proj <- loading %*% t(loading) / p
      tx <- trunc.tnsr((x[, int] - center) / sc, trunc.param.seq[tt])
    } else if(ll == 2){
      out <- rob.tfa(x[, int], r = r, standardize = 'mean', 
                     tau = max.tau * 1.1, kappa = NULL, nfold = 3)
      
      loading <- out$loading[[1]]
      eigval <- out$eigval[[1]][1:r]
      ff <- out$f[, dim(out$f)[2]]
      center <- out$center
      sc <- out$sc
      
      proj <- loading %*% t(loading) / p
      tx <- (x[, int] - center) / sc
    } else if(ll == 3){
      center <- out$center
      sc <- out$sc
      tx <- (x[, int] - center) / sc
      
      tryCatch(hpca <- HDRFA::HPCA(t(tx), Method = 'E', r = r), warning = function(w){ print(tt) })
      
      loading <- hpca$Lhat
      eigval <- diag(t(hpca$Fhat) %*% hpca$Fhat / m)
      ff <- hpca$Fhat[dim(hpca$Fhat)[1], ]
      proj <- loading %*% t(loading) / p
      
    } 
    
    for(hh in 1:24){
      h <- hh
      
      z <- (x[, tt + h, drop = FALSE] - center) / sc
      Gamma_x <- tx[, 1:(m - h)] %*% t(tx[, 1:(m - h) + h])/m
      
      if(ll == 1){
        fc <- (proj %*% trunc.tnsr(z, trunc.param.seq[tt])) * sc + center
        err[, tt, 1, ll, hh] <- fc[ind] - x[ind, tt + h, drop = FALSE]
      } else{
        fc <- (proj %*% z) * sc + center
        err[, tt, 1, ll, hh] <- fc[ind] - x[ind, tt + h, drop = FALSE]
      } 
      
      fc <- proj %*% t(Gamma_x) %*% t(t(loading) /eigval) %*% ff
      fc <- fc * sc + center
      err[, tt, 2, ll, hh] <- fc[ind] - x[ind, tt + h, drop = FALSE]
    }
  }
}

ls <- list(err = err, trunc.param.seq = trunc.param.seq)
save(ls, file = paste("~/downloads/fredmd_forecasting_r", r, ".RData", sep = ""))

## results / run fluctation_test

jj <- 2 # 1: in-sample estimation 2: forecasting error
h <- 24 # number of lags to examine
kk <- 2 #  1: INDPRO 2: CPIAUCSL
nm <- c("INDPRO", "CPIAUCSL")[kk]

ll <- 1 # trunc
mm <- 2 # pca

loss1 <- c(apply(abs(ls$err[kk,, jj, ll, 1:h, drop = FALSE]), c(1, 2, 3, 4), mean))[-c(1:m, n - 1:24 + 1)] 
loss2 <- c(apply(abs(ls$err[kk,, jj, mm, 1:h, drop = FALSE]), c(1, 2, 3, 4), mean))[-c(1:m, n - 1:24 + 1)]

par(mfrow = c(1, 1), mar = c(2.5, 2.5, 2.5, .5))
plot(x = 1:(n - 24 - m), y = loss1 - loss2, 
     bty = "n", ylab = "", xlab = "Time", 
     type = "l", col = "cornflowerblue", lwd = 2.5, axes = FALSE,
     main = paste(nm, ': difference in the loss', sep = ''))
axis(2)
abline(h = 0, lwd = 1.8, lty = 2)
time_labels <- data$date[(m + 1):(n - 24)]
inds <- floor(seq(from = 1, to = length(time_labels), length.out = 5))
axis(1, at = inds, labels = time_labels[inds], cex.axis = 1)

for(mm in 1:3){
  mu <- c(.2, .3, .4)[mm]
  ft <- fluctuation_test(loss1, loss2, mu = mu, lag_truncate = 0, 
                         conf_level = .1, dmv_fullsample = TRUE,
                         time_labels = data$date[(m + 1):(n - 24)], 
                         main = paste(nm, ': mu = ', mu, sep = ''))
  # print(ft$df$time[ft$df$dmstat < ft$CV[1]])
  # print(ft$df$time[ft$df$dmstat > ft$CV[2]])
}
