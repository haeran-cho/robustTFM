source('codes.R')

library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)

# only for comparison of factor number estimator
library(TensorPreAve)
library(HDRFA)

load('euro_HT.RData')

# data
x <- ls$x
dimnames(x)[[2]]

pn <- dim(x)
K <- length(pn) - 1
p <- pn[1:K]
n <- pn[K + 1]

center <- apply(x, 1:K, median)
sc <- apply(x, 1:K, function(z){ ifelse(mad(z) > 0, mad(z), sd(z)) })
sx <- apply(x, 1:K, function(z){ z <- z - median(z); z / ifelse(mad(z) > 0, mad(z), sd(z)) })
sx <- aperm(sx, c(2:(K + 1), 1))

## number of factors

# only for comparison 
ax <- aperm(sx, c(3, 1, 2))

proj <- TensorPreAve::rank_factors_est(rTensor::as.tensor(ax))
proj$rank

rtfa <- RTFA::TFM_FN(ax, method = 'HUBER')
rtfa$factor.num
#

tau.seq <- exp(seq(log(max(abs(sx))), log(median(abs(sx))), length.out = 50))
tau.ind <- unique(c(1, which(- diff(tau.seq) > diff(range(sx)) * .001), 50))
tau.seq <- tau.seq[tau.ind]

r.hat.array <- array(0, dim = c(length(tau.seq), K))
for(ll in 1:length(tau.seq)){  
  tau <- tau.seq[ll]
  rtfa <- rob.tfa(sx, r = NULL, standardize = 'none', tau = tau)
  r.hat.array[ll, ] <- rtfa$r
}

par(mfrow = c(1, 2), mar = c(4, 4, 2, .5))
for(kk in 1:2) plot(tau.seq, r.hat.array[, kk], type = 'b', 
                    xlab = expression(tau), log = "x", ylab = expression(r), main = paste("k = ", kk, sep = ""))
dev.off()

r <- c(1, 3)
out <- rob.tfa(sx, r = r, standardize = 'none', nfold = 3)
out$tau

par(mar = c(4, 4, .5, .5))
plot(as.numeric(dimnames(out$tau.cv)[[1]]), apply(out$tau.cv[,2,,], 1, mean), type = 'b', 
     log = "x", xlab = expression(tau), ylab = "CV")

out$f <- - out$f
out$loading[[1]] <- - out$loading[[1]]

# plotting f

f <- aperm(out$f, c(3, 1, 2))
df <- as.data.frame(f)
names(df) <- c("f.1", "f.2", "f.3")[1:r[2]]
df$time <- as.Date(ls$time)

ggplot(df, aes(x = time, y = f.1)) + geom_line() +
  xlab("time") + ylab("") +
  theme_minimal() + 
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) -> plt1

ggplot(df, aes(x = time, y = f.2)) + geom_line() +
  xlab("time") + ylab("") +
  theme_minimal() + 
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) -> plt2

ggplot(df, aes(x = time, y = f.3)) + geom_line() +
  xlab("time") + ylab("") +
  theme_minimal() + 
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) -> plt3

ga <- grid.arrange(plt1, plt2, plt3, layout_matrix = matrix(c(1, 2, 3), nrow = 3))

# plotting loading

rownames(out$loading[[1]]) <- dimnames(x)[[1]]
rownames(out$loading[[2]]) <- dimnames(x)[[2]]

rn <- range(out$loading)

df1 <- as.data.frame(as.table(out$loading[[1]]))
names(df1) <- c('country', 'j', 'value')
levels(df1$j) <- 1

ggplot(df1, aes(x = j, y = country, fill = value)) +
  geom_tile() + 
  scale_fill_viridis_c(option = 'H', limits = rn) +
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.title.y = element_blank(), legend.position = "none") -> plt1

plt1

vr <- varimax(out$loading[[2]], FALSE)
vr$rotmat

df2 <- as.data.frame(as.table(vr$loading))

names(df2) <- c('indicator', 'j', 'value')
levels(df2$j) <- 1:nlevels(df2$j)

ggplot(df2, aes(x = j, y = indicator, fill = value)) +
  geom_tile() + 
  scale_fill_viridis_c(option = 'H', limits = rn) +
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0, hjust = 1),
        axis.title.y = element_blank(), legend.position = "right") -> plt2

plt2

ga <- grid.arrange(plt1, plt2, layout_matrix = matrix(c(1, 2, 2), nrow = 1))

## plot loading kronecker

kk <- 1; jj <- 2

ll <- kronecker(out$loading[[1]], vr$loading, make.dimnames = TRUE)

dimnames(ll)[[1]][(1:p[kk] - 1) * p[jj] + 1] <- dimnames(x)[[kk]]
dimnames(ll)[[2]] <- 1:3

ldf <- as.data.frame(as.table(ll))
names(ldf) <- c('both', 'j', 'value')

ggplot(ldf, aes(x = j, y = both, fill = value)) +
  geom_tile() + 
  scale_fill_viridis_c(option = 'H') +
  scale_y_discrete(breaks = function(x){x[c(TRUE, rep(FALSE, p[jj] - 1))]}) +
  theme(strip.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "right") -> plt

plt
