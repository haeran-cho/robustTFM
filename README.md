# robustTFM
Codes accompanying "Tail-robust factor modelling of vector and tensor time series in high dimensions" by Barigozzi, Cho and Maeng (2024).

## file description

### codes.R
Contains the main routine `rob.tfa` for the proposed robust factor analysis.

### dgp.R
Contains the codes for generating the data used in simulation studies.

### fredmd.csv
Contains the FREMD-MD data analysed in Section 6.1; accessed from https://research.stlouisfed.org/econ/mccracken/fred-databases/ on 01-03-2024.
We provide the csv file to ensure consistency with the data analysis reported in the paper.

### fredmd_submit.R
Contains codes for reproducing the rolling window-based forecasting results in Section 6.1

### euro_HT.RData
Contains the EA-MD data analysed in Section 6.2.

### euro_submit.R
Contains codes for reproducing the results in Section 6.2

## usage example

```
source('codes.R')
source('dgp.R')

td <- tensor_dgp(n, case = 2, dist = 't', perc = .05)
rtfa <- rob.tfa(td$x, r = c(3, 3, 3))

plot(as.numeric(dimnames(rtfa$tau.cv)[[1]]), apply(rtfa$tau.cv[,2,,], 1, mean), type = 'b', 
     log = "x", xlab = expression(tau), ylab = "CV")
abline(v = rtfa$tau, col = 2)
```
