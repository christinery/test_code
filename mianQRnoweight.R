############################################################################
###This is the main code for the paper "Copula-based M-estimate for FVCM"###
###under the quantile level without weight for the objective function#######
############################################################################

rm(list = ls())
set.seed(123456789)
tstart = Sys.time()
setwd("/Users/wangyafei/Dropbox/paper_UA/proposal 1_wyf/modularcode/QR")
library(MASS)
library(mvtnorm)
n = 100
tau.fix = 0.5
m = 50
p = 3
sgrid = (0:(m-1))/(m-1) ##grid point
lambda = 0.1
eta = 2
S = matrix(1,n,m)
B = 200
REP = 500
##"BETATOTAL" store all \hat\beta(s) with num. REP
BETATOTAL = array(0, c(p, m, REP))

##"BETABOOTVAR" used to store variance of \hat\beta^(sim)(s) with sample size B
BETABOOTVAR = array(0, c(p, m, REP))

for(sim in 1:REP){
  source("dataQR.R")
  source("algorithmQR.R")
  
  ################################################
  ###estimate \beta(s) under no weight case#######
  ################################################
  beta.ini = update.fun(y, x, tau.fix, eta, lambda, phi, rk, S)
  BETATOTAL[,,sim] = beta.ini
  
  ###############################################################################
  ###generate bootstrap sample and give the corresponding estimate of \beta(s)###
  ###within the sub sample#######################################################
  ##BETABOOTSTRAP used to store all \hat\beta(s) under bootstrap sample for REP=sim
  BETABOOTSTRAP = array(0, c(p, m, B))
  ##set betabootsum to compute mean of \hat\beta(s) with the sum times B
  betabootsum = matrix(0, p, m)
  for (l in 1:B) {
    ##the sub sample index 
    subsampleindex = sample(1:n, n, replace = TRUE)
    ##based on the index to find the resulting sub sample
    yb = y[subsampleindex,]
    xb = x[subsampleindex,]
    betabootstrap = update.fun(yb, xb, tau.fix, eta, lambda, phi, rk, S)
    betabootsum = betabootsum + betabootstrap 
    BETABOOTSTRAP[,,l] = betabootstrap
  }
  betabootmean = betabootsum/B ##p*m matrix
  source("arraycomfun.R")
  ##"betabootvar" is actually the variance of \hat\beta^(sim)(s), return p*m matrix
  betabootvar = matrixdiff(BETABOOTSTRAP, betabootmean)
  BETABOOTVAR[,,sim] = betabootvar
} 
##################################################################################
###the end of using bootstrap method to estimate variance of \hat\beta^{sim}(s)###
##################################################################################

##################################################################################
### compute different measure value ##############################################
##################################################################################
  source("arraycomfun.R")

  ##the average value for \hat\beta(s), return p*m matrix
  BETAMEAN = matrixsqmean(BETATOTAL, 1)
  
  ##compute bias for \hat\beta(s), return p*m matrix
  BETABIAS = matrixsqmean(BETATOTAL, 1) - beta.true
  write.table(BETABIAS, "BIASQR.txt")
  
  ##compute standard error, return p*m matrix
  BETASE = sqrt(matrixdiff(BETATOTAL, BETAMEAN))
  write.table(BETASE, "SEQR.txt")
  
  ##compute standard deviation, return p*m matrix
  BETASD = matrixsqmean(BETABOOTVAR, 0)
  write.table(BETASD, "SDQR.txt")
  
  ##compute the upper and low value for 95% confidence interval, return p*m matrix
  ##firstly compute all(num.REP) CI (lower and upper) 
  cil = BETATOTAL - 1.96*sqrt(BETABOOTVAR)
  cir = BETATOTAL + 1.96*sqrt(BETABOOTVAR)
  BETACIL = matrixsqmean(cil, 1) ##lower bound
  BETACIR = matrixsqmean(cir, 1) ##upper bound
  
  
  ##compute empirical coverage probability, return p*m matrix
  BETAECP = matrixecp(BETATOTAL, BETABOOTVAR, beta.true)
  write.table(BETAECP, "ECPQR.txt")
  
  ##compute root mean integrated square error, return 1*p vector
  BETA = array(beta.true, c(p, m, REP))
  ##for each REP=sim, compute "((REP=sim)-(beta.true))^2"
  BETADIFSQ = (BETATOTAL - BETA)^2
  ##for each "((REP=sim)-(beta.true))^2", compute the sqrt mean of row sum, 
  ##and then combine these num. "REP" vector
  COMREP = sqrt(arraymatrixmean(BETADIFSQ))
  RMISE = apply(COMREP, 2, mean)
  write.table(RMISE, "RMISE.txt")
  ##########################################################################
  ############################# end ########################################
  ##########################################################################
  
  
  ##########################################################################
  ### plot kinds of figures ################################################
  ##########################################################################
  ##plot the \hat\beta(s)
  pdf("betahat.pdf")
  par(mfrow = c(1,p))
  for (k in 1:p){
    plot(beta.true[k,],type="l")
    lines(BETAMEAN[k,],type = "l",col="red")
  }
  dev.off()
  
  ##plot \hat\beta bias
  pdf("bias.pdf")
  par(mfrow = c(1,p))
  for (k in 1:p){
    plot(BETABIAS[k,],type="p")
  }
  dev.off()
  
  ##plot SE and SD
  pdf("SDE.pdf")
  par(mfrow = c(1, p))
  for (k in 1:p){
    plot(BETASE[k,],type="p", col = "blue")
    lines(BETASD[k,],type = "p",col="red")
  }
  dev.off()
  
  ##plot upper and lower bound line
  pdf("CI.pdf")
  par(mfrow = c(1, p))
  for (k in 1:p){
    plot(beta.true[k,],type="l")
    lines(BETACIL[k,], type = "l", col = "red")
    lines(BETACIR[k,], type = "l", col = "red")
  }
  dev.off()
  
  ##plot ECP
  pdf("ECP.pdf")
  par(mfrow = c(1, p))
  for (k in 1:p){
    plot(BETAECP[k,], type = "p")
    abline(h = 0.95, col = "red", lty=1, lwd=3)
  }
  dev.off()
  
  
  
  
  
  
  
  
  
  
  