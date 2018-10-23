############################################################################
###This file is used to generate random data for the paper "Copula-based"###
###M-estimate for FVCM"#####################################################
############################################################################

library(MASS)

## basis of null space and reproducing kernel of RKHS
phi = cbind(rep(1/m,m),sgrid)
rk = matrix(0,m,m)
for (j in 1:m) {
  rk[j,] = exp(-(sgrid[j]-sgrid)^2/(0.5)^2)
}

## set the true coefficient parameter \beta  p*m matrix
beta_1 = function(s){
  (1-s)^2
}
beta_2 = function(s){
  4*s*(1-s)
}
beta_3 = function(s){
  sin(s) + s^3
}
beta.true = rbind(beta_1(sgrid), beta_2(sgrid), beta_3(sgrid))
p = nrow(beta.true)

## generate covariate terms
x1 = rbinom(n, 1, 0.5)
x2 = rnorm(n, 0, 1)
x = cbind(1, x1, x2)

## generate random number for error term
mu1 = rep(0, m)
sigm.cov1 = matrix(0,m,m)
ro = 0.3
for (k in 1:(m-1)) {
  for (l in (k+1):m) {
    sigm.cov1[k,l] = ro^{l-k}
  }
}
sigm.cov1 = sigm.cov1 + t(sigm.cov1) + diag(1, m, m)
error1 = mvrnorm(n, mu1, sigm.cov1, empirical = TRUE)
Error = t(t(as.matrix(error1))-apply(error1, 2, function(x){quantile(x,tau.fix)}))

## generate response y  n*m matrix
y = x%*%beta.true+Error
