##########################################################################
###This file is used to show adopted algorithm for paper "Copula-based"###
###M-estimate for FVCM"###################################################
##########################################################################

## quantile loss function
quan.lf = function(tau,x){
  x*(tau-ifelse(x<0,1,0))
}

## compute all g_{ij},f_{ij} terms
gf.ini = function(x,phi,rk){
  ##x: n*p matrix covariate variable
  ##phi: m*n_0 matrix basis (number is n_0) of null space of RKHS
  ##rk: m*m matrix discreted reproducing kenel of RKHS
  g.com = f.com = NULL 
  for (i in 1:nrow(x)){
    g = kronecker(as.matrix(x[i,]),t(phi)) ##(pn_0)*m matrix constant term of coefficient b
    f = kronecker(as.matrix(x[i,]),rk) ##(pm)*m matrix contant term of coefficient c
    g.com = cbind(g.com,g) ##combine all term g_{ij} (p*n0)*(m*n) matrix
    f.com = cbind(f.com,f) ##combine all term f_{ij} (p*m)*(m*n) matrix
  }
  return(rbind(g.com,f.com))
}

## compute term related to covariate i.e. g_{ij}^Tb,f_{ij}^Tc
gf.fun = function(x,phi,rk,b,c){
  ##x: n*p matrix covariate variable
  ##phi: m*n_0 matrix basis (number is n_0) of null space of RKHS
  ##rk: m*m matrix discreted reproducing kenel of RKHS
  ##b,c: coefficients need to be estimated
  G = F1 = NULL
  for (i in 1:nrow(x)) {
    g = kronecker(as.matrix(x[i,]),t(phi)) ##(pn_0)*m matrix constant term of coefficient b
    f = kronecker(as.matrix(x[i,]),rk) ##(pm)*m matrix contant term of coefficient c
    g.tran = crossprod(g,b) ##m*1 matrix term g_{ij}^T %*%b
    f.tran = crossprod(f,c) ##m*1 matrix term f_{ij}^T%*%c
    G = cbind(G,g.tran) ##m*n matrix (g_{ji}^T b)_{j=1,...m;i=1,...,n}
    F1 = cbind(F1,f.tran) ##m*n matrix (f_{ij}^T c)_{j=1,...m;i=1,...,n}
  }
  gf.com = rbind(G,F1)
  return(gf.com) ##return value G and F arrange by row 2m*n matrix
}

## shrinkage operator 
shrin.ope = function(A,tau,u,lambda,v){
  ifelse(v>(u+lambda*tau*A),v-lambda*tau*A,
         ifelse(v<(u-lambda*(1-tau)*A),v+lambda*(1-tau)*A,0))
}

update.fun = function(y,x,tau,eta,lambda,phi,rk,S){
  ##y: n*m matrix functional response
  ##x: n*p matrix covariate variable
  ##tau: quantile leval
  ##eta: augmented lagrangian parameter
  ##lambda: penalty parameter
  ##phi: m*n_0 matrix basis (number is n_0) of null space of RKHS
  ##rk: m*m matrix discreted reproducing kenel of RKHS
  ##S: n*m matrix constant for shrinkage operator 
  max.iter = 500
  abs.error = 10^(-4) ##absolute error
  rel.error = 10^(-2) ##relative error
  n = nrow(y)
  p = ncol(x)
  n0 = ncol(phi)
  m = ncol(y)
  sigma.star = kronecker(diag(1,p),rk)
  GF = gf.ini(x,phi,rk)
  G.ini = GF[1:(p*n0),] ##all g_{ij} term (p*n0)*(n*m) matrix
  F.ini = GF[(p*n0+1):nrow(GF),] ##all f_{ij} term (p*m)*(n*m) matrix
  g.sum = matrix(0,n0*p,n0*p)
  f.sum = matrix(0,p*m,p*m)
  for (i in 1:n) {
    g.coe = G.ini[,((i-1)*m+1):(i*m)] ##(p*n0)*m matrix
    f.coe = F.ini[,((i-1)*m+1):(i*m)] ##(p*m)*m matrix
    g.sum = g.sum+g.coe%*%t(g.coe) 
    f.sum = f.sum+f.coe%*%t(f.coe)
  }
  f.sum = 2/eta*lambda*sigma.star+f.sum +0.01*diag(1,p*m,p*m)
  
  ## update step
  iter = 1
  b.iter = rep(0,n0*p)
  c.iter = rep(0,p*m)
  U.iter = matrix(0,n,m)
  W.iter = matrix(0,n,m) ##initial value of w_{ij}
  while (iter<=max.iter) {
    #if(iter == 63) browser()
    G0 = gf.fun(x,phi,rk,b.iter,c.iter)[1:m,] ##m*n matrix (g_{ji}^T b)_{j=1,...m;i=1,...,n}
    F0 = gf.fun(x,phi,rk,b.iter,c.iter)[(m+1):(2*m),] ##m*n matrix (f_{ij}^T c)_{j=1,...m;i=1,...,n}
    G = t(as.matrix(G0)) ## n*m matrix 
    F1 = t(as.matrix(F0)) ## n*m matrix
    ## update u_{ij}
    for (i in 1:n) {
      for (j in 1:m) {
        v = G[i,j]+F1[i,j]-W.iter[i,j]
        U.iter[i,j] = shrin.ope(S[i,j],1-tau,y[i,j],1/eta,v)
      }
    }
    
    for (b in 1:3) {
      ## update b
      b.pre = b.iter
      b1 = matrix(0,n0*p,1)
      for (i in 1:n) {
        b1 = b1+apply(t(t(G.ini[,((i-1)*m+1):(i*m)])*(as.vector(U.iter[i,]-F1[i,]+W.iter[i,]))),1,sum)
      }
      b.iter = solve(g.sum)%*%b1
      
      ## update c
      c.pre = c.iter
      G0 = gf.fun(x,phi,rk,b.iter,c.iter)[1:m,] ##m*n matrix (g_{ji}^T b)_{j=1,...m;i=1,...,n}
      G = t(as.matrix(G0)) ## n*m matrix 
      c1 = matrix(0,p*m,1)
      for (i in 1:n) {
        c1 = c1+apply(t(t(F.ini[,((i-1)*m+1):(i*m)])*(as.vector(U.iter[i,]-G[i,]+W.iter[i,]))),1,sum)
      }
      c.iter = solve(f.sum)%*%c1
    }
    ## update w_{ij}
    G0 = gf.fun(x,phi,rk,b.iter,c.iter)[1:m,] ##m*n matrix (g_{ji}^T b)_{j=1,...m;i=1,...,n}
    F0 = gf.fun(x,phi,rk,b.iter,c.iter)[(m+1):(2*m),] ##m*n matrix (f_{ij}^T c)_{j=1,...m;i=1,...,n}
    G = t(as.matrix(G0)) ## n*m matrix 
    F1 = t(as.matrix(F0)) ## n*m matrix
    for (i in 1:n) {
      for (j in 1:m) {
        W.iter[i,j] = W.iter[i,j]+(U.iter[i,j]-G[i,j]-F1[i,j])
      }
    }
    ## compute the estimate of functional coefficient \beta(s) based on \hat b and \hat c
    ## discreted form \beta(s_j), return p*m matrix
    betaes = kronecker(diag(1,p,p),phi)%*%b.iter+kronecker(diag(1,p,p),rk)%*%c.iter
    betatilde = t(matrix(betaes,m,p))
    
    ## compute primal and dual residual
    G = gf.fun(x,phi,rk,b.iter,c.iter)[1:m,] ##m*n matrix (g_{ji}^T b)_{j=1,...m;i=1,...,n}
    F1 = gf.fun(x,phi,rk,b.iter,c.iter)[(m+1):(2*m),] ##m*n matrix (f_{ij}^T c)_{j=1,...m;i=1,...,n}
    F.pre = gf.fun(x,phi,rk,b.iter,c.pre)[(m+1):(2*m),]
    G.pre = gf.fun(x,phi,rk,b.pre,c.iter)[1:m,]
    ## primal residual
    pri.res = sqrt(mean((U.iter-t(as.matrix(F1)))^2))
    ## dual residual
    #dua.res = sqrt(eta^2*mean((F1-F.pre)^2))
    dua.res = eta*sqrt(mean((F1-F.pre)+(G-G.pre))^2)
    #cat("dua.res=",dua.res)
    
    ## stopping criterion
    err.pri = sqrt(m*n)*abs.error+rel.error*max(sqrt(mean(U.iter^2)),sqrt(mean(F1^2)))
    err.dua = sqrt(m*n)*abs.error+rel.error*sqrt(mean(W.iter^2))
    #if(iter == 41) browser()
    #if(!is.finite(err.pri)) browser()
    #if(is.na(pri.res<=err.pri & dua.res<=err.dua)) browser()
    if(pri.res<=err.pri & dua.res<=err.dua) break
    iter=iter+1
  }
  return(betatilde)
}


