#===================================================================================================
# This script is to simulate data based on a SUR model with random effects
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 08-June-2021
#===================================================================================================

sim.ssur.re <- function(n, s, p, t0=0, seed=123){
  # set seed to fix coefficients
  set.seed(7193) 
  sd_b = 1
  mu_b = 1
  b = matrix(rnorm(p*s,mu_b,sd_b),p,s)
  #b = matrix(runif(p*s,-2,2),p,s)
  
  snr = 10
  # design groups and pathways of Gamma matrix
  gamma = matrix(FALSE,p,s)
  #gamma[1,] = TRUE
  gamma[1:5, 1:5] = TRUE
  gamma[10:20, 6:12] = TRUE
  gamma[30:50, 1:5] = TRUE
  gamma[30:50, 13:20] = TRUE
  gamma[51:60, 1:12] = TRUE
  gamma[70:90, 6:20] = TRUE
  gamma[110:120, 1:20] = TRUE
  
  G_kron = matrix(0,s*p,s*p)
  G_m = bdiag(matrix(1,ncol=5,nrow=5), matrix(1,ncol=7,nrow=7), matrix(1,ncol=8,nrow=8))
  G_p = bdiag(matrix(1,ncol=5,nrow=5), diag(3), matrix(1,ncol=11,nrow=11), diag(9), matrix(1,ncol=21,nrow=21), matrix(1,ncol=10,nrow=10), diag(9), matrix(1,ncol=21,nrow=21), diag(19), matrix(1,ncol=11,nrow=11), diag(181))
  G_kron = kronecker(G_m, G_p)
  
  combn11 = combn(rep((1:5-1)*p,each=length(1:5)) + rep(1:5,times=length(1:5)), 2)
  combn12 = combn(rep((1:5-1)*p,each=length(30:60)) + rep(30:60,times=length(1:5)), 2)
  combn13 = combn(rep((1:5-1)*p,each=length(110:120)) + rep(110:120,times=length(1:5)), 2)
  combn21 = combn(rep((6:12-1)*p,each=length(10:20)) + rep(10:20,times=length(6:12)), 2)
  combn22 = combn(rep((6:12-1)*p,each=length(51:60)) + rep(51:60,times=length(6:12)), 2)
  combn23 = combn(rep((6:12-1)*p,each=length(70:90)) + rep(70:90,times=length(6:12)), 2)
  combn24 = combn(rep((6:12-1)*p,each=length(110:120)) + rep(110:120,times=length(6:12)), 2)
  combn31 = combn(rep((13:20-1)*p,each=length(30:50)) + rep(30:50,times=length(13:20)), 2)
  combn32 = combn(rep((13:20-1)*p,each=length(70:90)) + rep(70:90,times=length(13:20)), 2)
  combn33 = combn(rep((13:20-1)*p,each=length(110:120)) + rep(110:120,times=length(13:20)), 2)
  
  combnAll = rbind(t(combn11), t(combn12), t(combn13), t(combn21), t(combn22), t(combn23), t(combn24), t(combn31), t(combn32), t(combn33))
  
  set.seed(seed+7284)
  sd_x=1  # I later re-scale it so it doesn't matter that much...
  x = matrix(rnorm(n*p,0,sd_x),n,p)
  
  xb = matrix(NA,n,s)
  
  for(i in 1:s){
    if(sum(gamma[,i])>1){
      xb[,i] = x[,gamma[,i]] %*% b[gamma[,i],i] 
    }else{
      xb[,i] = sapply(1:s,function(i) rep(1,n) * b[1,i])
    }
  }
  
  corr_param = 0.5 # in 0.3 , 0.6 , 0.9
  M = matrix(corr_param,s,s)
  diag(M) = rep(1,s)
  
  ## wanna make it decomposable
  Prime = list(c(1:(s*.4),(s*.8):s),c((s*.4):(s*.6)),c((s*.65):(s*.75)),c((s*.8):s))
  G = matrix(0,s,s)
  for(i in 1:length(Prime))
    G[Prime[[i]],Prime[[i]]] = 1  
  
  #check
  dimnames(G) = list(1:s,1:s)
  #length( gRbase::mcsMAT(G - diag(s)) ) > 0
  
  var = solve(BDgraph::rgwish(n=1,adj=G,b=3,D=M))
  
  #control signal-to-noise ratio
  factor = 10 ; factor_min = 0.01; factor_max = 1000
  count = 0 ; maxit = 10000
  
  factor_prev = 1
  
  set.seed(seed+8493)
  repeat{
    
    var = var / factor * factor_prev
    
    ### Sample the errors and the Ys
    cVar = chol(as.matrix(var))
    err = matrix(rnorm(n*s),n,s) %*% cVar
    
    if(t0 == 0){
      b.re <- 0
      z <- 0
      y <- xb + err
    }else{
      # add random effects
      z <- t(rmultinom(n, size = 1, prob = c(.1,.2,.3,.4)))
      # use fixed random seed for b.re
      b.re <- matrix(rnorm(t0*s,0,2), t0, s)
      y <- z%*%b.re + xb + err
    }
    
    ## Reparametrisation ( assuming PEO is 1:s )
    cVar = t(cVar) # make it lower-tri
    S = diag(diag(cVar))
    sigma = S*S
    L = cVar %*% solve(S)
    rho =  diag(s) - solve(L)
    
    ### S/N Ratio
    emp_snr = mean( diag( var(xb) %*% solve(sigma) ))
    emp_g_snr = mean( diag( var( (err)%*%t(rho) ) %*% solve(sigma) ))
    
    ##############
    
    if( abs(emp_snr - snr) < (snr/10) | count > maxit ){
      break
    }else{
      if( emp_snr < snr ){ # increase factor
        factor_min = factor
      }else{ # decrease factor
        factor_max = factor
      }
      factor_prev = factor
      factor = (factor_min + factor_max)/2
    }
    count = count+1
  }
  return( list(y=y, x=x, b=b, gamma=gamma, z=z, b.re=b.re, Gy=G, mrfG=combnAll ) )
  
}
