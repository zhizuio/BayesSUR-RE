#===================================================================================================
# This script is to simulate data based on a SUR model without random effects
#
# Author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# Date: 01-June-2022
#===================================================================================================

sim.ssur <- function(n, s, p, t0=0, seed=123, mv=TRUE){
  # set seed to fix coefficients
  set.seed(7193) 
  sd_b = 1
  mu_b = 1
  b = matrix(rnorm((p+ifelse(t0==0,1,0))*s,mu_b,sd_b),p+ifelse(t0==0,1,0),s)
  
  # design groups and pathways of Gamma matrix
  gamma = matrix(FALSE,p+ifelse(t0==0,1,0),s)
  if(t0==0) gamma[1,] = TRUE
  gamma[2:6-ifelse(t0==0,0,1), 1:5] = TRUE
  gamma[11:21-ifelse(t0==0,0,1), 6:12] = TRUE
  gamma[31:51-ifelse(t0==0,0,1), 1:5] = TRUE
  gamma[31:51-ifelse(t0==0,0,1), 13:20] = TRUE
  gamma[52:61-ifelse(t0==0,0,1), 1:12] = TRUE
  gamma[71:91-ifelse(t0==0,0,1), 6:20] = TRUE
  gamma[111:121-ifelse(t0==0,0,1), 1:20] = TRUE
  
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
  
  #x = (x-matrix( colMeans(x), n,p,byrow=TRUE)) %*% diag(1/apply(x,2,sd))
  if(t0==0) x = cbind(rep(1,n),x)
  
  xb = matrix(NA,n,s)
  
  if(mv){
    for(i in 1:s){
      if(sum(gamma[,i])>1){
        xb[,i] = x[,gamma[,i]] %*% b[gamma[,i],i] 
      }else{
        xb[,i] = sapply(1:s,function(i) rep(1,n) * b[1,i])
      }
    }
  }else{
    if(sum(gamma)>1){
      xb = x[,gamma] %*% b[gamma,]
    }else{
      xb = sapply(1:s,function(i) rep(1,n) * b[1,i])
    }
  }
  
  corr_param = 0.9 # in 0.3 , 0.6 , 0.9
  M = matrix(corr_param,s,s)
  diag(M) = rep(1,s)
  
  ## wanna make it decomposable
  Prime = list(c(1:(s*.4),(s*.8):s),c((s*.4):(s*.6)),c((s*.65):(s*.75)),c((s*.8):s))
  G = matrix(0,s,s)
  for(i in 1:length(Prime))
    G[Prime[[i]],Prime[[i]]] = 1  
  
  #check
  dimnames(G) = list(1:s,1:s)
  length( gRbase::mcsMAT(G - diag(s)) ) > 0
  
  var = solve(BDgraph::rgwish(n=1,adj=G,b=3,D=M))
  
  # change seeds to add randomness on error
  set.seed(seed+8493)
  sd_err = 0.5 
  err = matrix(rnorm(n*s, 0, sd_err),n,s) %*% chol(as.matrix(var))
  
  if(t0 == 0){
    b.re = NA
    z = NA
    y = xb + err
    return( list(y=y, x=x, b=b, gamma=gamma, z=z, b.re=b.re, Gy=G, mrfG=combnAll) )
  }else{
    # add random effects
#    z = model.matrix( ~ factor(rMultinom(matrix(c(.1,.2,.3,.4),nrow=1),n)) - 1)
    z = t(rmultinom(n, size = 1, prob = c(.1,.2,.3,.4)))
    # use fixed random seed for b.re
#    set.seed(1683)
#    b.re = matrix(rnorm(t0*s,0.5,sd_b), t0, s)
#    b.re = matrix(rnorm(t0*s,0,2), t0, s)
#    y = z%*%b.re + xb + err
    z <- sample(1:t0, n, replace=T, prob=rep(1/t0,t0))
    set.seed(1683)
    b.re <- rnorm(t0, 0, 2)
    y <- matrix(b.re[z],nrow=n,ncol=s) + xb + err

    return( list(y=y, x=x, b=b, gamma=gamma, z=model.matrix(~ factor(z) + 0)[,],
                 b.re=b.re, Gy=G, mrfG=combnAll) )
  }

}
