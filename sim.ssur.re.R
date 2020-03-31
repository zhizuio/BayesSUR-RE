
library(Matrix, lib.loc="/home/zhiz/RLibrary/")
library(coda, lib.loc="/home/zhiz/RLibrary/")
library(MCMCpack, lib.loc="/home/zhiz/RLibrary/")
library(gRbase, lib.loc="/home/zhiz/RLibrary/")
library(BDgraph, lib.loc="/home/zhiz/RLibrary/")
library(Formula, lib.loc="/home/zhiz/RLibrary/")
library(Hmisc, lib.loc="/home/zhiz/RLibrary/")
sim.ssur <- function(n, s, p, t0=0, snr=0, seed=123, mv1=FALSE, mv2=FALSE, mv3=FALSE, mv4=TRUE){
  # set seed to fix coefficients
  set.seed(seed) 
  sd_b = .2
  b = matrix(rnorm((p)*s,0.5,sd_b),p,s)
  
  ## generate sparsity 20% of gamma matrix, but randomly without specific patterns
  #spars_lvl = 0.2
  #gamma = c(TRUE, (1:p) %in% sample(1:p,p*spars_lvl,replace = FALSE) ) #c(TRUE,rep(FALSE,p)) 
  
  # sim2.1
  mv1 = mv1
  # sim2.2
  mv2 = mv2
  # sim4
  mv3 = mv3
  mv4 = mv4
  if(mv1){
    gamma = matrix(NA,p+1,s)
    for(i in 1:s)
      gamma[,i] = c(TRUE, (1:p) %in% sample(1:p,p*spars_lvl,replace = FALSE) ) #c(TRUE,rep(FALSE,p)) 
  }
  if(mv2){
    # design hotspot Gamma matrix
    gamma = matrix(FALSE,p+1,s)
    gamma[1,] = TRUE
    gamma[floor(p/20):floor(p/16),floor(s/7):floor(s/5)] = TRUE
    gamma[2+floor(p/15):floor(p/13),floor(s/3):floor(s/3)] = TRUE
    gamma[4+floor(p/12):floor(p/10.5),floor(s/7):floor(s/3)] = TRUE
    gamma[6+floor(p/11):floor(p/10),floor(s/2.5):floor(s/2.5)] = TRUE
    gamma[8+floor(p/11):floor(p/10),floor(s/10):floor(s/10)] = TRUE
    gamma[10+floor(p/9.3):floor(p/8.5),floor(s/s):floor(s/1)] = TRUE
    gamma[22+floor(p/8):floor(p/7.5),floor(s/s):floor(s/4)] = TRUE
    gamma[24+floor(p/8):floor(p/7.5),floor(s-1):floor(s/1)] = TRUE
    gamma[26+floor(p/7):floor(p/6.6),floor(s/7):floor(s/2.2)] = TRUE
    gamma[28+floor(p/6.3):floor(p/5.9),floor(s/2):floor(s/1.7)] = TRUE
    gamma[30+floor(p/5.7):floor(p/5.4),floor(s/1.6):floor(s/1)] = TRUE
    gamma[44+floor(p/4.8):floor(p/4.6),floor(s/1.4):floor(s/1)] = TRUE
    gamma[44+floor(p/4.8):floor(p/4.6),floor(s/s):floor(s/4)] = TRUE
    gamma[46+floor(p/4.45):floor(p/4.25),floor(s/3.5):floor(s/2.2)] = TRUE
    gamma[48+floor(p/4.16):floor(p/4),floor(s/2):floor(s/1.8)] = TRUE
    gamma[50+floor(p/3.8):floor(p/3.7),floor(s/1.7):floor(s/1.5)] = TRUE
    
    # summarize the corresponding G_m & G_p
    sig.loc <- which(gamma[-1,]!=0,arr.ind=T)
    combn.x <- combn(sort(unique(sig.loc[,1])), 2)
    combn.y <- combn(sort(unique(sig.loc[,2])), 2)
    G_m_loc <- matrix(NA, ncol=3, nrow=dim(combn.y)[2])
    G_m_loc[,1:2] <- t(combn.y)
    for(i in 1:dim(combn.y)[2]){
      temp <- combn.y[,i]
      # count the number of the common features between two response variables
      G_m_loc[i,3] <- sum(sig.loc[sig.loc[,2]==temp[1],1] %in% sig.loc[sig.loc[,2]==temp[2],1])
    }
    G_m_loc <- G_m_loc[G_m_loc[,3]!=0,]
    G_m_loc[,3] <- round(G_m_loc[,3]/max(G_m_loc[,3]), dig=3)
    
    G_p_loc <- matrix(NA, ncol=3, nrow=dim(combn.x)[2])
    G_p_loc[,1:2] <- t(combn.x)
    for(i in 1:dim(combn.x)[2]){
      temp <- combn.x[,i]
      # count the number of the common response variables between two features
      G_p_loc[i,3] <- sum(sig.loc[sig.loc[,1]==temp[1],2] %in% sig.loc[sig.loc[,1]==temp[2],2])
    }
    G_p_loc <- G_p_loc[G_p_loc[,3]!=0,]
    G_p_loc[,3] <- round(G_p_loc[,3]/max(G_p_loc[,3]), dig=3)
    
    GEdge <- cbind(rbind(G_m_loc,matrix(0,ncol=3,nrow=dim(G_p_loc)[1]-dim(G_m_loc)[1])), G_p_loc)
  }
  if(mv3){
    gamma = matrix(FALSE,p+1,s)
    gamma[1,] = TRUE
    gamma_v0 = as.matrix( read_table(file="/Users/zhiz/Documents/PhD_project/London/gamma_v.txt",col_names=FALSE))
    gamma_v = matrix(FALSE,p,s)
    gamma_v[gamma_v0 >= 0.1067] = TRUE
    gamma_v[gamma_v0 < 0.1067] = FALSE
    gamma_v = matrix(gamma_v, nrow=p, ncol=s)
    gamma[-1,] = gamma_v
    
  }
  if(mv4){
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
    # combn1_5 = combn(1:5, 2)
    # combn29_59 = combn(29:59, 2)
    # combn109_119 = combn(109:119, 2)
    # combn9_19 = combn(9:19, 2)
    # combn44_59 = combn(44:59, 2)
    # combn69_89 = combn(69:89, 2)
    # combn29_49 = combn(29:49, 2)
    
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
    ##combnAll = combnAll[!duplicated(data.frame(combnAll)), ]
    #for(i in 1:dim(combnAll)[1]) G_kron[combnAll[i,1], combnAll[i,2]] = 1
    
    ##image(G_kron,col=grey((0:32)/32), axes = FALSE,main=expression(G))
  }
  
  set.seed(seed+7284)
  sd_x=1  # I later re-scale it so it doesn't matter that much...
  x = matrix(rnorm(n*p,0,sd_x),n,p)
  
  #x = (x-matrix( colMeans(x), n,p,byrow=TRUE)) %*% diag(1/apply(x,2,sd))
  #x = cbind(rep(1,n),x)
  
  xb = matrix(NA,n,s)
  
  if(mv1 | mv2 | mv3 | mv4){
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
  
  v_r = mean(diag(var(xb))) / snr
  
  nu = s
  # sim3.1
  M = matrix(5,s,s) + .8* diag(s) 
  
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
  # change seeds to add randomness on error
  set.seed(seed+8493)
  sd_err = 1 # 2
  
  if(!snr){
    err = matrix(rnorm(n*s, 0, sd_err),n,s) %*% chol(as.matrix(var))
    
    if(t0 == 0){
      b.re <- 0
      z <- 0
      y <- xb + err
    }else{
      # add random effects
      z <- model.matrix( ~ factor(rMultinom(matrix(c(.1,.2,.3,.4),nrow=1),n)) - 1)
      # use fixed random seed for b.re
      b.re <- matrix(rnorm(t0*s,0,2), t0, s)
      y <- z%*%b.re + xb + err
    }
  }else{
    
    ### S/N Ratio
    
    factor = 10 ; factor_min = 0.01; factor_max = 1000
    count = 0 ; maxit = 10000
    
    factor_prev = 1
    
    repeat{
      
      var = var / factor * factor_prev
      
      ### Sample the errors and the Ys
      cVar = chol(as.matrix(var))
      #err = matrix(rnorm(n*s),n,s) %*% cVar
      err = matrix(rnorm(n*s,sd=sd_err),n,s) %*% cVar
      if(t0 == 0){
        b.re <- 0
        z <- 0
        y <- xb + err
      }else{
        # add random effects
        z <- model.matrix( ~ factor(rMultinom(matrix(c(.1,.2,.3,.4),nrow=1),n)) - 1)
        # use fixed random seed for b.re
        b.re <- matrix(rnorm(t0*s,0.5,sd_b), t0, s)
        zb.re <- z%*%b.re
        y <- zb.re + xb + err
      }
      
      ## Reparametrisation ( assuming PEO is 1:s )
      cVar = t(cVar) # make it lower-tri
      S = diag(diag(cVar))
      sigma = S*S
      L = cVar %*% solve(S)
      rho =  diag(s) - solve(L)
      
      ### S/N Ratio
      if(t0 ==0){
        emp_snr = mean( diag( var(xb) %*% solve(sigma) ))
      }else{
        emp_snr = mean( diag( var(zb.re+xb) %*% solve(sigma) ))
      }
      
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
  }

  return( list(y=y, x=x, b=b, gamma=gamma, z=z, b.re=b.re, Gy=G, mrfG=combnAll ))#, SNR=c(emp_snr, emp_g_snr)) )
  
}
