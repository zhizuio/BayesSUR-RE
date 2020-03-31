#######################
## R-code for all results in paper 
## "Structured Bayesian variable selection for multiple related response variables and high-dimensional predictors"
## 
## author: Zhi Zhao (zhi.zhao@medisin.uio.no)
## date: 31-Mar-2020
#######################


library(Matrix)
library(coda)
library(MCMCpack)
library(gRbase)
library(BDgraph)
library(Formula)
library(Hmisc)
library(BayesSUR)

#################################
#################################
##  Simulation studies  #########
#################################
#################################

#################################
## Simulation Algorithm 1
#################################

# Simulate data without random effects
n=250;s=20;p=300;t0=0
source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv4=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv4=TRUE)

# fit a SSUR model with hyper-inverse Wishart and hotspot priors
hyperpar=list(a_o=2, b_o=298, a_pi=2, b_pi=1, a_tau=0.1, b_tau=10)
fit1.hotspot <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_hotspot/", hyperpar = hyperpar, 
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "hotspot", maxThreads = 16)
summary(fit1.hotspot)

# plot posterior mean of gamma and Gy
plotEstimator(fit1.hotspot, estimator = c("gamma","Gy"), fig.tex = T)

# compute accuracy of variable selection and prediction performance
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# MSE & MSPE
gamma <- getEstimator(fit1.hotspot, "gamma")
beta <- getEstimator(fit1.hotspot, "beta")
beta <- (gamma>0.5)*beta/gamma # MPM estimator
beta[is.na(beta)] <- 0
MSE <- sum((scale(sim1$y)[,]- scale(sim1$x)[,-1]%*%beta)^2)/prod(dim(beta))
MSPE <- sum((scale(sim1.val$y)[,]- scale(sim1.val$x)[,-1]%*%beta)^2)/prod(dim(beta))

# run a SSUR model with MRF prior
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit1.mrf)

# plot posterior mean of gamma and Gy
plotEstimator(fit1.mrf, estimator = c("gamma","Gy"), fig.tex = T)

# compute accuracy of variable selection and prediction performance
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# MSE & MSPE
gamma <- getEstimator(fit1.mrf, "gamma")
beta <- getEstimator(fit1.mrf, "beta")
beta <- (gamma>0.5)*beta/gamma # MPM estimator
beta[is.na(beta)] <- 0
MSE <- sum((scale(sim1$y)[,]- scale(sim1$x)[,-1]%*%beta)^2)/prod(dim(beta))
MSPE <- sum((scale(sim1.val$y)[,]- scale(sim1.val$x)[,-1]%*%beta)^2)/prod(dim(beta))
summary(fit)

## Sensitivity analysis: case 1

# delete 1% edges uniformly
sim1 <- sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))*100*c(0:(0.01*nrow(sim1$mrfG))))),]
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_delete1/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(a) delete 1\\% edges uniformly"), output="delete1",cex.main=1.2)
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.gamma = paste("(a) delete 1\\% edges uniformly"), output="delete1G",cex.main=1.2)

## delete 10% edges uniformly
sim1 <- sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))/0.1*c(0:(0.1*nrow(sim1$mrfG))))),]
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_delete10/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(b) delete 10\\% edges uniformly"), output="delete10",cex.main=1.2)
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.gamma = paste("(b) delete 10\\% edges uniformly"), output="delete10G",cex.main=1.2)

## delete 50% edges uniformly
sim1 <- sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))/0.5*c(0:(0.5*nrow(sim1$mrfG))))),]
hyperpar=list(mrf_d=-2, mrf_e=0.2, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_delete50-2/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(c) delete 50\\% edges uniformly"), output="delete50",cex.main=1.2)
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.gamma = paste("(c) delete 50\\% edges uniformly"), output="delete50G",cex.main=1.2)

## delete 90% edges uniformly
sim1 <- sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))/0.9*c(0:(0.9*nrow(sim1$mrfG))))),]
hyperpar=list(mrf_d=-1, mrf_e=0.5, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_delete90/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(d) delete 90\\% edges uniformly"), output="delete90",cex.main=1.2)
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.gamma = paste("(d) delete 90\\% edges uniformly"), output="delete90G",cex.main=1.2)

## Sensitivity analysis: case 2

## delete last 1%
sim1 <- sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.99)+1):nrow(sim1$mrfG)),]
library("BayesSUR", lib.loc="RLibrary")
# fit a SSUR model with MRF prior
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_deleteL1/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
save(fit, file="sim1_mrf_deleteL1/fit.rda")

## delete last 10%
sim1 = sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.90)+1):nrow(sim1$mrfG)),]
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_deleteL10/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(b) delete 10\\% block edges"), output="deleteL10")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(b) delete 10\\% block edges"), output="deleteL10G")

## delete last 50%
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv4=TRUE)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.50)+1):nrow(sim1$mrfG)),]
hyperpar=list(mrf_d=-1, mrf_e=0.5, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_deleteL50/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(c) delete 50\\% block edges"), output="deleteL50")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(c) delete 50\\% block edges"), output="deleteL50G")

## delete last 90%
sim1 = sim.ssur(n, s, p, t0, seed=7193)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.10)+1):nrow(sim1$mrfG)),]
hyperpar=list(mrf_d=-2, mrf_e=2, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_deleteL90/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(d) delete 90\\% block edges"), output="deleteL90")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(d) delete 90\\% block edges"), output="deleteL90G")

## delete last 100%
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv4=TRUE)
library("BayesSUR", lib.loc="RLibrary")
hyperpar=list(mrf_d=-0.5, mrf_e=1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_deleteL100/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(d) delete 100\\% block edges"), output="deleteL100")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(d) delete 100\\% block edges"), output="deleteL100G")

## Sensitivity analysis: case 3

# add noise 0.1%
sim1 = sim.ssur(n, s, p, t0, seed=7193)
set.seed(2610)
combn.all <- data.frame( t(combn(rep((1:s)*p,each=p) + rep(1:p,times=s), 2)) )
combn.sub <- data.frame( combn.all[sample(1:nrow(combn.all), 20*300*(20*300-1)/2*0.002),] )
mrfG.noise <- setdiff(combn.sub, data.frame(sim1$mrfG))[20*300*(20*300-1)/2*0.001,]
sim1$mrfG <- rbind(sim1$mrfG, data.matrix(mrfG.noise)[,])
colnames(sim1$mrfG) <- NULL
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_add0.1/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(e) delete 100\\% block edges"), output="deleteL100")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(e) delete 100\\% block edges"), output="deleteL100G")

## add noise 0.5%
sim1 = sim.ssur(n, s, p, t0, seed=7193)
set.seed(2610)
combn.all <- data.frame( t(combn(rep((1:s)*p,each=p) + rep(1:p,times=s), 2)) )
combn.sub <- data.frame( combn.all[sample(1:nrow(combn.all), 20*300*(20*300-1)/2*0.006),] )
mrfG.noise <- setdiff(combn.sub, data.frame(sim1$mrfG))[20*300*(20*300-1)/2*0.005,]
sim1$mrfG <- rbind(sim1$mrfG, data.matrix(mrfG.noise)[,])
colnames(sim1$mrfG) <- NULL
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_add0.5/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(b) add 0.5\\% noise edges"), output="add0.5")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(b) add 0.5\\% noise edges"), output="add0.5G")

## add noise 1%
sim1 = sim.ssur(n, s, p, t0, seed=7193)
set.seed(2610)
combn.all <- data.frame( t(combn(rep((1:s)*p,each=p) + rep(1:p,times=s), 2)) )
combn.sub <- data.frame( combn.all[sample(1:nrow(combn.all), 20*300*(20*300-1)/2*0.011),] )
mrfG.noise <- setdiff(combn.sub, data.frame(sim1$mrfG))[20*300*(20*300-1)/2*0.01,]
sim1$mrfG <- rbind(sim1$mrfG, data.matrix(mrfG.noise)[,])
colnames(sim1$mrfG) <- NULL
hyperpar=list(mrf_d=-2, mrf_e=0.1, b_w=3)
fit <- BayesSUR(data = cbind(sim1$y,sim1$x[,-1]), Y = 1:s,
                X = s+1:p, outFilePath = "sim1_mrf_add1/", hyperpar = hyperpar,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 10)
summary(fit)
plotEstimator(fit,estimator = "gamma",fig.tex = T,title.gamma = paste("(c) add 1\\% noise edges"), output="add1")
plotEstimator(fit,estimator = "Gy",fig.tex = T,title.Gy = paste("(c) add 1\\% noise edges"), output="add1G")


#################################
## Simulation Algorithm 2
#################################

# Simulate data with random effects
n=250;s=20;p=300;t0=4
sim2 = sim.ssur.re(n, s, p, t0, seed=7193)

# fit a SSUR model with hotspot prior
hyperpar=list(mrf_d=-2, mrf_e=0.1, a_w0=60, b_w0=30, a_w=300, b_w=300)
fit2 <- BayesSUR(data = cbind(sim2$y,sim2$z,sim2$x), Y = 1:s, X_0 = s+1:t0, hyperpar = hyperpar,
                X = s+t0+1:p, outFilePath = "sim2_mrf/", betaPrior = "reGroup", standardize=F, standardize.response=F,
                nIter = 500000, nChains = 5, covariancePrior = "HIW", burnin=300000,
                gammaPrior = "MRF",  mrfG = sim2$mrfG, maxThreads = 16)
summary(fit2)
plotEstimator(fit, estimator = c("gamma","Gy"), fig.tex=T, ylab="Responses", xlab="Predictors")

# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit2, "gamma")
beta <- getEstimator(fit2, "beta")
gamma.re <- rbind(matrix(1, nrow=nrow(beta)-nrow(gamma),ncol=ncol(beta)), gamma)
beta.MPM <- (gamma.re>0.5)*beta/gamma.re
beta.MPM[is.na(beta.MPM)] <- 0
g.re <- getEstimator(fit, "Gy")
g.accuracy <- sum((g.re>0.5)==sim1$Gy)/prod(dim(g.re))
g.sensitivity <- sum(((g.re>0.5)==sim1$Gy)[sim1$Gy==1])/sum(sim1$Gy==1)
g.specificity <- sum(((g.re>0.5)==sim1$Gy)[sim1$Gy==0])/sum(sim1$Gy==0)

beta.accuracy <- sum((gamma>0.5)==sim1$gamma)/prod(dim(gamma))
beta.sensitivity <- sum(((gamma>0.5)==sim1$gamma)[sim1$gamma==1])/sum(sim1$gamma==1)
beta.specificity <- sum(((gamma>0.5)==sim1$gamma)[sim1$gamma==0])/sum(sim1$gamma==0)

beta.l2 <- sum((beta.MPM[-c(1:4),]-sim1$b)^2)/prod(dim(sim1$b))
beta0.l2 <- sum((beta.MPM[1:4,]-sim1$b.re)^2)/prod(dim(sim1$b.re))

#################################
#################################
##  GDSC data analysis  #########
#################################
#################################
rm(list = ls())
library(BayesSUR)
# load the GDSC data set that is generated from the file GDSC_mrfG2.R
load("example_GDSC.rda")
example_GDSC$data[,7+13+1:343] <- log(2^(example_GDSC$data[,7+13+1:343]))
hyperpar=list(mrf_d=-2.7, mrf_e=0.1, a_w0=54.6, b_w0=400, a_w=11.718, b_w=70.308, a_eta=0.1)
example_GDSC$data <- cbind(example_GDSC$data[,1:7], rep(1,nrow(example_GDSC$data)), example_GDSC$data[,8:857])
fit <- BayesSUR(data = example_GDSC$data, Y = example_GDSC$blockList[[1]], X_0 = c(8,1+example_GDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+example_GDSC$blockList[[3]], outFilePath = "GDSC/", hyperpar=hyperpar,
                nIter = 800000, nChains = 25, covariancePrior = "HIW", burnin=500000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = example_GDSC$mrfG, maxThreads = 16)
summary(fit)

# plot estimators
plotEstimator(fit, fig.tex=T, xlab="Predictors", ylab="Drugs")

# plot network representation
pdf("GDSC_network.pdf",width=10, height=10)
plotNetwork(fit, lineup=1.5,label.predictor = "", name.predictors = "Genes", name.responses = "Drugs", nodesizePredictor = 2)
dev.off()

# MCMC diagnostics
pdf("GDSC_MCMCdiag.pdf",width=10, height=10)
plotMCMCdiag(fit)
dev.off()

# plot barplot for estimated tissue effects
name.responses <- names(read.table(paste(fit$output$outFilePath,fit$output$Y,sep=""),header=T))
name.genes <- names(read.table(paste(fit$output$outFilePath,fit$output$X,sep=""),header=T))
name.tissues <- c("digestive","urogenital","blood","kidney","nervous","thyroid","skin","soft","aerodigestive","lung","pancreas","breast","bone")
name.predictors <- c("intercept", name.tissues, name.genes)
beta <- getEstimator(fit, "beta")
gamma <- getEstimator(fit)
rownames(gamma) <- name.genes
rownames(gamma)[rowSums(gamma>0.5)==7]
gamma <- rbind(matrix(1,nrow=nrow(beta)-nrow(gamma),ncol=ncol(beta)),gamma)
beta.MPM <- (gamma>0.5)*beta/gamma
beta.MPM[is.na(beta.MPM)] <- 0
colnames(beta.MPM) <- name.responses
rownames(beta.MPM) <- name.predictors
pdf("TissueEffects.pdf",width=8,height=15)
par(mar=c(5+2, 4, 4, 2) + 0.1)
layout(matrix(c(1,0,2:7),ncol=2,byrow=T))
for(i in 1:7){
  ext.idx <- 2:14
  barplot(beta.MPM[ext.idx,i], ylim=c(-6,6),ylab=name.responses[i], las=2, border=NA, cex.names=1)
  axis(side=2,at=(-3:3),labels=(-3:3), las=2)
  abline(h=0, lty=3)
}
title("\nTissue effects (MPM)", outer=T)
dev.off()

