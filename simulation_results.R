#######################
## This script is for all simulation results in paper 
## "Structured multivariate Bayesian variable selection for pharmacogenomic studies"
##
## Author: Zhi Zhao (zhi.zhao@medisin.uio.no)
## Date: 01-June-2022
#######################

#################################
## Simulation Algorithm 1
#################################

# Simulate data without random effects
rm(list = ls())
library(Matrix)

n=250;s=20;p=300;t0=0
source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)

####
## fit a SSUR-hotspot model
####
set.seed(1038)
hyperpar=list(a_w=15, b_w=60, a_o=100, a_pi=50)
fit1.hotspot <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                         X = s+1+1:p, outFilePath = "sim1_hotspot", hyperpar = hyperpar, gammaInit="0",
                         nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                         gammaPrior = "hotspot",  maxThreads = 1, output_CPO = T)
save(fit1.hotspot, file="sim1_hotspot/fit.rda")

# Show results for Figure 1(a) and 1(b)
plot(fit1.hotspot, estimator=c("gamma","Gy"), type="heatmap",fig.tex = TRUE, output="sim1_hotspot",xlab = "Predictors", ylab = "Responses")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 1
summary(fit1.hotspot) 

## Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 1
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.hotspot)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.hotspot, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

####
## run a SSUR-MRF model
####

rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf/", hyperpar = hyperpar, gammaInit="0",
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, output_CPO = T)
save(fit1.mrf, file="sim1_mrf/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 1
summary(fit1.mrf) 

# Show results for Figure 1(c) and 1(d)
#plot(fit1.mrf, estimator=c("gamma","Gy"), type="heatmap") # Run this instead of the following line if not need LaTeX math symbols
plot(fit1.mrf, estimator=c("gamma","Gy"), type="heatmap", fig.tex=T, output="sim1_mrf", ylab="Responses", xlab="Predictors")

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 1
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

######
# Sensitivity analysis: case 1
######

# delete 1% edges uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))*100*c(0:(0.01*nrow(sim1$mrfG))))),]
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1.4, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_delete1/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_delete1/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

## delete 10% edges uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))/0.1*c(0:(0.1*nrow(sim1$mrfG))))),]
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1.4, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_delete10/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_delete10/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

## delete 50% edges uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))/0.5*c(0:(0.5*nrow(sim1$mrfG))))),]
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=0.8, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_delete50/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_delete50/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

## delete 90% edges uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c(1+floor((1-1/nrow(sim1$mrfG))/0.9*c(0:(0.9*nrow(sim1$mrfG))))),]
set.seed(1038)
hyperpar=list(mrf_d=-0.8, mrf_e=0.4, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_delete90/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_delete90/fit.rda") 

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

#######
## Sensitivity analysis: case 2
#######

## delete 1% edges non-uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.99)+1):nrow(sim1$mrfG)),]
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=0.8, a_w=15, b_w=60) 
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_deleteL1/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_deleteL1/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

# Show Figure 4(a)
plotEstimator(fit1.mrf,estimator = "gamma",fig.tex = T,title.gamma = paste("(a) delete 1\\% block edges"), output="deleteL1")

## delete 10% edges non-uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.90)+1):nrow(sim1$mrfG)),]
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1, a_w=15, b_w=60) 
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_deleteL10/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_deleteL10/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

# Show Figure 4(b)
plotEstimator(fit1.mrf,estimator = "gamma",fig.tex = T,title.gamma = paste("(b) delete 10\\% block edges"), output="deleteL10")

## delete 50% edges non-uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.50)+1):nrow(sim1$mrfG)),]
set.seed(1038)
hyperpar=list(mrf_d=-0.8, mrf_e=0.8, a_w=15, b_w=60) 
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_deleteL50/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_deleteL50/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

# Show Figure 4(c)
plotEstimator(fit1.mrf,estimator = "gamma",fig.tex = T,title.gamma = paste("(c) delete 50\\% block edges"), output="deleteL50")

## delete 90% edges non-uniformly
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- sim1$mrfG[-c((floor(nrow(sim1$mrfG)*0.10)+1):nrow(sim1$mrfG)),]
set.seed(1038)
hyperpar=list(mrf_d=-0.4, mrf_e=0.8, a_w=15, b_w=60) 
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_deleteL90/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_deleteL90/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

# Show Figure 4(d)
plotEstimator(fit1.mrf,estimator = "gamma",fig.tex = T,title.gamma = paste("(d) delete 90\\% block edges"), output="deleteL90")


########
## delete 100% edges (uniformly)
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
sim1$mrfG <- matrix(1, ncol=2, nrow=1)
set.seed(1038)
hyperpar=list(mrf_d=-0.1, mrf_e=0, a_w=15, b_w=60) 
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_delete100/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_delete100/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

# Show Figure 4(e)
plotEstimator(fit1.mrf,estimator = "gamma",fig.tex = T,title.gamma = paste("(e) delete 100\\% block edges"), output="deleteL100")


#############
## Sensitivity analysis: case 3
###########
# add noise 0.1%
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
set.seed(2610)
combn.all <- data.frame( t(combn(rep((1:s)*p,each=p) + rep(1:p,times=s), 2)) )
combn.sub <- data.frame( combn.all[sample(1:nrow(combn.all), 20*300*(20*300-1)/2*0.002),] )
mrfG.noise <- setdiff(combn.sub, data.frame(sim1$mrfG))[20*300*(20*300-1)/2*0.001,]
sim1$mrfG <- rbind(sim1$mrfG, data.matrix(mrfG.noise)[,])
colnames(sim1$mrfG) <- NULL
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=0.8, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_add0.1/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_add0.1/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

## add noise 0.5%
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
set.seed(2610)
combn.all <- data.frame( t(combn(rep((1:s)*p,each=p) + rep(1:p,times=s), 2)) )
combn.sub <- data.frame( combn.all[sample(1:nrow(combn.all), 20*300*(20*300-1)/2*0.006),] )
mrfG.noise <- setdiff(combn.sub, data.frame(sim1$mrfG))[20*300*(20*300-1)/2*0.005,]
sim1$mrfG <- rbind(sim1$mrfG, data.matrix(mrfG.noise)[,])
colnames(sim1$mrfG) <- NULL
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1.6, a_w=15, b_w=60) 
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_add0.5/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_add0.5/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

## add noise 1%
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
set.seed(2610)
combn.all <- data.frame( t(combn(rep((1:s)*p,each=p) + rep(1:p,times=s), 2)) )
combn.sub <- data.frame( combn.all[sample(1:nrow(combn.all), 20*300*(20*300-1)/2*0.011),] )
mrfG.noise <- setdiff(combn.sub, data.frame(sim1$mrfG))[20*300*(20*300-1)/2*0.01,]
sim1$mrfG <- rbind(sim1$mrfG, data.matrix(mrfG.noise)[,])
colnames(sim1$mrfG) <- NULL
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=0.8, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_add1/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, maxThreads = 1, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_add1/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE

#######
## use Kroneker product to construct edge potentials
####
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R") 
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
# use Kroneker product to construct edge potentials
Idx_Y <- list(g1=1:5, g2=6:12, g3=13:20)
Idx_X <- list(g1=1:5, g2=10:20, g3=30:50, g4=51:60, g5=70:90, g6=110:120)
mrfG <- rbind(t(combn(rep(Idx_X[[1]],each=length(Idx_Y[[1]])) + rep((Idx_Y[[1]]-1)*sum(p),times=length(Idx_X[[1]])), 2)),
              t(combn(rep(Idx_X[[2]], each = length(Idx_Y[[2]])) + rep((Idx_Y[[2]] - 1) * sum(p), times = length(Idx_X[[2]])), 2)),
              t(combn(rep(Idx_X[[3]], each = length(Idx_Y[[1]])) + rep((Idx_Y[[1]] - 1) * sum(p), times = length(Idx_X[[3]])), 2)),
              t(combn(rep(Idx_X[[3]], each = length(Idx_Y[[3]])) + rep((Idx_Y[[3]] - 1) * sum(p), times = length(Idx_X[[3]])), 2)),
              t(combn(rep(Idx_X[[4]], each = length(Idx_Y[[1]])) + rep((Idx_Y[[1]] - 1) * sum(p), times = length(Idx_X[[4]])), 2)),
              t(combn(rep(Idx_X[[4]], each = length(Idx_Y[[2]])) + rep((Idx_Y[[2]] - 1) * sum(p), times = length(Idx_X[[4]])), 2)),
              t(combn(rep(Idx_X[[5]], each = length(Idx_Y[[2]])) + rep((Idx_Y[[2]] - 1) * sum(p), times = length(Idx_X[[5]])), 2)),
              t(combn(rep(Idx_X[[5]], each = length(Idx_Y[[3]])) + rep((Idx_Y[[3]] - 1) * sum(p), times = length(Idx_X[[5]])), 2)),
              t(combn(rep(Idx_X[[6]], each = length(Idx_Y[[1]])) + rep((Idx_Y[[1]] - 1) * sum(p), times = length(Idx_X[[6]])), 2)),
              t(combn(rep(Idx_X[[6]], each = length(Idx_Y[[2]])) + rep((Idx_Y[[2]] - 1) * sum(p), times = length(Idx_X[[6]])), 2)),
              t(combn(rep(Idx_X[[6]], each = length(Idx_Y[[3]])) + rep((Idx_Y[[3]] - 1) * sum(p), times = length(Idx_X[[6]])), 2)))

set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=0.2, a_w=15, b_w=60)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "sim1_mrf_kroneker5/", hyperpar = hyperpar, 
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = mrfG, output_CPO = T)
save(fit1.mrf, file="sim1_mrf_kroneker/fit.rda")
# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 2
summary(fit1.mrf)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE


#################################
## Simulation Algorithm 2
#################################

# Simulate heterogeneous data
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=4; source("sim.ssur.R")
sim2 = sim.ssur(n, s, p, t0, seed=7193)
sim2.val = sim.ssur(n, s, p, t0, seed=822)

# fit a SSUR-MRF without random effects, but with intercept
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=10, a_w=15, b_w=60, a_eta=0.001, b_eta=1000) 
fit2 <- BayesSUR(data = cbind(sim2$y,rep(1,n),sim2$x), Y = 1:s, X_0 = s+1, 
                 X = s+1+1:p, outFilePath = "sim2_mrf1/", hyperpar = hyperpar, 
                 nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                 gammaPrior = "MRF",  mrfG = sim2$mrfG, maxThreads = 1, output_CPO = T)
save(fit2, file="sim2_mrf1/fit.rda")

# Show results for Figure 5(c) and 5(d)
plot(fit2, estimator = c("gamma","Gy"), type="heatmap", fig.tex = T, output="mrfRE1", ylab="Responses", xlab="Predictors")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 3
summary(fit2)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit2)
(accuracy <- sum(data.matrix(gamma>0.5) == sim2$gamma)/prod(dim(gamma)))
(sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim2$gamma==1))/sum(sim2$gamma==1))
(specificity <- sum((data.matrix(gamma>0.5)==0) & (sim2$gamma==0))/sum(sim2$gamma==0))
# compute RMSE & RMSPE
beta <- getEstimator(fit2, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
(MSE <- sqrt(sum((sim2$y - cbind(rep(1,n),sim2$x)%*%beta)^2)/prod(dim(sim2$y))))
(RMSPE <- sqrt(sum((sim2.val$y - cbind(rep(1,n),sim2.val$x)%*%beta)^2)/prod(dim(sim2.val$y))))
# compute bias of beta estimates
b <- sim2$b; b[sim2$gamma==0] <- 0
(beta.l2 <- sqrt(sum((beta[-1,]-b)^2)/prod(dim(b))))

g.re <- getEstimator(fit2, estimator="Gy")
(g.accuracy <- sum((g.re>0.5)==sim2$Gy)/prod(dim(g.re)))
(g.sensitivity <- sum(((g.re>0.5)==sim2$Gy)[sim2$Gy==1])/sum(sim2$Gy==1))
(g.specificity <- sum(((g.re>0.5)==sim2$Gy)[sim2$Gy==0])/sum(sim2$Gy==0))


# fit a SSUR-MRF with random effects
rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=4; source("sim.ssur.R")
sim2 = sim.ssur(n, s, p, t0, seed=7193)
sim2.val = sim.ssur(n, s, p, t0, seed=822)
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1.2, a_w0=100, b_w0=500, a_w=15, b_w=60)
fit2 <- BayesSUR(data = cbind(sim2$y,sim2$z,sim2$x), Y = 1:s, X_0 = s+1:t0,
                 X = s+t0+1:p, outFilePath = "sim2_mrf2/", hyperpar = hyperpar, betaPrior = "reGroup", 
                 nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                 gammaPrior = "MRF",  mrfG = sim2$mrfG, maxThreads = 1, output_CPO = T)
save(fit2, file="sim2_mrf2/fit.rda")

# Show results for Figure 5(a) and 5(b)
plot(fit2, estimator = c("gamma","Gy"), type="heatmap", fig.tex = T, output="mrfRE2", ylab="Responses", xlab="Predictors")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 3
summary(fit2)

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 2
# compute accuracy of variable selection and prediction performance
gamma <- getEstimator(fit2)
(accuracy <- sum(data.matrix(gamma>0.5) == sim2$gamma)/prod(dim(gamma)))
(sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim2$gamma==1))/sum(sim2$gamma==1))
(specificity <- sum((data.matrix(gamma>0.5)==0) & (sim2$gamma==0))/sum(sim2$gamma==0))
# compute RMSE & RMSPE
beta <- getEstimator(fit2, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
(MSE <- sqrt(sum((sim2$y - cbind(sim2$z,sim2$x)%*%beta)^2)/prod(dim(sim2$y))))
(RMSPE <- sqrt(sum((sim2.val$y - cbind(sim2.val$z,sim2.val$x)%*%beta)^2)/prod(dim(sim2.val$y))))
# compute bias of beta estimates
b <- sim2$b; b[sim2$gamma==0] <- 0
beta0.l1 <- sqrt(sum((beta[1:4,]-sim2$b.re)^2)/prod(dim(sim2$b.re)))

g.re <- getEstimator(fit2, estimator="Gy")
(g.accuracy <- sum((g.re>0.5)==sim2$Gy)/prod(dim(g.re)))
(g.sensitivity <- sum(((g.re>0.5)==sim2$Gy)[sim2$Gy==1])/sum(sim2$Gy==1))
(g.specificity <- sum(((g.re>0.5)==sim2$Gy)[sim2$Gy==0])/sum(sim2$Gy==0))

####
## run supplementary simulations for additional sensitivity analysis
####

rm(list = ls())
library(BayesSUR)
library(Matrix)
n=250;s=20;p=300;t0=0; source("sim.ssur.R")
sim1 = sim.ssur(n, s, p, t0, seed=7193, mv=TRUE)
sim1.val = sim.ssur(n, s, p, t0, seed=822, mv=TRUE)
set.seed(1038)
hyperpar=list(mrf_d=-2, mrf_e=1, a_w=15, b_w=60, nu=100)
hyperpar=list(mrf_d=-2, mrf_e=1, a_w=15, b_w=60, a_tau=0.001)
hyperpar=list(mrf_d=-2, mrf_e=1, a_w=15, b_w=60, a_tau=100)
hyperpar=list(mrf_d=-2, mrf_e=1, a_w=15, b_w=60, a_eta=0.001)
fit1.mrf <- BayesSUR(data = cbind(sim1$y,sim1$x), Y = 1:s, X_0 = s+1,
                     X = s+1+1:p, outFilePath = "/data/zhiz/bayessur_sim/sim1_mrf_suppl4/", hyperpar = hyperpar, gammaInit="0",
                     nIter = 400000, nChains = 5, covariancePrior = "HIW", burnin=200000, standardize=F, standardize.response=F,
                     gammaPrior = "MRF",  mrfG = sim1$mrfG, output_CPO = T)
save(fit1.mrf, file="/data/zhiz/bayessur_sim/sim1_mrf_suppl4/fit.rda")

# Show results, elpd.loo and elpd.waic values from 'summary()' for Table 1
summary(fit1.mrf) 

# Show results for Figure 1(c) and 1(d)
#plot(fit1.mrf, estimator=c("gamma","Gy"), type="heatmap") # Run this instead of the following line if not need LaTeX math symbols
plot(fit1.mrf, estimator=c("gamma","Gy"), type="heatmap", fig.tex=T, output="sim1_mrf_suppl4", ylab="Responses", xlab="Predictors")

# Show results accuracy, sensitivity, specificity, RMSE and RMSPE for Table 1
# compute accuracy of variable selection and prediction performance
summary(fit1.mrf) 
gamma <- getEstimator(fit1.mrf)
accuracy <- sum(data.matrix(gamma>0.5) == sim1$gamma[-1,])/prod(dim(gamma))
sensitivity <- sum((data.matrix(gamma>0.5)==1) & (sim1$gamma[-1,]==1))/sum(sim1$gamma[-1,]==1)
specificity <- sum((data.matrix(gamma>0.5)==0) & (sim1$gamma[-1,]==0))/sum(sim1$gamma[-1,]==0)
# compute RMSE & RMSPE
beta <- getEstimator(fit1.mrf, estimator = "beta", Pmax=.5, beta.type = "conditional") # MPM estimator
RMSE <- sqrt(sum((sim1$y - sim1$x%*%beta)^2)/prod(dim(sim1$y)))
RMSPE <- sqrt(sum((sim1.val$y - sim1.val$x%*%beta)^2)/prod(dim(sim1.val$y)))
accuracy; sensitivity; specificity;
RMSE; RMSPE
