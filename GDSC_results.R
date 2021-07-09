#######################
## This script is for all results of real data analysis in paper 
## "Structured multivariate Bayesian variable selection for pharmacogenomic studies"
## 
## author: Zhi Zhao (zhi.zhao@medisin.uio.no)
## date: 08-July-2021
#######################

rm(list = ls())
library(BayesSUR)

#####
## SSUR-MRF for Feature set I 
#####

rm(list = ls())
set.seed(18372)
library("BayesSUR")
load("example_GDSC_kegg10%.RData")
load("targetGene_kegg_10%.RData")
exampleGDSC <- example_GDSC_kegg
exampleGDSC$data <- cbind(exampleGDSC$data[,1:7], rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)])
hyperpar <- list(mrf_d=-2.7, mrf_e=0.2, a_w0=55, b_w0=400, a_w=4, b_w=32)
fit <- BayesSUR(data = exampleGDSC$data, Y = exampleGDSC$blockList[[1]], X_0 = c(8,1+exampleGDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+exampleGDSC$blockList[[3]], outFilePath = "GDSC_kegg10/", hyperpar=hyperpar,
                nIter = 500000, nChains = 6, covariancePrior = "HIW", burnin=300000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = exampleGDSC$mrfG, maxThreads = 1, output_CPO = T)
summary(fit)

# Show Table 4 line SSUR-MRF: Feature set I
gamma <- getEstimator(fit); colSums(gamma>.5) 

# Show Table 5 line SSUR-MRF: Feature set I: #identified features
gene.names <- colnames(exampleGDSC$data)[-c(1:20)];  length(gene.names[rowSums(gamma[,2:5]>0.5)>0])
# Show Table 5 line SSUR-MRF: Feature set I: #identified targets
gene.sel10 <- gene.names[rowSums(gamma[,2:5]>0.5)>0]; length(gene.sel10[gene.sel10 %in% colnames(targetGene[[1]])])
save(gene.sel10, file="gene.sel10.rda")

# Show Figure 9(a) 
plotNetwork(fit, estimator = c("gamma","Gy"),includeResponse = c("RDEA119", "PD.0325901", "CI.1040", "AZD6244"),
            includePredictor = colnames(targetGene[[1]]), main="(a)")

#####
## SSUR-MRF for Feature set II 
#####
rm(list = ls())
set.seed(18372)
library("BayesSUR")
load("example_GDSC_kegg30%.RData")
load("targetGene_kegg_30%.RData")
exampleGDSC <- example_GDSC_kegg
exampleGDSC$data <- cbind(exampleGDSC$data[,1:7], rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)])
hyperpar <- list(mrf_d=-4, mrf_e=0.3, a_w0=55, b_w0=400, a_w=4, b_w=32)
fit <- BayesSUR(data = exampleGDSC$data, Y = exampleGDSC$blockList[[1]], X_0 = c(8,1+exampleGDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+exampleGDSC$blockList[[3]], outFilePath = "GDSC_kegg30/", hyperpar=hyperpar,
                nIter = 500000, nChains = 8, covariancePrior = "HIW", burnin=300000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = exampleGDSC$mrfG, maxThreads = 1, output_CPO = T)

# Show Table 4 line SSUR-MRF: Feature set II
gamma <- getEstimator(fit); colSums(gamma>.5) 

# Show Table 5 line SSUR-MRF: Feature set II: #identified features
gene.names <- colnames(exampleGDSC$data)[-c(1:20)];  length(gene.names[rowSums(gamma[,2:5]>0.5)>0])
# Show Table 5 line SSUR-MRF: Feature set II: #identified targets
gene.sel30 <- gene.names[rowSums(gamma[,2:5]>0.5)>0]; length(gene.sel30[gene.sel30 %in% colnames(targetGene[[1]])])
save(gene.sel30, file="gene.sel30.rda")

# Show Figure 9(b)
plotNetwork(fit, estimator = c("gamma","Gy"),includeResponse = c("RDEA119", "PD.0325901", "CI.1040", "AZD6244"),
            includePredictor = colnames(targetGene[[1]]), main="(b)")

#####
# SSUR-MRF for Feature set III 
#####
rm(list = ls())
set.seed(18372)
library("BayesSUR")
load("example_GDSC_kegg50%.RData")
load("targetGene_kegg_50%.RData")
exampleGDSC <- example_GDSC_kegg
exampleGDSC$data <- cbind(exampleGDSC$data[,1:7], rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)])
hyperpar <- list(mrf_d=-4.6, mrf_e=0.5, a_w0=55, b_w0=400, a_w=4, b_w=33)
fit <- BayesSUR(data = exampleGDSC$data, Y = exampleGDSC$blockList[[1]], X_0 = c(8,1+exampleGDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+exampleGDSC$blockList[[3]], outFilePath = "GDSC_kegg50/", hyperpar=hyperpar,
                nIter = 500000, nChains = 10, covariancePrior = "HIW", burnin=300000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = exampleGDSC$mrfG, maxThreads = 1, output_CPO = T)
save(fit, file="GDSC_kegg50/fit.rda")

# Show Table 6 SSUR-MRF:elpd.loo and elpd.waic in the paper
summary(fit)

# Show Figure 7 in the paper
plotGraph(fit, estimator = "Gy")

# Show Table 4 line SSUR-MRF: Feature set III
gamma <- getEstimator(fit); colSums(gamma>.5) 

# Show Table 5 line SSUR-MRF: Feature set III: #identified features
gene.names <- colnames(exampleGDSC$data)[-c(1:20)];  length(gene.names[rowSums(gamma[,2:5]>0.5)>0])
# Show Table 5 line SSUR-MRF: Feature set III: #identified targets
gene.sel50 <- gene.names[rowSums(gamma[,2:5]>0.5)>0]; length(gene.sel50[gene.sel50 %in% colnames(targetGene[[1]])])
save(gene.sel50, file="gene.sel50.rda")

# Show Figure 9(c)
plotNetwork(fit, estimator = c("gamma","Gy"),includeResponse = c("RDEA119", "PD.0325901", "CI.1040", "AZD6244"),
            includePredictor = colnames(targetGene[[1]]), main="(c)")

# Show Table 6 SSUR-MRF:RMSE and RMSPE
beta <- getEstimator(fit, estimator = "beta", Pmax=.5, beta.type = "conditional")
rmse <- sqrt(sum((cbind(rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)]) %*% beta - exampleGDSC$data[,1:7])^2)/prod(dim(exampleGDSC$data[,1:7])));rmse
# load the independent validation data which was produced by file 'GDSC_preprocess2.R'
load("example_GDSC_kegg50%_val.RData")
rmspe <- sqrt(sum((cbind(rep(1,nrow(example_GDSC_kegg50.val)), example_GDSC_kegg50.val[,8:ncol(example_GDSC_kegg50.val)]) %*% beta - example_GDSC_kegg50.val[,1:7])^2)/prod(dim(example_GDSC_kegg50.val[,1:7])));rmspe

# Show Figure 10
## barplot for the tissue effects
beta0 <- beta[2:14,]
### error bars for the tissue effects
betaSD <- read.table("GDSC_kegg50/data_SSUR_betaSD_out.txt")[1:13,]
#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

ze_barplot <- barplot(beta0[,1], border=NA, ylim=c(-4,2.5),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,1], betaSD[,1], col="darkgray")
ze_barplot <- barplot(beta0[,2], border=NA, ylim=c(-3.5,3),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,2], betaSD[,1], col="darkgray")
ze_barplot <- barplot(beta0[,3], border=NA, ylim=c(-3.5,3),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,3], betaSD[,1], col="darkgray")
ze_barplot <- barplot(beta0[,4], border=NA, ylim=c(-3.5,3),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,4], betaSD[,1], col="darkgray")
ze_barplot <- barplot(beta0[,5], border=NA, ylim=c(-3.5,3),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,5], betaSD[,1], col="darkgray")
ze_barplot <- barplot(beta0[,6], border=NA, ylim=c(-4,2.5),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,6], betaSD[,1], col="darkgray")
ze_barplot <- barplot(beta0[,7], border=NA, ylim=c(-4,2.5),las=2);box();abline(h=0, lty=3)
error.bar(ze_barplot,beta0[,7], betaSD[,1], col="darkgray")

# Show Figure 8(b)
load("gene.sel10.rda");load("gene.sel30.rda");load("gene.sel50.rda")
library("VennDiagram")
venn.diagram(
  x = list(gene.sel10, gene.sel30, gene.sel50),
  height = 1700, width = 1700,
  category.names = c("Feature set I" , "Feature set II" , "Feature set III"),
  filename = 'venn_diagram1.png',
  output=TRUE)



#####
# SSUR-Ber for Feature set I 
#####
rm(list = ls())
set.seed(18372)
library("BayesSUR")
load("example_GDSC_kegg10%.RData")
load("targetGene_kegg_10%.RData")
exampleGDSC <- example_GDSC_kegg
exampleGDSC$mrfG <- matrix(1, ncol=2, nrow=1)
exampleGDSC$data <- cbind(exampleGDSC$data[,1:7], rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)])
hyperpar <- list(mrf_d=-2, mrf_e=0, a_w0=55, b_w0=400, a_w=4, b_w=32)
fit <- BayesSUR(data = exampleGDSC$data, Y = exampleGDSC$blockList[[1]], X_0 = c(8,1+exampleGDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+exampleGDSC$blockList[[3]], outFilePath = "GDSC_kegg_ber10/", hyperpar=hyperpar,
                nIter = 500000, nChains = 6, covariancePrior = "HIW", burnin=300000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = exampleGDSC$mrfG, maxThreads = 1, output_CPO = T)
summary(fit)

# Show Table 4 line SSUR-Ber: Feature set I
gamma <- getEstimator(fit); colSums(gamma>.5) 

# Show Table 5 line SSUR-Ber: Feature set I: #identified features
gene.names <- colnames(exampleGDSC$data)[-c(1:20)];  length(gene.names[rowSums(gamma[,2:5]>0.5)>0])
# Show Table 5 line SSUR-Ber: Feature set I: #identified targets
gene.sel_ber10 <- gene.names[rowSums(gamma[,2:5]>0.5)>0]; length(gene.sel_ber10[gene.sel_ber10 %in% colnames(targetGene[[1]])])
save(gene.sel_ber10, file="gene.sel_ber10.rda")


#####
# SSUR-Ber for Feature set II
#####
rm(list = ls())
set.seed(18372)
library("BayesSUR")
load("example_GDSC_kegg30%.RData")
load("targetGene_kegg_30%.RData")
exampleGDSC <- example_GDSC_kegg
exampleGDSC$mrfG <- matrix(1, ncol=2, nrow=1)
exampleGDSC$data <- cbind(exampleGDSC$data[,1:7], rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)])
hyperpar <- list(mrf_d=-4, mrf_e=0, a_w0=55, b_w0=400, a_w=4, b_w=32)
fit <- BayesSUR(data = exampleGDSC$data, Y = exampleGDSC$blockList[[1]], X_0 = c(8,1+exampleGDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+exampleGDSC$blockList[[3]], outFilePath = "GDSC_kegg_ber30/", hyperpar=hyperpar,
                nIter = 500000, nChains = 8, covariancePrior = "HIW", burnin=300000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = exampleGDSC$mrfG, maxThreads = 1, output_CPO = T)
summary(fit)

# Show Table 4 line SSUR-Ber: Feature set II
gamma <- getEstimator(fit); colSums(gamma>.5) 

# Show Table 5 line SSUR-Ber: Feature set II: #identified features
gene.names <- colnames(exampleGDSC$data)[-c(1:20)];  length(gene.names[rowSums(gamma[,2:5]>0.5)>0])
# Show Table 5 line SSUR-Ber: Feature set II: #identified targets
gene.sel_ber30 <- gene.names[rowSums(gamma[,2:5]>0.5)>0]; length(gene.sel_ber30[gene.sel_ber30 %in% colnames(targetGene[[1]])])
save(gene.sel_ber30, file="gene.sel_ber30.rda")


#####
# SSUR-Ber for Feature set III
#####
rm(list = ls())
set.seed(18372)
library("BayesSUR")
load("example_GDSC_kegg50%.RData")
load("targetGene_kegg_50%.RData")
exampleGDSC <- example_GDSC_kegg
exampleGDSC$mrfG <- matrix(1, ncol=2, nrow=1)
exampleGDSC$data <- cbind(exampleGDSC$data[,1:7], rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)])
hyperpar <- list(mrf_d=-3.5, mrf_e=0, a_w0=55, b_w0=400, a_w=4, b_w=32)
fit <- BayesSUR(data = exampleGDSC$data, Y = exampleGDSC$blockList[[1]], X_0 = c(8,1+exampleGDSC$blockList[[2]]), betaPrior = "reGroup",
                X = 1+exampleGDSC$blockList[[3]], outFilePath = "GDSC_kegg_ber50/", hyperpar=hyperpar,
                nIter = 500000, nChains = 8, covariancePrior = "HIW", burnin=300000, standardize=F, standardize.response=F,
                gammaPrior = "MRF", mrfG = exampleGDSC$mrfG, maxThreads = 1, output_CPO = T)
summary(fit)

# Show Table 4 line SSUR-Ber: Feature set IIII
gamma <- getEstimator(fit); colSums(gamma>.5) 

# Show Table 5 line SSUR-Ber: Feature set III: #identified features
gene.names <- colnames(exampleGDSC$data)[-c(1:20)];  length(gene.names[rowSums(gamma[,2:5]>0.5)>0])
# Show Table 5 line SSUR-Ber: Feature set III: #identified targets
gene.sel_ber50 <- gene.names[rowSums(gamma[,2:5]>0.5)>0]; length(gene.sel_ber50[gene.sel_ber50 %in% colnames(targetGene[[1]])])
save(gene.sel_ber50, file="gene.sel_ber50.rda")

# Show Table 6 SSUR-Ber:elpd.loo and elpd.waic in the paper
summary(fit)
# Show Table 6 SSUR-Ber:RMSE and RMSPE
beta <- getEstimator(fit, estimator = "beta", Pmax=.5, beta.type = "conditional")
rmse <- sqrt(sum((cbind(rep(1,nrow(exampleGDSC$data)), exampleGDSC$data[,8:ncol(exampleGDSC$data)]) %*% beta - exampleGDSC$data[,1:7])^2)/prod(dim(exampleGDSC$data[,1:7])));rmse
# load the independent validation data which was produced by file 'GDSC_preprocess2.R'
load("example_GDSC_kegg50%_val.RData")
rmspe <- sqrt(sum((cbind(rep(1,nrow(example_GDSC_kegg50.val)), example_GDSC_kegg50.val[,8:ncol(example_GDSC_kegg50.val)]) %*% beta - example_GDSC_kegg50.val[,1:7])^2)/prod(dim(example_GDSC_kegg50.val[,1:7])));rmspe

# Show Figure 8(a)
load("gene.sel_ber10.rda");load("gene.sel_ber30.rda");load("gene.sel_ber50.rda")
library("VennDiagram")
venn.diagram(
  x = list(gene.sel_ber10, gene.sel_ber30, gene.sel_ber50),
  height = 1700, width = 1700,
  category.names = c("Feature set I" , "Feature set II" , "Feature set III"),
  filename = 'venn_diagram2.png',
  output=TRUE)





