#######################
## This script is for producing independent validation GDSC data 
## "Structured multivariate Bayesian variable selection for pharmacogenomic studies"
## 
## author: Zhi Zhao (zhi.zhao@medisin.uio.no)
## date: 08-July-2021
#######################

rm(list = ls())
library("plyr")
library("data.table")

features <- data.frame(read.csv("gdsc_en_input_w5.csv", head=T))
names.fea <- strsplit(rownames(features), "")
features <- t(features)
p <- c(13321, 13747-13321, 13818-13747)
Cell.Line <- rownames(features)
features <- data.frame(Cell.Line, features)

library(PharmacoGx)
library(Biobase)
availablePSets()
## run on Linux
GDSC <- downloadPSet("GDSC")
ic50 <- summarizeSensitivityProfiles(pSet=GDSC, sensitivity.measure="ic50_published", summary.stat="median")
ic50 <- t(ic50[rownames(ic50) %in% c("Methotrexate","RDEA119","PD-0325901","CI-1040","AZD6244","Nilotinib", "Axitinib"),])

ic50 <- data.frame(Cell.Line=rownames(ic50), ic50)
ic50$Cell.Line <- gsub("-", ".", ic50$Cell.Line)
load("GDSC_BayesSUR.rda") # NOTE that 'GDSC_BayesSUR.rda' was saved from file 'GDSC_preprocess1.R'
# select cell lines not in the training data
y00 <- ic50[!ic50$Cell.Line %in% GDSC$names.cell.line,]
y00 <- y00[!is.na(rowSums(y00[,-1])),]
y00[,-1] <- log(y00[,-1])  

#===============
# combine drug sensitivity, tissues and molecular features
#===============
yx <- merge(y00, features, by="Cell.Line")

# numbers of gene expression features, copy number festures and muatation features
p <- c(13321, 13747-13321, 13818-13747) 
num.nonpen <- 13
yx <- data.matrix(yx[,-1])
y <- yx[,1:7]; y <- y[,c(1,3,6,4,7,2,5)]; colnames(y) <- c("Methotrexate","RDEA119","PD-0325901","CI-1040","AZD6244","Nilotinib", "Axitinib")
x0 <- yx[,dim(y00)[2]-1+sum(p)+1:num.nonpen]
x1 <- yx[,dim(y00)[2]-1+1:p[1]]; x1 <- log2(x1)
x2 <- yx[,dim(y00)[2]-1+p[1]+1:p[2]]; colnames(x2) <- paste(substr(colnames(x2), 1, nchar(colnames(x2))-3), ".CNV", sep="")
x3 <- yx[,dim(y00)[2]-1+p[1]+p[2]+1:p[3]]

ic50.val <- cbind(y,x0,x1,x2,x3)

# NOTE that 'example_GDSC_kegg10%.RData', 'example_GDSC_kegg30%.RData' and 'example_GDSC_kegg50%.RData' was saved from file 'GDSC_preprocess1.R'
load("example_GDSC_kegg10%.RData")
example_GDSC_kegg10.val <- ic50.val[,colnames(ic50.val) %in% colnames(example_GDSC_kegg$data)]
save(example_GDSC_kegg10.val,file="example_GDSC_kegg10%_val.RData")
load("example_GDSC_kegg30%.RData")
example_GDSC_kegg30.val <- ic50.val[,colnames(ic50.val) %in% colnames(example_GDSC_kegg$data)]
save(example_GDSC_kegg30.val,file="example_GDSC_kegg30%_val.RData")
load("example_GDSC_kegg50%.RData")
example_GDSC_kegg50.val <- ic50.val[,colnames(ic50.val) %in% colnames(example_GDSC_kegg$data)]
save(example_GDSC_kegg50.val,file="example_GDSC_kegg50%_val.RData")




