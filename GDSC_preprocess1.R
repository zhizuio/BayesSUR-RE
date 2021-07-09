#######################
## This script is for preprocessing GDSC data 
## "Structured multivariate Bayesian variable selection for pharmacogenomic studies"
## 
## author: Zhi Zhao (zhi.zhao@medisin.uio.no)
## date: 08-July-2021
#######################

rm(list = ls())
library("plyr")
library("data.table")

#===============
# The user needs load the datasets from https://www.cancerrxgene.org/downloads/bulk_download archived data folder 'release-5.0'. 
# Downloading the three datasets used for our analysis.
#===============

features <- data.frame(read.csv("gdsc_en_input_w5.csv", head=T))
names.fea <- strsplit(rownames(features), "")
features <- t(features)
p <- c(13321, 13747-13321, 13818-13747)
Cell.Line <- rownames(features)
features <- data.frame(Cell.Line, features)

ic50_00 <- data.frame(read.csv("gdsc_drug_sensitivity_fitted_data_w5.csv", head=T))
ic50_0 <- ic50_00[,c(1,4,7)]
drug.id <- data.frame(read.csv("gdsc_tissue_output_w5.csv", head=T))[,c(1,3)]
drug.id2 <- drug.id[!duplicated(drug.id$drug.id),]
# delete drug.id=1066 since ID1066 and ID156 both correspond drug AZD6482, 
# and no ID1066 in the "suppl.Data1" by Garnett et al. (2012)
drug.id2 <- drug.id2[drug.id2$drug.id!=1066,] 
drug.id2$drug.name <- as.character(drug.id2$drug.name)
drug.id2$drug.name <- substr(drug.id2$drug.name, 1, nchar(drug.id2$drug.name)-6)
drug.id2$drug.name <- gsub(" ", "-", drug.id2$drug.name)

ic50 <- ic50_0
# mapping the drug_id to drug names in drug sensitivity data set
ic50$drug_id <- plyr::mapvalues(ic50$drug_id, from = drug.id2[,2], to = drug.id2[,1])
colnames(ic50) <- c("Cell.Line", "compound", "IC50")

# transform drug sensitivity overall cell lines to a data matrix
y0 <- reshape(ic50, v.names="IC50", timevar="compound", idvar="Cell.Line", direction="wide")
y0$Cell.Line <- gsub("-", ".", y0$Cell.Line)

#===============
# select nonmissing pharmacological data
#===============
y00 <- y0
m0 <- dim(y0)[2]-1
eps <- 0.05
# r1.na is better to be not smaller than r2.na
r1.na <- 0.3
r2.na <- 0.2
k <- 1
while(sum(is.na(y0[,2:(1+m0)]))>0){
  r1.na <- r1.na - eps/k
  r2.na <- r1.na - eps/k
  k <- k + 1
  ## select drugs with <30% (decreasing with k) missing data overall cell lines
  na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y<r1.na)<m0){
    y0 <- y0[,-c(1+which(na.y>=r1.na))]
    m0 <- sum(na.y<r1.na)
    na.y <- apply(y0[,2:(1+m0)], 2, function(xx) sum(is.na(xx))/length(xx))
  }
  
  ## select cell lines with treatment of at least 80% (increasing with k) drugs
  na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  while(sum(na.y0<r2.na)<(dim(y0)[1])){
    y0 <- y0[na.y0<r2.na,]
    na.y0 <- apply(y0[,2:(1+m0)], 1, function(xx) sum(is.na(xx))/length(xx))
  }
  num.na <- sum(is.na(y0[,2:(1+m0)]))
  message("#{NA}=", num.na, "\n", "r1.na =", r1.na, ", r2.na =", r2.na, "\n")
}

#===============
# combine drug sensitivity, tissues and molecular features
#===============
yx <- merge(y0, features, by="Cell.Line")
names.cell.line <- yx$Cell.Line
names.drug <- colnames(yx)[2:(dim(y0)[2])]
names.drug <- substr(names.drug, 6, nchar(names.drug))
# numbers of gene expression features, copy number festures and muatation features
p <- c(13321, 13747-13321, 13818-13747) 
num.nonpen <- 13
yx <- data.matrix(yx[,-1])
y <- yx[,1:(dim(y0)[2]-1)]
x <- cbind(yx[,dim(y0)[2]-1+sum(p)+1:num.nonpen], yx[,dim(y0)[2]-1+1:sum(p)])

## preselect gene expression explaining 10%/30%/50% variations over all samples
var.x1 <- apply(log(x[,num.nonpen+1:p[1]]), 2, var)
var.sort <- sort(var.x1, decreasing=TRUE)
sum.a <- cumsum(var.sort)
#half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]*0.5))[1]] # explaining 50% variations
#half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]*0.3))[1]] # explaining 30% variations
half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]*0.1))[1]] # explaining 10% variations
x1 <- x[,num.nonpen+1:p[1]][,var.x1>=half.a]
x <- cbind(x[,1:num.nonpen], x1, x[,num.nonpen+p[1]+1:(p[2]+p[3])])
p[1] <- dim(x1)[2]

# delete genes with only one mutated cell line
x <- x[,-c(num.nonpen+p[1]+p[2]+which(colSums(x[,num.nonpen+p[1]+p[2]+1:p[3]])<=1))]
p[3] <- ncol(x) - num.nonpen - p[1] - p[2]

GDSC <- list(y=y, x=x, p=p, num.nonpen=num.nonpen, names.cell.line=names.cell.line, 
             names.drug=names.drug)
save(GDSC, "GDSC_BayesSUR.rda")

##================
##================
## select a small set of drugs
##================
##================

name_drugs <- c("Methotrexate","RDEA119","PD-0325901","CI-1040",
                "AZD6244","Nilotinib", "Axitinib")

# extract the drugs' pharmacological profiling and tissue dummy
col_filter <- colnames(GDSC$y) %in% paste("IC50.", name_drugs,sep="")
YX0 <- cbind(
  GDSC$y[, col_filter][, c(1, 3, 6, 4, 7, 2, 5)],
  GDSC$x[, 1:GDSC$num.nonpen]
)
colnames(YX0) <- c(name_drugs, colnames(GDSC$x)[1:GDSC$num.nonpen])
# extract the genetic information of CNV & MUT
X23 <- GDSC$x[, GDSC$num.nonpen+GDSC$p[1]+1:(p[2]+p[3])]
colnames(X23)[1:p[2]] <- paste(substr(colnames(X23)[1:p[2]], 1, 
                                      nchar(colnames(X23)[1:p[2]] )-3), ".CNV", sep="")

# locate all genes with GEX, CNV or MUT information
name_genes_duplicate <- c( colnames(x1),
  substr(colnames(X23)[1:p[2]], 1, nchar(colnames(X23)[1:p[2]])-4),
  substr(colnames(X23)[p[2]+1:p[3]], 1, nchar(colnames(X23)[p[2]+1:p[3]])-4)
)
name_genes <- name_genes_duplicate[!duplicated(name_genes_duplicate)]

# log2 transfermation for the expression levels
X1 <- log2(x1)

# summary the data information
example_GDSC <- list( data=cbind( YX0, X1, X23 ) )
p <- ncol(X1)+ncol(X23)
example_GDSC$blockList <- list(1:length(name_drugs),
                               length(name_drugs)+1:GDSC$num.nonpen, 
                               ncol(YX0)+1:p)

#========================
# construct the G matrix: edge potentials in the MRF prior
#========================

# edges between drugs: Group1 ("RDEA119","17-AAG","PD-0325901","CI-1040", "AZD6244") 
# indexed as (2:5)

# # The MAPK_pathway.txt file is originally from KEGG or GSEA database at
# # http://software.broadinstitute.org/gsea/msigdb/cards/KEGG_MAPK_SIGNALING_PATHWAY
pathway_genes <- read.table("MAPK_pathway.txt")[[1]]
Idx_Pathway1 <- which(name_genes_duplicate %in% pathway_genes)
rep1 <- rep(Idx_Pathway1, each=length(2:5))
rep2 <- rep((2:5-1) * sum(p), times=length(Idx_Pathway1))
rep3 <- rep1 + rep2
Gmrf_Group1Pathway1 <- t(combn(rep3, 2))

# edges between drugs: Group2 ("Nilotinib","Axitinib") indexed as (6:7)
# delete gene ABL2 which is not related to BCR-ABL fusion gene
#Idx_Pathway2 <- which(name_genes_duplicate %like% "BCR" | name_genes_duplicate %like% "ABL")[-c(1:2,4)] # for top 50% cumVar[GEX]
#Idx_Pathway2 <- which(name_genes_duplicate %like% "BCR" | name_genes_duplicate %like% "ABL")[-c(1,3)] # for top 30% cumVar[GEX]
Idx_Pathway2 <- which(name_genes_duplicate %like% "BCR" | name_genes_duplicate %like% "ABL")[-2] # for top 10% cumVar[GEX]
Gmrf_Group2Pathway2 <- t(combn(rep(Idx_Pathway2,each=length(6:7)) + 
                                 rep((6:7-1)*sum(p),times=length(Idx_Pathway2)), 2))

# edges between the common gene in different data sources
Gmrf_CommonGene <- NULL
list_CommonGene <- list(0)
k <- 1
for(i in 1:length(name_genes)){
  Idx_CommonGene <- which(name_genes_duplicate == name_genes[i])
  if(length(Idx_CommonGene) > 1){
    Gmrf_CommonGene <- rbind(Gmrf_CommonGene,t(combn(rep(Idx_CommonGene,each=length(name_drugs))
                                                     + rep((1:length(name_drugs)-1)*sum(p),times=length(Idx_CommonGene)), 2)))
    k <- k+1
  }
}
Gmrf <- rbind(  Gmrf_Group1Pathway1, Gmrf_Group2Pathway2, Gmrf_CommonGene )
example_GDSC$mrfG <- Gmrf
example_GDSC_kegg <- example_GDSC

#save(example_GDSC_kegg, file="example_GDSC_kegg50%.RData") # save dataset with top 30% cumVar[GEX]
#save(example_GDSC_kegg, file="example_GDSC_kegg30%.RData") # save dataset with top 30% cumVar[GEX]
save(example_GDSC_kegg, file="example_GDSC_kegg10%.RData") # save dataset with top 10% cumVar[GEX]

# create the target gene names of the two groups of drugs
targetGenes1 <- matrix(Idx_Pathway1,nrow=1)
colnames(targetGenes1) <- colnames(example_GDSC$data)[7+13+Idx_Pathway1]
targetGenes2 <- matrix(Idx_Pathway2,nrow=1)
colnames(targetGenes2) <- colnames(example_GDSC$data)[7+13+Idx_Pathway2]

targetGene <- list(group1=targetGenes1, group2=targetGenes2)
#save(targetGene, file="targetGene_kegg_50%.RData") # save target genes of MAPK pathway from the dataset of top 50% cumVar[GEX]
#save(targetGene, file="targetGene_kegg_30%.RData") # save target genes of MAPK pathway from the dataset of top 30% cumVar[GEX]
save(targetGene, file="targetGene_kegg_10%.RData") # save target genes of MAPK pathway from the dataset of top 10% cumVar[GEX]

