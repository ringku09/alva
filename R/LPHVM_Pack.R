
##--------------------------------------------------------------------------------------------------------------##
##-------------------------------------------------Real Data----------------------------------------------------##
source("PHVM_Package_1.R")
source("BioMCoClustr_Package2.R")
library(ggplot2)
library(factoextra)
library(NbClust)

get_k <- function(DtMat, seed = 123) {
  set.seed(seed)
  GCmat <- as.matrix(round(100*(1/(1+exp(-abs(DtMat))))))
  df_TH_D2 <- scale(t(GCmat))
  NclasTHD2 <- fviz_nbclust(df_TH_D2, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
  g <- NclasTHD2$data[,"gap"]
  d <- NclasTHD2$data[,"gap"]-NclasTHD2$data[,"SE.sim"]
  rr <- cbind(g,d)
  kk <- rep(0,9)
  for (i in 1:nrow(rr)-1)
  {
    kk[i] <- g[i]>=d[i+1]
  }
  k <- which(kk>0)[1]
  return(k)
}





LPHVM <- function (DtMat, k)
{
GCmat <- as.matrix(round(100*(1/(1+exp(-abs(DtMat))))))
##--------------------------------------------------------------------------------------------------------------##
CoCls <- k
ResultPHVM <- PHVM(GCmat,k)                           #names(ResultPHVM)
GCjointProb <- ResultPHVM$GCjointProb
GclsMem <- ResultPHVM$Gclust
CclsMem <- ResultPHVM$Cclust
BioMCls <- BioMCocls(GCmat,CoCls,GCjointProb,GclsMem,CclsMem)      #names(BioMCls)
Co_ClsLvl <- BioMCls$ClsLvl
Co_ClsMean <- BioMCls$CoClsMean
# bioMgene <- sapply(strsplit(rownames(as.matrix(unlist(GclsMem[-which(Co_ClsMean==min(Co_ClsMean))]))), ".", fixed=TRUE),function(x) x[2])
# bioMCC <- sapply(strsplit(rownames(as.matrix(unlist(CclsMem[-which(Co_ClsMean==min(Co_ClsMean))]))), ".", fixed=TRUE),function(x) x[2])
Co_cls_JointProb <- BioMCls$Co_cls_JointProb
# bioMfcDt <- DtMat[bioMgene,bioMCC]
# UpRgene <- rownames(as.matrix(which(rowMeans(bioMfcDt)>0)))
# DnRgene <- rownames(as.matrix(which(rowMeans(bioMfcDt)<0)))
bioMjntProb <- Co_cls_JointProb[bioMgene,bioMCC]
ToxCCrank <- as.matrix(sort((colMeans(bioMjntProb)/max(colMeans(bioMjntProb)))*100, decreasing = T))
ToxCC <- rep(colnames(bioMjntProb), each = nrow(bioMjntProb))
BioGene <- rep(rownames(bioMjntProb),ncol(bioMjntProb))
JprobScore <- (as.vector(bioMjntProb)/max(as.vector(bioMjntProb)))*100
GeneCCJprob <- cbind(ToxCC,BioGene,JprobScore)  #class(GeneCCJprob)
CC_Gene_JproRank <- GeneCCJprob [order(as.numeric(GeneCCJprob[,3]), decreasing = T),]

return(list(k = k, GclsMem=GclsMem, Co_ClsMean=Co_ClsMean, CclsMem=CclsMem, UpRgene=UpRgene,
	DnRgene=DnRgene, ToxCCrank=ToxCCrank, CC_Gene_JproRank=CC_Gene_JproRank, Co_cls_JointProb=Co_cls_JointProb))
}

##----------------------------------------------Simulation Structure----------------------------------------------------------##

