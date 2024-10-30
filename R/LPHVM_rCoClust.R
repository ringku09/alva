

RbstCoClust <- function(DtMat)
	{
	 source("LPHVM_Pack.R")
	 LPHVM_Res <- LPHVM(DtMat)
	 nCoCluster <- LPHVM_Res$k
	 GeneCluster <- LPHVM_Res$GclsMem
	 CompCluster <- LPHVM_Res$CclsMem
	 CoClustMean <- LPHVM_Res$Co_ClsMean
	 UpDEgene <- LPHVM_Res$UpRgene
	 DownDEgene <- LPHVM_Res$DnRgene
	 RnkToxCC <- LPHVM_Res$ToxCCrank
	 RnkGCCrel <- LPHVM_Res$CC_Gene_JproRank
	 GCjointProb <- LPHVM_Res$Co_cls_JointProb

return(list(nCoCluster=nCoCluster, GeneCluster=GeneCluster, CompCluster=CompCluster,
	 CoClustMean=CoClustMean, UpDEgene=UpDEgene,DownDEgene=DownDEgene, RnkToxCC=RnkToxCC,
	 RnkGCCrel=RnkGCCrel,  GCjointProb=GCjointProb))
	}
