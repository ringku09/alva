BioMCocls <- function (GCmat, CoCls, GCjointProb, GclsMem, CclsMem)
	{
	 CoClsMean <- rep(0,CoCls)
	 for (cls in 1:CoCls)
		{
	 	 CoClsMean[cls] <- mean(GCjointProb[GclsMem[[cls]],CclsMem[[cls]]])
		}
	 ClsLvl <- paste("Co-Cluster", 1:CoCls, sep=".")
	 Co_clsGMem <- sapply(strsplit(rownames(as.matrix(unlist(GclsMem))),".",fixed=TRUE),function(x)x[2])
	 Co_clsCMem <- sapply(strsplit(rownames(as.matrix(unlist(CclsMem))),".",fixed=TRUE),function(x)x[2])
	 Co_cls_JointProb <- GCjointProb[unique(Co_clsGMem), Co_clsCMem]

return(list(Co_cls_JointProb = Co_cls_JointProb, ClsLvl = ClsLvl, CoClsMean=CoClsMean, Co_clsGMem=Co_clsGMem))
}



