

PHVM <- function(GCmat,k){
PCCVarVec <- NULL
PGVarVec <- NULL
CC_G_AveVarVec <- NULL

GCcoOccuNewVec = list()
clusterCCvec = list()
clusterGvec = list()
pCgivenZnewVec = list()
pGgivenZnewVec = list()
CompPriorVec = list()
GenePriorVec = list()
Hclass <- k
for (i in 1:10)
{
#PZ <- runif(Hclass)
#ProbZ <- PZ/sum(PZ)
ProbZ <- as.vector(gtools::rdirichlet(1, sample(10:150,Hclass)))
CompPrior <- gtools::rdirichlet(Hclass, c(colSums(GCmat)))
GenePrior <- gtools::rdirichlet(Hclass, c(rowSums(GCmat)))
CompPriorVec[[i]] = CompPrior
GenePriorVec[[i]] = GenePrior

##---------------------------- Structuring Initial Probability for pHSA --------------------------##
pZ <- diag(ProbZ)
pCgivenZ <- t(CompPrior)
pGgivenZ <- t(GenePrior)
##------------------------------------------------------------------------------------------------##

##------------------- EM Algorithm for Estimating Parameters of pHSA stability -------------------##
d <- 0.5
iter <- 0
while (d > 1e-05)
   	{

################################## E-Step ####################################
	  nHiClass <- length(ProbZ)
	  pZgGC <- list()
		for (h in 1:nHiClass)
    			{
      		  pZgGC[[h]] <- as.matrix(pGgivenZ[,h])%*%as.matrix(pZ[h,h])%*%t(pCgivenZ[,h])
    			}
	  pZgivenGC <- rapply(pZgGC,function(x){x/Reduce('+',pZgGC)},how="list")
#----------------------------------------------------------------------------------------------------#
	  GCcoOccuOld <- pGgivenZ%*%pZ%*%t(pCgivenZ)
	  logLGCold <- log(GCcoOccuOld)
	  totalLogLGCold <- sum(GCmat*logLGCold) #dim(GCmat) #dim(logLGCold)

################################## M-Step ####################################
   	  GCmatpZgGC <- rapply(pZgivenGC,function(x,GCmat){GCmat*x},how="list",GCmat=GCmat)
   	  pZn <- rapply(GCmatpZgGC,function(x){sum(x)/sum(Reduce('+', GCmatpZgGC))},how="list")
   	  pZnew <- diag(unlist(pZn), ncol=ncol(pZ))
   	  pCgivenZn <- rapply(GCmatpZgGC,function(x){colSums(x)/sum(x)},how="list")
   	  pCgivenZnew <- matrix(unlist(pCgivenZn), ncol =ncol(pCgivenZ))
   	  pGgivenZn <- rapply(GCmatpZgGC,function(x){rowSums(x)/sum(x)},how="list")
   	  pGgivenZnew <- matrix(unlist(pGgivenZn), ncol =ncol(pGgivenZ))
#-----------------------------------------------------------------------------------------------------#
   	  GCcoOccuNew <- pGgivenZnew%*%pZnew%*%t(pCgivenZnew)
	  dimnames(GCcoOccuNew) <- list(rownames(GCmat),colnames(GCmat))
	  logLGCNew <- log(GCcoOccuNew)
   	  totalLogLGCnew <- sum(GCmat*logLGCNew)

        d <- totalLogLGCnew-totalLogLGCold

        totalLogLGCold <- totalLogLGCnew
        pZ <- pZnew
        pCgivenZ <- pCgivenZnew
        pGgivenZ <- pGgivenZnew
        iter = iter+1
      }
##--------------------------------------- End of EM --------------------------------------------------##

dimnames(pCgivenZnew)<-list(colnames(GCmat),paste('Class:',1:ncol(pCgivenZnew),sep=''))
dimnames(pGgivenZnew)<-list(rownames(GCmat),paste('Class:',1:ncol(pGgivenZnew),sep=''))

##--------------------------- Chemical Compound (CC)Clustering----------------------------------------##

CCrankOclass <- apply(pCgivenZnew,1,function(x)which(x==max(x)))
clusterCC <- lapply(1:Hclass, function(x)NULL)
names(clusterCC) <- paste("Class",1:Hclass,sep="")
for( j in 1:ncol( pCgivenZnew))
	{
	  clusterCC[[j]] <- which(CCrankOclass==j)[order(pCgivenZnew[which(CCrankOclass==j),j],decreasing =F)]
	}
clusterCCvec[[i]] = clusterCC

GrankOclass <- apply(pGgivenZnew,1,function(x)which(x==max(x)))
clusterG <- lapply(1:Hclass, function(x)NULL)
names(clusterG) <- paste("Class",1:Hclass,sep="")
for( g in 1:ncol(pGgivenZnew))
	{
	  clusterG[[g]] <- which(GrankOclass==g)[order(pGgivenZnew[which(GrankOclass==g),g],decreasing =F)]
	}
clusterGvec[[i]] = clusterG

CCClassESS = NULL
GClassESS = NULL
nCC = NULL
nG = NULL
poolVar = NULL
for (k in 1:nHiClass)
	{
	  CCClassESS[k] = var(as.vector(GCmat[,clusterCC[[k]]]))
	  GClassESS[k] = var(as.vector(GCmat[clusterG[[k]],]))
	  nCC[k] = length(as.vector(GCmat[,clusterCC[[k]]]))
	  nG[k] = length(as.vector(GCmat[clusterG[[k]],]))
	}
InitialOpt <- function(CCvar, Gvar, nCC, nG)
	{
	 CCV <- NULL
	 GV <- NULL
	 CC_G_comVar <- NULL
	 for (v in 1:nHiClass)
		{
		 CCV[v] <- CCvar[v]*(nCC[v]-1)
		 GV[v] <-  Gvar[v]*(nG[v]-1)
		}
	 PCCVar <- sum(CCV,na.rm=T)/(sum(nCC)-nHiClass)
	 PGVar <- sum(GV,na.rm=T)/(sum(nG)-nHiClass)
	 CC_G_comVar[1] <-  PCCVar
	 CC_G_comVar[2] <-  PGVar
	 CC_G_AveVar <- mean(CC_G_comVar,na.rm=T)
	 return(list(PCCVar = PCCVar,PGVar =PGVar ,CC_G_AveVar = CC_G_AveVar))
	}
resultCC_G_CCGvar <- InitialOpt(CCClassESS, GClassESS, nCC, nG)
PCCVarVec <- c(PCCVarVec,resultCC_G_CCGvar$PCCVar)
PGVarVec <- c(PGVarVec,resultCC_G_CCGvar$PGVar)
CC_G_AveVarVec <- c(CC_G_AveVarVec,resultCC_G_CCGvar$CC_G_AveVar)

GCcoOccuNewVec[[i]] = GCcoOccuNew
pCgivenZnewVec[[i]] = pCgivenZnew
pGgivenZnewVec[[i]] = pGgivenZnew
OptIter <- which(CC_G_AveVarVec == min(CC_G_AveVarVec))
OptIterG <- which(PGVarVec == min(PGVarVec))
}
return(list(GCjointProb = GCcoOccuNewVec[[OptIter[1]]], GzProb = pGgivenZnewVec[[OptIter[1]]], CzProb = pCgivenZnewVec[[OptIter[1]]],
Gclust = clusterGvec[[OptIter[1]]], Cclust = clusterCCvec[[OptIter[1]]]))
}
