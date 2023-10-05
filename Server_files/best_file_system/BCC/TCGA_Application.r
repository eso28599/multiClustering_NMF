

###Load the data.  Each file can be downloaded from the TCGA Data Portal at link https://tcga-data.nci.nih.gov/docs/publications/brca_2012/ (accessed 12/14/2012)
#For formatting reasons, we convert the GE data (BRCA.exp.348.med) to a .csv file be for loading. Open the file in Microsoft Excel and save it as .csv. 
#All other files can be read as they are givein in the Data Portal.  
GE = read.csv("BRCA.exp.348.med.csv", header = TRUE)
miRNA = read.csv("BRCA.348.precursor.txt", header = TRUE)
Protein = read.table("rppaData-403Samp-171Ab-Trimmed.txt", header = TRUE)
Meth = read.table("BRCA.Methylation.574probes.802.txt", header = TRUE)

################Match columns (samples) between sources##############################
namesExp = names(GE)[2:349]
namesmiRNA = names(miRNA)[2:349]
namesProtein = names(Protein)[2:404]
namesMeth = names(Meth)

namesExp = substr(namesExp,1,16)
namesmiRNA = substr(namesmiRNA,1,16)
namesProtein= substr(namesProtein,1,16)

namesmiRNA[order(namesmiRNA)] == namesExp ##all match
MatchProt = match(namesExp,namesProtein,nomatch=0)
namesProtein[MatchProt] == namesExp ##all match
MatchMeth = match(namesExp,namesMeth,nomatch=0)
namesMeth[MatchMeth] == namesExp ##all match

miRNA = miRNA[,2:349]
miRNA.mat = as.matrix(miRNA[,order(namesmiRNA)])

Protein.mat = Protein[,2:404]
Protein.mat = as.matrix(Protein.mat[,MatchProt])

Meth.mat = as.matrix(Meth[,MatchMeth])

###################Data processing#############################

###Impute missing expression values
#load package 'impute' from CRAN
Exp.mat = as.matrix(GE[,2:349])
Exp.mat = impute.knn(Exp.mat) ##Impute missing values via KNN (K=10)
Exp.mat = Exp.mat$data

rowSums(miRNA.mat==0)
##Remove miRNAs with > 50% missing values
miRNA.mat = miRNA.mat[rowSums(miRNA.mat==0) < 348*0.5,]

X1 = Exp.mat[apply(Exp.mat,1,sd)>1.5,] ###Filter to select only most variable genes
X2 = sqrt(Meth.mat) ##take square root of methylation data
X3 = log(1+miRNA.mat) ##take log of miRNA data
X4 = scale(Protein.mat,center=TRUE,scale=TRUE) #Column center/scale protein

X = list(X1,X2,X3,X4)
#X1 = GE, X2 = Meth, X3 = miRNA, X4 = Protein

#######################Perform BCC###############################################
##Load the BCC function(e.g., enter source("BCC.r") if the BCC.r file is in your working directory)
##Also need to load the 'gtools' package

 
###Fit the model for K=2,...,10 to compare the mean adjusted adherence
##This can take a while (30 min-1 hour for each K)
##To reduce computation time, lower the number of MCMC draws
Results = list()
for(i in 2:10) Results[[i-1]] = BCC(X,K=i,IndivAlpha = TRUE,NumDraws=10000)

meanAlphaStar = c()
upperAlphaStar = c()
lowerAlphaStar = c()
for(i in 1:9){
 K= i+1
 AlphaStar = apply((Results[[i]]$AlphaVec[,2000:10000]-1/K)/(1-1/K),2,mean)
 meanAlphaStar[i] = mean(AlphaStar)
 upperAlphaStar[i] = quantile(AlphaStar,0.975)
 lowerAlphaStar[i] = quantile(AlphaStar,0.025)
}
###Plot AlphaStar (mean adjusted adherence) values for each K
plotCI(c(2:10),meanAlphaStar,ui = upperAlphaStar,li = lowerAlphaStar, xlab = 'Number of clusters (K)',ylab ='Mean adjusted adherence')# , gap = 0.2,xlab = expression(paste("True ", alpha)),ylab = expression(paste("Estimated ", alpha)), ylim = c(0.5,1), xlim =c(0.5,1), sfrac = 0.005,main = "Point estimation")
#Maximized when K=3

#Fit model for K=3 (can take 30 min - 1 hour)
  
K3 = BCC(X,K=3,IndivAlpha = TRUE,NumDraws=10000)

#Get alpha values 
K3$Alpha

#Make principal component plots of clusters 
par(mfrow=c(2,2))
PCs = prcomp(t(X1))
plot(PCs$x[,1],PCs$x[,2], cex = 0, xlab = "PC 1", ylab = "PC 2", main = expression(paste("GE  (", alpha, " = 0.91)")))
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[1]][,1]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[1]][,1]==1,2],pch=16)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[1]][,2]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[1]][,2]==1,2],pch = 3)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[1]][,3]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[1]][,3]==1,2],pch = 8)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[1]][,1]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[1]][,1]==1,2], col = 'red',pch=16)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[1]][,2]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[1]][,2]==1,2], col = 'red',pch = 3)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[1]][,3]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[1]][,3]==1,2], col = 'red',pch = 8)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[1]][,1]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[1]][,1]==1,2], col = 'blue',pch=16)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[1]][,2]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[1]][,2]==1,2], col = 'blue',pch = 3)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[1]][,3]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[1]][,3]==1,2], col = 'blue',pch = 8)

PCs = prcomp(t(X2))
plot(PCs$x[,1],PCs$x[,2], cex = 0, xlab = "PC 1", ylab = "PC 2", main = expression(paste("ME  (", alpha, " = 0.69)")))
m=2
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,1]==1,2],pch=16)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,2]==1,2],pch = 3)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,3]==1,2],pch = 8)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,1]==1,2], col = 'red',pch=16)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,2]==1,2], col = 'red',pch = 3)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,3]==1,2], col = 'red',pch = 8)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,1]==1,2], col = 'blue',pch=16)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,2]==1,2], col = 'blue',pch = 3)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,3]==1,2], col = 'blue',pch = 8)

PCs = prcomp(t(X3))
plot(PCs$x[,1],PCs$x[,2], cex = 0, xlab = "PC 1", ylab = "PC 2", main = expression(paste("miRNA  (", alpha, " = 0.56)")))
m=3
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,1]==1,2],pch=16)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,2]==1,2],pch = 3)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,3]==1,2],pch = 8)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,1]==1,2], col = 'red',pch=16)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,2]==1,2], col = 'red',pch = 3)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,3]==1,2], col = 'red',pch = 8)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,1]==1,2], col = 'blue',pch=16)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,2]==1,2], col = 'blue',pch = 3)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,3]==1,2], col = 'blue',pch = 8)

PCs = prcomp(t(X4))
plot(PCs$x[,1],PCs$x[,2], cex = 0, xlab = "PC 1", ylab = "PC 2", main = expression(paste("RPPA  (", alpha, " = 0.70)")))
m=4
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,1]==1,2],pch=16)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,2]==1,2],pch = 3)
points(PCs$x[K3$Cbest[,1]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,1]==1& K3$Lbest[[m]][,3]==1,2],pch = 8)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,1]==1,2], col = 'red',pch=16)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,2]==1,2], col = 'red',pch = 3)
points(PCs$x[K3$Cbest[,2]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,2]==1& K3$Lbest[[m]][,3]==1,2], col = 'red',pch = 8)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,1]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,1]==1,2], col = 'blue',pch=16)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,2]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,2]==1,2], col = 'blue',pch = 3)
points(PCs$x[K3$Cbest[,3]==1 & K3$Lbest[[m]][,3]==1,1],PCs$x[K3$Cbest[,3]==1& K3$Lbest[[m]][,3]==1,2], col = 'blue',pch = 8)

###Find cluster matching matrices
###The TCGA clusterings were obtained from Supplemental Table 1 in the 2012 Nature article "Comprehensive molecular portraits of human breast tumours"
###The table can be downloaded at this link: http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html#/supplementary-information *accessed 12/14/2012
###The table should be saved as a 'csv' file before loading
Table = read.csv("Table1Nature.csv",header = TRUE)
Subtype = Table$Integrated.Clusters..with.PAM50.##Comprehensive subtypes
###Match IDs
TableIDs = Table$Complete.TCGA.ID 
TableIDs = gsub('-','.', as.character(TableIDs))
namesExp2 = substr(namesExp,1,12)
MatchSub = match(namesExp2,TableIDs,nomatch=0)
sum(TableIDs[MatchSub] == namesExp2)
Subtype = Subtype[MatchSub]

###get and match source-specific subtypes
GESubtype = Table$PAM50.mRNA
GESubtype = GESubtype[MatchSub]

MethSubtype = Table$methylation.Clusters
MethSubtype = MethSubtype[MatchSub]

miRNASubtype = Table$miRNA.Clusters
miRNASubtype = miRNASubtype[MatchSub]

RPPASubtype = Table$RPPA.Clusters
RPPASubtype = RPPASubtype[MatchSub]
##Group "LumA' and "LumA/B"
RPPASubtype[as.numeric(RPPASubtype) == 3] = levels(RPPASubtype)[4]

Clusters = rep(1,348)
Clusters[K3$Cbest[,2]==1] = 2
Clusters[K3$Cbest[,3]==1] = 3
##ConfusionMat - All
confmat = matrix(nrow = 4, ncol = 3)
for(i in 1:4){
for(j in 1:3){
confmat[i,j] = sum(Subtype == i & Clusters ==j) }} 

##ConfusionMat - Exp
ClustersGE = rep(1,348)
ClustersGE[K3$Lbest[[1]][,2]==1] = 2
ClustersGE[K3$Lbest[[1]][,3]==1] = 3
confmat = matrix(nrow = 5, ncol = 3)
for(i in 1:5){
for(j in 1:3){
confmat[i,j] = sum(as.numeric(GESubtype) == i & ClustersGE ==j) }} 

##ConfusionMat - Meth
ClustersMeth = rep(1,348)
ClustersMeth[K3$Lbest[[2]][,2]==1] = 2
ClustersMeth[K3$Lbest[[2]][,3]==1] = 3
confmat = matrix(nrow = 5, ncol = 3)
for(i in 1:5){
for(j in 1:3){
confmat[i,j] = sum(MethSubtype == i & ClustersMeth ==j) }} 

##ConfusionMat - miRNA
ClustersmiRNA = rep(1,348)
ClustersmiRNA[K3$Lbest[[3]][,2]==1] = 2
ClustersmiRNA[K3$Lbest[[3]][,3]==1] = 3
confmat = matrix(nrow = 7, ncol = 3)
for(i in 1:7){
for(j in 1:3){
confmat[i,j] = sum(miRNASubtype == i & ClustersmiRNA ==j) }} 

#ConfusionMat - RPPA
ClustersRPPA = rep(1,348)
ClustersRPPA[K3$Lbest[[4]][,2]==1] = 2
ClustersRPPA[K3$Lbest[[4]][,3]==1] = 3
confmat = matrix(nrow = 6, ncol = 3)
for(i in c(1:6)){ ###ignore 'Lum A' and  'X' factors
for(j in 1:3){
confmat[i,j] = sum(as.numeric(RPPASubtype) == i & ClustersRPPA ==j) }} 


######### Heatmaps - Expression

p = dim(X1)[1]
pre.gene = t(X1)

subtype = ClustersGE
n = length(subtype)

############################
### 2.centering datasets ###
############################

gene = pre.gene-matrix(rep(apply(pre.gene,2,mean,na.rm=T),each=n),nrow=n)

#################################################
###   Functions for hierarchical clustering   ###
#################################################

row.cluster<-function(data){
	d.data<-dist(data)
	h.data<-hclust(d.data)
	ordered.data<-data[h.data$order,]
	return(ordered.data)
}

col.cluster<-function(data){
	d.data<-dist(t(data))
	h.data<-hclust(d.data)
	ordered.data<-data[,h.data$order]
	return(ordered.data)
}

######################################################
### 3.Gene data ordered by subtypes and clustering ###
######################################################

types<-levels(as.factor(subtype))
K = length(types)
ordered.x = c()
sub.n = c()
for(i in 1:K){
	typedx = gene[subtype==types[i],]
	sub.n = c(sub.n,dim(typedx)[1])
	ordered.typedx = row.cluster(typedx)
	ordered.x = rbind(ordered.x,ordered.typedx)
}


clustered.gene = col.cluster(ordered.x)
final.x = clustered.gene


## thresholding 
final.x[final.x>4]=4 
final.x[final.x<(-4)]=-4 

###########################################
### 	   4.Heatmap of gene	    	    ###
###########################################

layout(matrix(c(1:2),2,1), heights = c(1,6), respect = FALSE)

### color key ###
grid=seq(-4,4,by=8/100)
par(mar=c(1.5,0.5,2,20))
image(x=grid,y=1:5,matrix(rep(grid,5),length(grid),5),col=greenred(100),axes=FALSE,main="Color Key",xlab="",ylab="",cex.main=1)
axis(1,seq(-4,4,by=1),seq(-4,4,by=1))

### title of heatmap ###
mtext("GE with subtypes",side=4,line=4,las=1,cex=1.7,font=2)

### Heatmap ###

par(mar=c(5,0.5,1,3))
image(x=1:n, y=1:p, final.x[,p:1], axes=FALSE,col=greenred(100),xlab="",ylab="")		
axis(1, cumsum(sub.n)-sub.n/2,types, tick=F, cex.axis=1.3,col.axis="blue",font.axis=3)
axis(1, c(0,cumsum(sub.n))+0.5, rep(" ",K+1))
abline(v=cumsum(sub.n)+0.5,col="white")
mtext("Samples",side=1,line=3,cex=1.5)
mtext("Genes",side=4,line=1,cex=1.5)






######### Heatmaps - Meth
hist(X2)
p = dim(X2)[1]
pre.gene = t(X2)

subtype = ClustersMeth
n = length(subtype)

############################
### 2.centering datasets ###
############################

gene = pre.gene-matrix(rep(apply(pre.gene,2,mean,na.rm=T),each=n),nrow=n)

######################################################
### 3.Gene data ordered by subtypes and clustering ###
######################################################

types<-levels(as.factor(subtype))
K = length(types)
ordered.x = c()
sub.n = c()
for(i in 1:K){
	typedx = gene[subtype==types[i],]
	sub.n = c(sub.n,dim(typedx)[1])
	ordered.typedx = row.cluster(typedx)
	ordered.x = rbind(ordered.x,ordered.typedx)
}


clustered.gene = col.cluster(ordered.x)
final.x = clustered.gene

## thresholding 
final.x[final.x>0.4]= 0.4
final.x[final.x<(-0.4)]=-0.4 

###########################################
### 	   4.Heatmap of gene	    	    ###
###########################################

layout(matrix(c(1:2),2,1), heights = c(1,6), respect = FALSE)

### color key ###
grid=seq(-0.4,0.4,by=0.08/100)
par(mar=c(1.5,0.5,2,20))
image(x=grid,y=1:5,matrix(rep(grid,5),length(grid),5),col=greenred(100),axes=FALSE,main="Color Key",xlab="",ylab="",cex.main=1)
axis(1,seq(-.4,0.4,by=0.1),seq(-.4,0.4,by=0.1))

### title of heatmap ###
mtext("ME with subtypes",side=4,line=4,las=1,cex=1.7,font=2)

### Heatmap ###

par(mar=c(5,0.5,1,3))
image(x=1:n, y=1:p, final.x[,p:1], axes=FALSE,col=greenred(100),xlab="",ylab="")		
axis(1, cumsum(sub.n)-sub.n/2,types, tick=F, cex.axis=1.3,col.axis="blue",font.axis=3)
axis(1, c(0,cumsum(sub.n))+0.5, rep(" ",K+1))
abline(v=cumsum(sub.n)+0.5,col="white")
mtext("Samples",side=1,line=3,cex=1.5)
mtext("Probes",side=4,line=1,cex=1.5)


######### Heatmaps - miRNA
p = dim(X3)[1]
pre.gene = t(X3)

subtype = ClustersmiRNA
n = length(subtype)

############################
### 2.centering datasets ###
############################

gene = pre.gene-matrix(rep(apply(pre.gene,2,mean,na.rm=T),each=n),nrow=n)

######################################################
### 3.Gene data ordered by subtypes and clustering ###
######################################################

types<-levels(as.factor(subtype))
K = length(types)
ordered.x = c()
sub.n = c()
for(i in 1:K){
	typedx = gene[subtype==types[i],]
	sub.n = c(sub.n,dim(typedx)[1])
	ordered.typedx = row.cluster(typedx)
	ordered.x = rbind(ordered.x,ordered.typedx)
}


clustered.gene = col.cluster(ordered.x)
final.x = clustered.gene

## thresholding 
final.x[final.x>2]= 2
final.x[final.x<(-2)]=-2 

###########################################
### 	   4.Heatmap of gene	    	    ###
###########################################

layout(matrix(c(1:2),2,1), heights = c(1,6), respect = FALSE)

### color key ###
grid=seq(-2,2,by=5*0.08/100)
par(mar=c(1.5,0.5,2,20))
image(x=grid,y=1:5,matrix(rep(grid,5),length(grid),5),col=greenred(100),axes=FALSE,main="Color Key",xlab="",ylab="",cex.main=1)
axis(1,seq(-2,2,by=5*0.1),seq(-2,2,by=5*0.1))

### title of heatmap ###
mtext("miRNA with subtypes",side=4,line=4,las=1,cex=1.7,font=2)

### Heatmap ###

par(mar=c(5,0.5,1,3))
image(x=1:n, y=1:p, final.x[,p:1], axes=FALSE,col=greenred(100),xlab="",ylab="")		
axis(1, cumsum(sub.n)-sub.n/2,types, tick=F, cex.axis=1.3,col.axis="blue",font.axis=3)
axis(1, c(0,cumsum(sub.n))+0.5, rep(" ",K+1))
abline(v=cumsum(sub.n)+0.5,col="white")
mtext("Samples",side=1,line=3,cex=1.5)
mtext("miRNAs",side=4,line=1,cex=1.5)



######### Heatmaps - RPPA
p = dim(X4)[1]
pre.gene = t(X4)

subtype = ClustersRPPA
n = length(subtype)

############################
### 2.centering datasets ###
############################

gene = pre.gene-matrix(rep(apply(pre.gene,2,mean,na.rm=T),each=n),nrow=n)

######################################################
### 3.Gene data ordered by subtypes and clustering ###
######################################################

types<-levels(as.factor(subtype))
K = length(types)
ordered.x = c()
sub.n = c()
for(i in 1:K){
	typedx = gene[subtype==types[i],]
	sub.n = c(sub.n,dim(typedx)[1])
	ordered.typedx = row.cluster(typedx)
	ordered.x = rbind(ordered.x,ordered.typedx)
}


clustered.gene = col.cluster(ordered.x)
final.x = clustered.gene

## thresholding 
final.x[final.x>2]= 2
final.x[final.x<(-2)]=-2 

###########################################
### 	   4.Heatmap of gene	    	    ###
###########################################

layout(matrix(c(1:2),2,1), heights = c(1,6), respect = FALSE)

### color key ###
grid=seq(-2,2,by=5*0.08/100)
par(mar=c(1.5,0.5,2,20))
image(x=grid,y=1:5,matrix(rep(grid,5),length(grid),5),col=greenred(100),axes=FALSE,main="Color Key",xlab="",ylab="",cex.main=1)
axis(1,seq(-2,2,by=5*0.1),seq(-2,2,by=5*0.1))

### title of heatmap ###
mtext("RPPA with subtypes",side=4,line=4,las=1,cex=1.7,font=2)

### Heatmap ###

par(mar=c(5,0.5,1,3))
image(x=1:n, y=1:p, final.x[,p:1], axes=FALSE,col=greenred(100),xlab="",ylab="")		
axis(1, cumsum(sub.n)-sub.n/2,types, tick=F, cex.axis=1.3,col.axis="blue",font.axis=3)
axis(1, c(0,cumsum(sub.n))+0.5, rep(" ",K+1))
abline(v=cumsum(sub.n)+0.5,col="white")
mtext("Samples",side=1,line=3,cex=1.5)
mtext("Proteins",side=4,line=1,cex=1.5)













