spectralClustering <- function(affinity, K, type=3) {
  
  ###This function implements the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
  ###THe inputs are as follows:
  
  #affinity: the similarity matrix;
  #K: the number of clusters
  # type: indicators of variants of spectral clustering 
  
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1 / d)
    NL = Di %*% L
  } else if(type == 3) {
    Di = diag(1 / sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  
  
  
  return(labels)
}


.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}

.dominateset <- function(xx,KK=20) {
  ###This function outputs the top KK neighbors.	
  
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
    
  }
  
  
  return(normalize(A))
}
# Calculate the mutual information between vectors x and y.
.mutualInformation <- function(x, y) {
  classx <- unique(x)
  classy <- unique(y)
  nx <- length(x)
  ncx <- length(classx)
  ncy <- length(classy)
  
  probxy <- matrix(NA, ncx, ncy)
  for (i in 1:ncx) {
    for (j in 1:ncy) {
      probxy[i, j] <- sum((x == classx[i]) & (y == classy[j])) / nx
    }
  }
  
  probx <- matrix(rowSums(probxy), ncx, ncy)
  proby <- matrix(colSums(probxy), ncx, ncy, byrow=TRUE)
  result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm=TRUE)
  return(result)
}

# Calculate the entropy of vector x.
.entropy <- function(x) {
  class <- unique(x)
  nx <- length(x)
  nc <- length(class)
  
  prob <- rep.int(NA, nc)
  for (i in 1:nc) {
    prob[i] <- sum(x == class[i])/nx
  }
  
  result <- -sum(prob * log(prob, 2))
  return(result)
}

.repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  if (is.null(dim(X))) {
    mx = length(X)
    nx = 1
  } else {
    mx = dim(X)[1]
    nx = dim(X)[2]
  }
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

normalize <- function(X) X / rowSums(X)


result<-function(Express=Express,group1,group2,group3,group4,group5,group6){
  collect <- matrix(rep(0,6*ncol(Express)),6,ncol(Express))
  collect[1,]=group1
  collect[2,]=group2
  collect[3,]=group3
  collect[4,]=group4
  collect[5,]=group5
  collect[6,]=group6
  #View(collect)
  tran <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y1=collect[1,]
  x1=t(y1)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y1[j]==x1[i])
        tran[i,j]=1
      else
        tran[i,j]=0
    }
  }
  
  tran2 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y2=collect[2,]
  x2=t(y2)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y2[j]==x2[i])
        tran2[i,j]=1
      else
        tran2[i,j]=0
    }
  }
  
  tran3 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y3=collect[3,]
  x3=t(y3)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y3[j]==x3[i])
        tran3[i,j]=1
      else
        tran3[i,j]=0
    }
  }
  
  tran4 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y4=collect[4,]
  x4=t(y4)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y4[j]==x4[i])
        tran4[i,j]=1
      else
        tran4[i,j]=0
    }
  }
  
  tran5 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y5=collect[5,]
  x5=t(y5)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y5[j]==x5[i])
        tran5[i,j]=1
      else
        tran5[i,j]=0
    }
  }
  
  tran6 <- matrix(rep(0,ncol(collect)*ncol(collect)),ncol(collect),ncol(collect))
  y6=collect[6,]
  x6=t(y6)
  
  for(i in 1:ncol(collect)){
    for(j in 1:ncol(collect)){
      if(y6[j]==x6[i])
        tran6[i,j]=1
      else
        tran6[i,j]=0
    }
  }
  
  
  
  trantotal=tran+tran2+tran3+(tran4+tran5+tran6)*0.5
  #trantotal=tran+tran2+tran3+(tran4+tran5+tran6)*0.4
  #View(trantotal)
  for (i in 1:ncol(collect)) {
    for (j in 1:ncol(collect)) {
      if(trantotal[i,j]>3){
        trantotal[i,j]=1
      }
      else{
        trantotal[i,j]=0
      }
    }
  }
  return(trantotal)
}



library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library("Seurat")
library("ggplot2")
library("dplyr")
library("hdf5r")



#cos的mat没有进行log2
hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5')
pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
#pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
mat<-as.matrix(pbmc[["RNA"]]@data)
a <- VariableFeatures(pbmc)
mat=mat[rownames(mat) %in% a,]
dim(mat)

# group1=read.csv("C:/non_graphST_label/Zero/e500_w05_4500_k4/151676_label_4000.csv",row.names = 1)
# group2=read.csv("C:/non_graphST_label/Zero/e500_w05_5000_k4/151676_label_4000.csv",row.names = 1)
# group3=read.csv("C:/non_graphST_label/Zero/e500_w05_4000_k4/151676_label_4000.csv",row.names = 1)
# group4=read.csv("C:/data_five_label1/151676/label1_151676_4500.csv",row.names = 1)
# group5=read.csv("C:/data_five_label1/151676/label1_151676_5000.csv",row.names = 1)
# group6=read.csv("C:/data_five_label1/151676/label1_151676_4000.csv",row.names = 1)


group1=read.csv("./label_4000.csv",row.names = 1)
group2=read.csv("./label_4500.csv",row.names = 1)
group3=read.csv("./label_5000.csv",row.names = 1)
group4=read.table("./label1.txt",row.names = 1)
group5=read.table("./label2.txt",row.names = 1)
group6=read.table("./label3.txt",row.names = 1)




group1=t(group1)
group2=t(group2)
group3=t(group3)
group4=t(group4)
group5=t(group5)
group6=t(group6)
mm=result(mat,group1,group2,group3,group4,group5,group6)


system.time(
  label<-spectralClustering(mm, K = 7)
  #kmeans(metadata(sce)$sc3$transformations$pearson_pca,7)
)

pre_label=label
pre_label[1] 
pre_label=as.data.frame(pre_label)
rownames(pre_label)=colnames(mat)
true_label=read.csv('/home/cuiyaxuan/spatialLIBD/151673/cluster_labels_151673.csv',row.names = 1)
a <- rownames(pre_label)
aa <- rownames(true_label)
pre_label=pre_label[rownames(pre_label) %in% aa,]

library("mclust")
true_label=as.array(true_label[,1])
ari=adjustedRandIndex(pre_label, true_label)



write.csv(ari,'ARI.csv')




