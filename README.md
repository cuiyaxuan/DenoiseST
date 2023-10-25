# DenoiseST
##### Due to the protocol issues of various space technology platforms, the data format is very different, and various platforms do not provide morphological images. For the convenience of users, we have changed the way of reading data to make it easier to use.<br>
##### DenoiseST is used on spatial transcriptomics (ST) datasets. In essence, you can refer to the following examples: <br>

##### * _DenoiseST on DLPFC from 10x Visium._ <br>
##### Using python virtual environment with conda. Please create a Pytorch environment, install Pytorch and some other packages, such as "numpy","pandas", "scikit-learn" and "scanpy". See the requirements.txt file for an overview of the packages in the environment we used to produce our results. <br>

#####  We execute the Nonlinear model in the python environment and can refer to the document DenoiseST_DP_run.py.  <br>
##### First, cd /home/.../DenoiseST-main/Full <br>

```python
from DenoiseST import DenoiseST
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp


# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
os.environ['R_HOME'] = '/scbio4/tools/R/R-4.0.3_openblas/R-4.0.3'
def setup_seed(seed=41):
    import torch
    import os
    import numpy as np
    import random
    torch.manual_seed(seed)  
    np.random.seed(seed)  # Numpy module.
    random.seed(seed)  # Python random module.
    if torch.cuda.is_available():
        # torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.cuda.manual_seed(seed)  
        torch.cuda.manual_seed_all(seed) 

setup_seed(41)

device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')


for i in [4000, 4500, 5000]:
   n_clusters = 7  ###### the number of spatial domains.
   file_fold = '/home/cuiyaxuan/spatialLIBD/151673' #### to your path
   adata = sc.read_visium(file_fold, count_file='151673_filtered_feature_bc_matrix.h5', load_images=True) #### project name
   adata.var_names_make_unique()
   model = DenoiseST(adata,device=device,n_top_genes=i)
   adata = model.train()
   radius = 50
   tool = 'mclust' # mclust, leiden, and louvain
   from utils import clustering

   if tool == 'mclust':
      clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)
   elif tool in ['leiden', 'louvain']:
      clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

   adata.obs['domain']
   adata.obs['domain'].to_csv(f"label_{i}.csv")

```
##### We execute the Linear model in the R environment and can refer to the document DenoiseST_TR_run.py. We need to install the R language (version R>4.0) in the system's default environment, and if you're using it for the first time, you'll need to install some R packages.
```python

import rpy2.robjects as robjects


robjects.r('''
install.packages('Seurat')
install.packages("hdf5r")
install.packages('ggplot2')
install.packages('dplyr')
install.packages('foreach')
install.packages('parallel')
install.packages('doParallel')
install.packages('mclust')

           ''')

```

##### We run R code within a Python environment.
```python

import rpy2.robjects as robjects


robjects.r('''

library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library("Seurat")
library("ggplot2")
library("dplyr")
library("hdf5r")
library(foreach)
library(parallel)
library(doParallel)

select_feature<-function(hc1,initial,rate){
  feature_dat=vector(mode = 'list', length = 3)
  for (i in 1:3) {
    pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = initial)
    all.genes <- rownames(pbmc)
    mat<-as.matrix(pbmc[["RNA"]]@data)
    a <- VariableFeatures(pbmc)
    mat=mat[rownames(mat) %in% a,]
    feature_dat[[i]] = mat
    initial=initial+rate
  }
  return(feature_dat)
}




pearson_metric<-function(X,feature,k){
  library(SingleCellExperiment)
  library(SC3)
  library(foreach)
  library(parallel)
  library(doParallel)
  pearson_pre <- function(mat,k){
    library(SingleCellExperiment)
    library(SC3)
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(mat),
        logcounts = log2(as.matrix(mat) + 1)
      )
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
    sce <- sc3_prepare(sce)
    sce <- sc3_calc_dists(sce)
    metadata(sce)$sc3$distances[[1]]=metadata(sce)$sc3$distances[[2]]
    metadata(sce)$sc3$distances[[3]]=metadata(sce)$sc3$distances[[2]]
    sce <- sc3_calc_transfs(sce)
    sce <- sc3_kmeans(sce, ks = k)
    sce <- sc3_calc_consens(sce)
    return(sce@metadata$sc3$consensus[[1]]$consensus)
  }
  
  if(X==1){
    pearson_dat1=pearson_pre(feature[[1]],k)
    write.table(pearson_dat1,file="data1.txt",sep=" ",quote=TRUE)
  }
  if(X==2){
    pearson_dat2=pearson_pre(feature[[2]],k)
    write.table(pearson_dat2,file="data2.txt",sep=" ",quote=TRUE)
  }
  if(X==3){
    pearson_dat3=pearson_pre(feature[[3]],k)
    write.table(pearson_dat3,file="data3.txt",sep=" ",quote=TRUE)
  }
  # if(X==4){
  #   pearson_dat4=pearson_pre(feature[[4]],k)
  #   write.table(pearson_dat4,file="data4.txt",sep=" ",quote=TRUE)
  # }
  # if(X==5){
  #   pearson_dat5=pearson_pre(feature[[5]],k)
  #   write.table(pearson_dat5,file="data5.txt",sep=" ",quote=TRUE)
  # }
}



#View(feature[[1]])

construct_adj_matrix<-function(feature,tissue_local){
  mat=t(feature)
  aa <- rownames(mat)
  tissue_local=tissue_local[rownames(tissue_local) %in% aa,]
  DF1 <- mutate(tissue_local, id = rownames(tissue_local))
  class(tissue_local)
  class(mat)
  mat=as.data.frame(mat)
  DF2 <- mutate(mat, id = rownames(mat))
  dat=merge(DF1,DF2,by="id")
  dat1=dat[,1:4]
  hyperG=matrix(0,dim(mat)[1],dim(mat)[1])
  for (i in 1:dim(mat)[1]) {
    x = dat1[i,3]
    y = dat1[i,4]
    for (j in 1:dim(mat)[1]) {
      x1 = dat1[j,3]
      y1 = dat1[j,4]
      radius=(x-x1)^2+(y-y1)^2
      if(radius<=16)
        hyperG[i,j]=1
    }
  }
  return(hyperG)
}






library(foreach)
library(parallel)
library(doParallel)

spectral_nei<-function(X,K){
  W1=read.table("./data1.txt",header = T,quote = "",sep=' ')
  W2=read.table("./data2.txt",header = T,quote = "",sep=' ')
  W3=read.table("./data3.txt",header = T,quote = "",sep=' ')
  # W4=read.table("./data4.txt",header = T,quote = "",sep=' ')
  # W5=read.table("./data5.txt",header = T,quote = "",sep=' ')
  hyperG=read.table("./adj_matrix.txt",header = T,quote = "",sep=' ')
  iter<-function(hyperG,W,K){
    hyperG=hyperG*0.9
    hyperG[hyperG==0]=0.1
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
    W=normalize(W)
    W = (W+t(W))/2
    hyperG=normalize(hyperG)
    W=as.matrix(hyperG) %*% as.matrix(W) %*% as.matrix(t(hyperG))
    W = (W+t(W))/2
    write.csv(W,"SC3_hyper_metrics.csv")
    label<-spectralClustering(W, K)
    return(label)
  }
  
  
  
  if(X==1){
    label1=iter(hyperG,W1,K)
    write.table(label1,file="label1.txt",sep=" ",quote=TRUE)
  }
  if(X==2){
    label2=iter(hyperG,W2,K)
    write.table(label2,file="label2.txt",sep=" ",quote=TRUE)
  }
  if(X==3){
    label3=iter(hyperG,W3,K)
    write.table(label3,file="label3.txt",sep=" ",quote=TRUE)
  }
  # if(X==4){
  #   label4=iter(hyperG,W4,K)
  #   write.table(label4,file="label4.txt",sep=" ",quote=TRUE)
  # }
  # if(X==5){
  #   label5=iter(hyperG,W5,K)
  #   write.table(label5,file="label5.txt",sep=" ",quote=TRUE)
  # }
}

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


hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
feature<-select_feature(hc1,4000,500)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
k=7 # k represent the number of spatial domains.
parLapply(cl,1:3,feature=feature,k=k,pearson_metric) 
stopCluster(cl)

tissue_local=read.csv("/home/cuiyaxuan/spatialLIBD/151673/spatial/tissue_positions_list.csv",row.names = 1,header = FALSE)
adj_matrix=construct_adj_matrix(feature[[1]],tissue_local)
write.table(adj_matrix,file="adj_matrix.txt",sep=" ",quote=TRUE)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
parLapply(cl,1:3,K=k,spectral_nei)
stopCluster(cl)


# hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
mat<-as.matrix(pbmc[["RNA"]]@data)
a <- VariableFeatures(pbmc)
mat=mat[rownames(mat) %in% a,]
# dim(mat)

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

label<-spectralClustering(mm, K=k)
pre_label=label
pre_label[1] 
pre_label=as.data.frame(pre_label)
rownames(pre_label)=colnames(mat)
true_label=read.csv('/home/cuiyaxuan/spatialLIBD/151673/cluster_labels_151673.csv',row.names = 1) # the ground truth labe l#### to your path
a <- rownames(pre_label)
aa <- rownames(true_label)
pre_label=pre_label[rownames(pre_label) %in% aa,]

library("mclust")
true_label=as.array(true_label[,1])
ari=adjustedRandIndex(pre_label, true_label)
paste("ARI:",ari)

           ''')

```








##### If the Python environment cannot run Linear model, you can also use the R language environment to execute it. It's necessary to set up a virtual R environment. DenoiseST requires R 4.0+ and Bioconductor 3.12+. Specific package dependencies are defined in the package DESCRIPTION and are managed by the Bioconductor and devtools installers. <br>
##### Using virtual environment with conda
```R

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SC3")


install.packages('Seurat')
install.packages("hdf5r")
install.packages('ggplot2')
install.packages('dplyr')
install.packages('foreach')
install.packages('parallel')
install.packages('doParallel')
install.packages('mclust')

```

##### Then, we execute the Linear model in the R environment


```R


library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library("Seurat")
library("ggplot2")
library("dplyr")
library("hdf5r")
library(foreach)
library(parallel)
library(doParallel)

source('151673_cri4.R')
hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
feature<-select_feature(hc1,4000,500)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
k=7 # k represent the number of spatial domains.
parLapply(cl,1:3,feature=feature,k=k,pearson_metric) 
stopCluster(cl)

tissue_local=read.csv("/home/cuiyaxuan/spatialLIBD/151673/spatial/tissue_positions_list.csv",row.names = 1,header = FALSE)
adj_matrix=construct_adj_matrix(feature[[1]],tissue_local)
write.table(adj_matrix,file="adj_matrix.txt",sep=" ",quote=TRUE)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
parLapply(cl,1:3,K=k,spectral_nei)
stopCluster(cl)



source('GNN_Tradition_6.R')

# hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
mat<-as.matrix(pbmc[["RNA"]]@data)
a <- VariableFeatures(pbmc)
mat=mat[rownames(mat) %in% a,]
# dim(mat)

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

label<-spectralClustering(mm, K=k)
pre_label=label
pre_label[1] 
pre_label=as.data.frame(pre_label)
rownames(pre_label)=colnames(mat)
true_label=read.csv('/home/cuiyaxuan/spatialLIBD/151673/cluster_labels_151673.csv',row.names = 1) # the ground truth label
a <- rownames(pre_label)
aa <- rownames(true_label)
pre_label=pre_label[rownames(pre_label) %in% aa,]

library("mclust")
true_label=as.array(true_label[,1])
ari=adjustedRandIndex(pre_label, true_label)
paste("ARI:",ari)

```

![image](https://github.com/cuiyaxuan/DenoiseST/blob/main/Image/151673pic.jpg)







# To prevent the algorithm from overfitting, we propose a simplified version. It not only reduces the complexity of the algorithm, but also reduces the operation time. <br>
## Simplified version <br>

```python
from DenoiseST import DenoiseST
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp

def setup_seed(seed=41):
    import torch
    import os
    import numpy as np
    import random
    torch.manual_seed(seed)  
    np.random.seed(seed)  # Numpy module.
    random.seed(seed)  # Python random module.
    if torch.cuda.is_available():
        # torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.cuda.manual_seed(seed)  
        torch.cuda.manual_seed_all(seed) 
        #os.environ['PYTHONHASHSEED'] = str(seed)

setup_seed(41)

device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')


n_clusters = 10  ###### the number of spatial domains.
file_fold = '/home/cuiyaxuan/spatialLIBD/3.Human_Breast_Cancer' #### to your path
adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True) #### project name
adata.var_names_make_unique()
model = DenoiseST(adata,device=device,n_top_genes=5000)
adata = model.train()
radius = 50
tool = 'mclust' # mclust, leiden, and louvain
from utils import clustering

if tool == 'mclust':
   clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)
elif tool in ['leiden', 'louvain']:
   clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

adata.obs['domain']
adata.obs['domain'].to_csv("label.csv")

```


## High resolution data <br>

```python
from DenoiseST import DenoiseST
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp

def setup_seed(seed=41):
    import torch
    import os
    import numpy as np
    import random
    torch.manual_seed(seed)  
    np.random.seed(seed)  # Numpy module.
    random.seed(seed)  # Python random module.
    if torch.cuda.is_available():
        # torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.cuda.manual_seed(seed)  
        torch.cuda.manual_seed_all(seed) 
        #os.environ['PYTHONHASHSEED'] = str(seed)

setup_seed(41)

device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')


n_clusters = 10  ###### the number of spatial domains.
file_path = '/home/cuiyaxuan/spatialLIBD/6.Mouse_Hippocampus_Tissue/' #please replace 'file_path' with the download path
adata = sc.read_h5ad(file_path + 'filtered_feature_bc_matrix_200115_08.h5ad') #### project name
adata.var_names_make_unique()
model = DenoiseST(adata,datatype='Slide',device=device,n_top_genes=4000)
adata = model.train()
radius = 50
tool = 'louvain' # mclust, leiden, and louvain
from utils import clustering

if tool == 'mclust':
   clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)
elif tool in ['leiden', 'louvain']:
   clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

adata.obs['domain']
adata.obs['domain'].to_csv("label.csv")

```


# Estimated number of spatial transcriptome data clusters

##### Using R virtual environment with conda <br>
```R

install.packages("devtools")
devtools::install_github("shaoqiangzhang/DEGman")


install.packages('Seurat')
install.packages("hdf5r")
install.packages('dplyr')
install.packages('ClusterR')

```
##### First, cd /home/.../DenoiseST-main/Full <br>

```R
library("Seurat")
library("dplyr")
library("hdf5r")
library("ClusterR")
source('spatial_data_estimate.R')
hc1= Read10X_h5('/Users/cyx/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5')
estimate_spatial(hc1=hc1)
```










# FVG:identify functionally variable genes

##### Using R virtual environment with conda <br>
```R

install.packages("devtools")
devtools::install_github("shaoqiangzhang/DEGman")


install.packages('Seurat')
install.packages("hdf5r")
install.packages('philentropy')
install.packages('dplyr')
install.packages('foreach')
install.packages('parallel')
install.packages('doParallel')

```

##### Then, we execute the FVG model in the R environment <br>
##### First, cd /home/.../DenoiseST-main/FVG <br>

```R
library(DEGman)
library("Seurat")
library("dplyr")
library("hdf5r")
library(philentropy)
library(foreach)
library(parallel)


library(doParallel)
source('distribution.R')

hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
label=read.csv("/home/cuiyaxuan/metric_change/revise_R2/est_151673/conlabel.csv",header = T,row.names = 1) # cluster label
k=2 ##### Define the cluster to be analyzed

dis<-distri(hc1,label,k)












#################################Spatial gene value compute################################



library(DEGman)
library("Seurat")
library("dplyr")
library("hdf5r")
library(philentropy)
library(foreach)
library(parallel)
library(doParallel)


files<-dir(path = "./",
             full.names = T,
             pattern = ".csv")
  library(tidyverse)
  df<-map(files,read.csv)
  class(df)
  #df1<-reduce(df,full_join)
  df1<-reduce(df,inner_join)
  df1=df1[-1,]
  
  write.csv(df1,"df1.csv")

source('test_finally.R')

tissue_local=read.csv("/home/cuiyaxuan/spatialLIBD/151673/spatial/tissue_positions_list.csv",row.names = 1,header = FALSE) #### to your path and project name
hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
pbmc=CreateSeuratObject(counts = hc1, project = "HC_1", min.cells = 10)
pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 30000)
all.genes <- rownames(pbmc)
mat<-as.matrix(pbmc[["RNA"]]@data)
max(mat)
min(mat)
dim(mat)
a <- VariableFeatures(pbmc)
mat=mat[rownames(mat) %in% a,]
dim(mat)
mat=t(mat)
aa <- rownames(mat)
tissue_local=tissue_local[rownames(tissue_local) %in% aa,]

DF1 <- mutate(tissue_local, id = rownames(tissue_local))
class(tissue_local)
class(mat)
mat=as.data.frame(mat)
DF2 <- mutate(mat, id = rownames(mat))
dat=merge(DF1,DF2,by="id")
#View(dat)
###对基因进行处理，合并基因和xy坐标
x_y_list=dat[,3:4]
#View(x_y_list)
dim(x_y_list)
dat=t(dat)
#View(dat)
df1=read.csv("df1.csv",row.names = 1,header = T)
class(df1[,1])
class(rownames(dat))
dim(dat)
class(rownames(dat))
#View(rownames(dat))
class(df1[,1])


dat1=dat[rownames(dat) %in% df1[,1],]
#geneval=as.numeric(dat[rownames(dat) %in% df1[,2][8],])
dat1=t(dat1)
dim(dat1)
############################################################################################
n_cores=24
cls <- makeCluster(n_cores) ## call 24 cpu cores
registerDoParallel(cls)
crinum=foreach(q=1:dim(dat1)[2],.combine='rbind') %dopar% cripar(q,dat1,x_y_list)
stopCluster(cls)
write.csv(crinum,"mark_gene_cri_vec.csv")
write.csv(df1[,1],"vectmark.csv")


```


