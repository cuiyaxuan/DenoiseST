
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

source('151673_cri4.R')
hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151672/151672_filtered_feature_bc_matrix.h5') #### to your path and project name
feature<-select_feature(hc1,4000,500)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
k=7 # k represent the number of spatial domains.
parLapply(cl,1:3,feature=feature,k=k,pearson_metric) 
stopCluster(cl)

tissue_local=read.csv("/home/cuiyaxuan/spatialLIBD/151672/spatial/tissue_positions_list.csv",row.names = 1,header = FALSE)
adj_matrix=construct_adj_matrix(feature[[1]],tissue_local)
write.table(adj_matrix,file="adj_matrix.txt",sep=" ",quote=TRUE)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
parLapply(cl,1:3,K=k,spectral_nei)
stopCluster(cl)



source('GNN_Tradition_6.R')

# hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151672/151672_filtered_feature_bc_matrix.h5') #### to your path and project name
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
write.csv(label,"label.csv")
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
write.csv(ari,"ARI.csv")

           ''')
