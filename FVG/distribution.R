library(DEGman)
library("Seurat")
library("dplyr")
library("hdf5r")
library(philentropy)




distri<-function(hc1,label,k){
  pbmc=CreateSeuratObject(counts = hc1, project = "HC_1", min.cells = 10)
  pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 100000)

  mat<-pbmc@assays$RNA@data
  mat=as.matrix(mat)
  dim(mat)
  dim(label)
  label=t(label)
  mat1=rbind(label,mat)
  mat1=t(mat1)
  mat1<-mat1[order(mat1[,1]),]
  matlist=list()

    

  for (i in 1:k) {
    xx=t(mat1[mat1[,1]==i,])
    matlist[[i]]=xx[-1,]
  }

  n=1
  amat=matlist[[n]]
  for (i in 1:k) {
    if (i == n) {
      next
    }
    amat=cbind(amat,matlist[[i]])
  }


  dim(amat)

  result=DEGman(amat,dim(matlist[[n]])[2],(dim(amat)[2]-dim(matlist[[n]])[2]))

  file_name <- paste("result", k, ".csv", sep = "")
  write.csv(result,file = file_name)
}


hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151507/151507_filtered_feature_bc_matrix.h5')
label=read.csv("/home/cuiyaxuan/metric_change/revise_R2/est_151507/conlabel.csv",header = T,row.names = 1)
k=5

dis<-distri(hc1,label,k)







