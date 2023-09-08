library("Seurat")
library("dplyr")
library("hdf5r")



###################################################################################

library(foreach)
library(parallel)
library(doParallel)

cripar<-function(q,dat1,x_y_list){
  # Create the function.
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  # Create the vector with numbers.
  geneval=as.numeric(dat1[,q])
  gene1=cbind(x_y_list,geneval)
  A=gene1[gene1[,3]>0,]
  b=round(A[,3],1)
  # Calculate the mode using the user function.
  if(length(b)<=50){
    cri=0
  }
  else if(length(b)<500 && length(b)>50){
    library(philentropy)
    cri_dat=A[,1:2]
    xx=distance(cri_dat, method="euclidean")
    sum(xx<=2)
    dim(A)[1]
    cri=(sum(xx<=2)-dim(A)[1])/2
  }
  else{
    result <- getmode(b)
    if(floor(result)==ceiling(result)){
      A=A[A[,3]>(floor(result)-0.5),]
      A=A[A[,3]<(ceiling(result)+0.5),]
    }
    else{
      A=A[A[,3]>floor(result),]
      A=A[A[,3]<ceiling(result),]
    }
    library(philentropy)
    cri_dat=A[,1:2]
    xx=distance(cri_dat, method="euclidean")
    sum(xx<=2)
    dim(A)[1]
    cri=(sum(xx<=2)-dim(A)[1])/2
  }
  return(cri)
}


############################################################################################


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

# dim(dat[rownames(dat) %in% df1[,1],])
# str(dat[rownames(dat) %in% df1[,2][2],])
# View(dat[rownames(dat) %in% df1[,2][2],])
dat1=dat[rownames(dat) %in% df1[,1],]
#geneval=as.numeric(dat[rownames(dat) %in% df1[,2][8],])
dat1=t(dat1)
dim(dat1)





############################################################################################


n_cores=24
cls <- makeCluster(n_cores)
registerDoParallel(cls)
crinum=foreach(q=1:dim(dat1)[2],.combine='rbind') %dopar% cripar(q,dat1,x_y_list)
stopCluster(cls)
write.csv(crinum,"mark_gene_cri_vec.csv")
write.csv(df1[,1],"vectmark.csv")