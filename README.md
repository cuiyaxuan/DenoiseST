# DenoiseST
##### Due to the protocol issues of various space technology platforms, the data format is very different, and various platforms do not provide morphological images. For the convenience of users, we have changed the way of reading data to make it easier to use.<br>
##### DenoiseST is used on spatial transcriptomics (ST) datasets. In essence, you can refer to the following examples: <br>

##### * _DenoiseST on DLPFC from 10x Visium._ <br>
##### Using python virtual environment with conda
```python

anndata==0.8.0
backports.zoneinfo==0.2.1
brotlipy==0.7.0
certifi @ file:///croot/certifi_1665076670883/work/certifi
cffi @ file:///tmp/build/80754af9/cffi_1669355737025/work
charset-normalizer @ file:///tmp/build/80754af9/charset-normalizer_1630003229654/work
contourpy==1.0.6
cryptography @ file:///croot/cryptography_1665612644927/work
cycler==0.11.0
fonttools==4.38.0
geomloss==0.2.6
h5py==3.7.0
idna @ file:///croot/idna_1666125576474/work
igraph==0.10.5
importlib-metadata==5.1.0
Jinja2==3.1.2
joblib==1.2.0
kiwisolver==1.4.4
leidenalg==0.10.0
llvmlite==0.39.1
louvain==0.8.0
MarkupSafe==2.1.2
matplotlib==3.6.2
mkl-fft==1.3.1
mkl-random @ file:///tmp/build/80754af9/mkl_random_1626186064646/work
mkl-service==2.4.0
natsort==8.2.0
networkx==2.8.8
numba==0.56.4
numpy @ file:///croot/numpy_and_numpy_base_1668593735768/work
packaging==21.3
pandas==1.5.2
patsy==0.5.3
Pillow==9.2.0
POT==0.9.0
pycparser @ file:///tmp/build/80754af9/pycparser_1636541352034/work
pynndescent==0.5.8
pyOpenSSL @ file:///opt/conda/conda-bld/pyopenssl_1643788558760/work
pyparsing==3.0.9
PySocks @ file:///tmp/build/80754af9/pysocks_1605305779399/work
python-dateutil==2.8.2
pytz==2022.6
pytz-deprecation-shim==0.1.0.post0
requests @ file:///opt/conda/conda-bld/requests_1657734628632/work
rpy2==3.5.11
scanpy==1.9.1
scikit-learn==1.1.3
scikit-misc==0.1.4
scipy==1.9.3
seaborn==0.12.1
session-info==1.0.0
six @ file:///tmp/build/80754af9/six_1644875935023/work
statsmodels==0.13.5
stdlib-list==0.8.0
texttable==1.6.7
threadpoolctl==3.1.0
torch==1.13.0
torchaudio==0.13.0
torchvision==0.14.0
tqdm==4.64.1
typing_extensions @ file:///tmp/abs_ben9emwtky/croots/recipe/typing_extensions_1659638822008/work
tzdata==2023.3
tzlocal==4.3
umap-learn==0.5.3
urllib3 @ file:///croot/urllib3_1666298941550/work
zipp==3.10.0

```



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

##### Using R virtual environment with conda
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

```


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
hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5')
feature<-select_feature(hc1,4000,500)
detectCores()
cl <- makeCluster(3) # call 3 cpu cores
parLapply(cl,1:3,feature=feature,k=7,pearson_metric)
stopCluster(cl)
tissue_local=read.csv("/home/cuiyaxuan/spatialLIBD/151673/spatial/tissue_positions_list.csv",row.names = 1,header = FALSE)
adj_matrix=construct_adj_matrix(feature[[1]],tissue_local)
write.table(adj_matrix,file="adj_matrix.txt",sep=" ",quote=TRUE)
detectCores()
cl <- makeCluster(3) # call 5 cpu cores
parLapply(cl,1:3,K=7,spectral_nei)
stopCluster(cl)



source('GNN_Tradition_6.R')

hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5')
pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
mat<-as.matrix(pbmc[["RNA"]]@data)
a <- VariableFeatures(pbmc)
mat=mat[rownames(mat) %in% a,]
dim(mat)

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

label<-spectralClustering(mm, K = 7)
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
print(ari)

```
