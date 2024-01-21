# DenoiseST
##### Due to the protocol issues of various space technology platforms, the data format is very different, and various platforms do not provide morphological images. For the convenience of users, we have changed the way of reading data to make it easier to use.<br>
##### DenoiseST is used on spatial transcriptomics (ST) datasets. In essence, you can refer to the following examples: <br>

##### * _DenoiseST on DLPFC from 10x Visium._ <br>
##### Using python virtual environment with conda. Please create a Pytorch environment, install Pytorch and some other packages, such as "numpy","pandas", "scikit-learn" and "scanpy". See the requirements.txt file for an overview of the packages in the environment we used to produce our results. <br>

#### Estimated number of spatial transcriptome data clusters. We utilized the ClusterR package in the R language for estimation; however, it can also be executed in a Python environment with the prerequisite installation of specific R packages
##### to install some R packages. <br>

```R
robjects.r('''
install.packages('Seurat')
install.packages("hdf5r")
install.packages('dplyr')
install.packages('ClusterR')
           ''')
```
##### First, cd /home/.../DenoiseST-main/Full <br>

```R


import rpy2.robjects as robjects


robjects.r('''
library("Seurat")
library("dplyr")
library("hdf5r")
library("ClusterR")
source('spatial_data_estimate.R')
hc1= Read10X_h5('/Users/cyx/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5')
estimate_spatial(hc1=hc1)
           ''')


```


##### Compute dropout rate. If the dropout rate is less than 0.85, use the simplified version; otherwise, use the full version. <br>

```python
import matplotlib as mpl
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import dropout
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
warnings.filterwarnings('ignore')
file_fold = '/home/cuiyaxuan/spatialLIBD/151672/' # your path
adata = sc.read_visium(file_fold, count_file='151672_filtered_feature_bc_matrix.h5', load_images=True)
drop=dropout.dropout(adata)
dropout.dropout(adata)
```


##### Full Version. We execute the Nonlinear model in the python environment and can refer to the document DenoiseST_DP_run.py.  <br>
##### First, cd /home/.../DenoiseST-main/Full <br>

```python
from DenoiseST import DenoiseST
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp


# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path. Please install mclust package <install.packages('mclust')> 
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

source('label_ARI.R')
true_label=read.csv('/home/cuiyaxuan/spatialLIBD/151673/cluster_labels_151673.csv',row.names = 1)
conlabel(hc1,k,true_label)
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

source('label_ARI.R')
true_label=read.csv('/home/cuiyaxuan/spatialLIBD/151673/cluster_labels_151673.csv',row.names = 1)
conlabel(hc1,k,true_label)

```


##### Visualization data
```python
import matplotlib as mpl
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import visual
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
warnings.filterwarnings('ignore')
file_fold = '/home/cuiyaxuan/spatialLIBD/151672/' # your path
adata = sc.read_visium(file_fold, count_file='151672_filtered_feature_bc_matrix.h5', load_images=True)
df_label=pd.read_csv('/home/cuiyaxuan/dropout/label.csv', index_col=0) 
visual.visual(adata,df_label)

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
a <- VariableFeatures(pbmc)
mat=mat[rownames(mat) %in% a,]
mat=t(mat)
aa <- rownames(mat)
tissue_local=tissue_local[rownames(tissue_local) %in% aa,]

DF1 <- mutate(tissue_local, id = rownames(tissue_local))
mat=as.data.frame(mat)
DF2 <- mutate(mat, id = rownames(mat))
dat=merge(DF1,DF2,by="id")
x_y_list=dat[,3:4]
dat=t(dat)
df1=read.csv("df1.csv",row.names = 1,header = T)

dat1=dat[rownames(dat) %in% df1[,1],]
dat1=t(dat1)
############################################################################################
n_cores=24
cls <- makeCluster(n_cores) ## call 24 cpu cores
registerDoParallel(cls)
crinum=foreach(q=1:dim(dat1)[2],.combine='rbind') %dopar% cripar(q,dat1,x_y_list)
stopCluster(cls)
write.csv(crinum,"mark_gene_cri_vec.csv")
write.csv(df1[,1],"vectmark.csv")


```


