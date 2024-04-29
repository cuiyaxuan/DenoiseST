# DenoiseST : A dual-channel unsupervised deep learning-based denoising method to identify spatial domains and functionally variable genes in spatial transcriptomics
![image](https://github.com/cuiyaxuan/DenoiseST/blob/master/Image/%E5%B9%BB%E7%81%AF%E7%89%871.png)
## Tip: To facilitate researchers' usage, examples of our project can be run in the full folder's IPython notebooks (after configuring the environment dependencies as described in the README). We will soon optimize the README page for better usability. Additionally, we are developing a web tutorial for DenoiseST, which will be released soon. We apologize for any inconvenience this may cause. <br>


##### Due to the protocol issues of various space technology platforms, the data format is very different, and various platforms do not provide morphological images. For the convenience of users, we have changed the way of reading data to make it easier to use.<br>

##### DenoiseST is used on spatial transcriptomics (ST) datasets. In essence, you can refer to the following examples: <br>

##### * _DenoiseST on DLPFC from 10x Visium._ <br>
##### Using python virtual environment with conda. Please create a Pytorch environment, install Pytorch and some other packages, such as "numpy","pandas", "scikit-learn" and "scanpy". See the requirements.txt file for an overview of the packages in the environment we used to produce our results. Alternatively, you can install the environment dependencies in the following sequence to minimize environment conflicts. <br>

```R
conda create -n pipeline
source activate pipeline

conda search r-base
conda install r-base=4.2.0
conda install python=3.8

conda install conda-forge::gmp
conda install conda-forge::r-seurat==4.4.0
conda install conda-forge::r-hdf5r
conda install bioconda::bioconductor-sc3

conda install conda-forge::pot
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
pip install scanpy
pip install anndata==0.8.0
pip install pandas==1.4.2

pip install rpy2==3.5.1
pip install scikit-learn==1.1.1
pip install scipy==1.8.1
pip install tqdm==4.64.0

```

##### to install some R packages. <br>
```R
import rpy2.robjects as robjects
robjects.r('''
install.packages('ClusterR')
install.packages('foreach')
install.packages('doParallel')
install.packages('mclust')
           ''')
```

#### Estimated number of spatial transcriptome data clusters. We utilized the ClusterR package in the R language for estimation; however, it can also be executed in a Python environment with the prerequisite installation of specific R packages
##### First, cd /home/.../DenoiseST-main/Full <br>

```R

import rpy2.robjects as robjects

robjects.r('''
library("Seurat")
library("dplyr")
library("hdf5r")
library("ClusterR")
source('spatial_data_estimate.R')
hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151672/151672_filtered_feature_bc_matrix.h5')  #### to your path and file name
estimate_spatial(hc1=hc1)
           ''')

```
##### We choose the cluster number where the first occurrence of the numerical value of f(k) reaches 1. It will automatically generate a picture depicting the estimated number of clusters in the current directory.<br>
![image](https://github.com/cuiyaxuan/DenoiseST/blob/master/Image/est.png)
##### For specific details, please refer to the file "1_Example_lowresolve_test.ipynb".Instructions for running other types of spatial transcriptomics data can be found in the IPython notebooks under the "full" directory.<br>
##### If your current system environment does not support both R and Python, you can refer to the User Manual for more detailed instructions. Alternatively, you can create two separate environments to run the program. <br>
