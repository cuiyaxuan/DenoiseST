# -*- coding: utf-8 -*

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
    torch.manual_seed(seed)  # 为CPU设置随机种子
    np.random.seed(seed)  # Numpy module.
    random.seed(seed)  # Python random module.
    if torch.cuda.is_available():
        # torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.cuda.manual_seed(seed)  # 为当前GPU设置随机种子
        torch.cuda.manual_seed_all(seed)  # 为所有GPU设置随机种子
        #os.environ['PYTHONHASHSEED'] = str(seed)



setup_seed(41)

device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# n_clusters = 7
# file_fold = '/home/cuiyaxuan/spatialLIBD/151673'
# adata = sc.read_visium(file_fold, count_file='151673_filtered_feature_bc_matrix.h5', load_images=True)
# adata.var_names_make_unique()
# run three times
for i in [4000, 4500, 5000]:
   n_clusters = 7
   file_fold = '/home/cuiyaxuan/spatialLIBD/151673'
   adata = sc.read_visium(file_fold, count_file='151673_filtered_feature_bc_matrix.h5', load_images=True)
   adata.var_names_make_unique()
   model = DenoiseST(adata,device=device,n_top_genes=i)
   adata = model.train()
   radius = 50
   tool = 'mclust' # mclust, leiden, and louvain
   from utils import clustering

   if tool == 'mclust':
      clustering(adata, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
   elif tool in ['leiden', 'louvain']:
      clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

   adata.obs['domain']
   adata.obs['domain'].to_csv(f"label_{i}.csv")
