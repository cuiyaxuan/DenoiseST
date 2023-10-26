def dropout(adata):
    import matplotlib as mpl
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import warnings
    from anndata import AnnData
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams["font.sans-serif"] = "Arial"
    warnings.filterwarnings('ignore')

    adata.var_names_make_unique()
    adata
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=5000)
    adata=adata[:,adata.var['highly_variable']]
    adata.to_df()
    df=pd.DataFrame(adata.to_df())
    zero_count = np.count_nonzero(df.values == 0)
    rate = zero_count/(df.shape[0]*df.shape[1])
    print(rate)
    return rate