import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mode
import scanpy as sc
import sklearn
import warnings

sys.path.insert(0, "../")
#import scgpt as scg

# extra dependency for similarity search
try:
    import faiss

    faiss_imported = True
except ImportError:
    faiss_imported = False
    print(
        "faiss not installed! We highly recommend installing it for fast similarity search."
    )
    print("To install it, see https://github.com/facebookresearch/faiss/wiki/Installing-Faiss")

warnings.filterwarnings("ignore", category=ResourceWarning)


model_dir = Path("/mnt/pixstor/dbllab/suli/tools_related/scgpt_data_model/scgpt_human")
adata = sc.read_h5ad("/mnt/pixstor/dbllab/suli/tools_related/scgpt_data_model/demo_train.h5ad")
cell_type_key = "Celltype"
gene_col = "index"


wdir = '/mnt/pixstor/dbllab/suli/tools_related/scGPT/'
os.chdir( wdir )
sys.path.insert(0, wdir)
import scgpt as scg

ref_embed_adata = scg.tasks.embed_data(
    adata,
    model_dir,
    cell_type_key=cell_type_key,
    gene_col=gene_col,
    batch_size=64,
    return_new_adata=True,
)

# Optional step to visualize the reference dataset using the embeddings
sc.pp.neighbors(ref_embed_adata, use_rep="X")
sc.tl.umap(ref_embed_adata)
sc.pl.umap(ref_embed_adata, color=cell_type_key, frameon=False, wspace=0.4)
