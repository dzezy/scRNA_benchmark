#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import time
import timeit

import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import celltypist as ct
import decoupler as dc

# conda install -c conda-forge scanpy python-igraph leidenalg
# conda install scikit-image
# conda install celltypist decoupler-py


# In[2]:


# scRNA analysis params

min_genes = 100 # per cell
min_cells = 3 # per gene
n_top_genes = 2000 # number of top highly variable genes to keep

cell_annot_model = "Immune_All_Low.pkl"  # for celltypist annotation


# In[3]:


## IMPORTANT: first argument is .h5 file path, second argument is sample name

# h5_path = 'melanoma_sample_filtered_feature_bc_matrix.h5'
# plot_path = 'melanoma_plots.pdf'
# log_path = 'melanoma_benchmarks.log'

try:
    if len(sys.argv) != 3:
        raise ValueError("Two string arguments required: input filepath and sample name.")

    h5_path = sys.argv[1]
    sample = sys.argv[2]
except ValueError as ve:
    print(f"Error: {ve}")
    print("Usage: python scanpy_analysis.py <input_path> <sample_name>")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    sys.exit(1)

plot_path = 'scanpy_' + sample + '_plots.pdf'
log_path = 'scanpy_' + sample + '_benchmarks.log'

benchmarks = ["Scanpy Benchmarks:\n"]  # list to hold benchmark strings
# initialize pdf file to save plots to
pdf = PdfPages(plot_path)


# In[4]:


# Define time benchmark decorator
def time_usage(func):
    def wrapper(*args, **kwargs):
        start = timeit.default_timer()
        retval = func(*args, **kwargs)
        elapsed = timeit.default_timer() - start
        
        benchmark_str = f"Function '{func.__name__}' took {elapsed:.2f} seconds\n"
        print(benchmark_str)

        # Save the result to benchmarks list
        benchmarks.append(benchmark_str)
        return retval
    return wrapper


# # Reading in Data

# In[5]:

@time_usage
def read_input_data(h5_path):
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    print("Reading in %s cells, %s genes." % (adata.X.shape[0], adata.X.shape[1]))
    return adata

adata = read_input_data(h5_path=h5_path)


# # QC

# In[6]:


@time_usage
def quality_control(adata, min_genes, min_cells):
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    # plot QC metrics for mt, ribo, hb genes 
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    # QC plots, save to pdf
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False
    )
    pdf.savefig()
    plt.close()


    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
    pdf.savefig()
    plt.close()


    # Hard Filtering cells and genes based on minimums
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    print("After QC and filtering, %s cells and %s genes remaining." % (adata.X.shape[0], adata.X.shape[1]))

    return adata


# In[7]:


adata = quality_control(adata, min_genes, min_cells)


# # Doublet Detection

# In[8]:


@time_usage
def doublet_detection(adata):
    # Predict doublets using nearest neighbor classifier of observed transcriptomes and simulated doublets
    sc.pp.scrublet(adata)
    # remove predicted doublets from dataset
    adata = adata[adata.obs.predicted_doublet == False]
    print("After doublet removal, %s cells and %s genes remaining." % (adata.X.shape[0], adata.X.shape[1]))
    
    return adata


# In[9]:


adata = doublet_detection(adata)


# # Normalization

# In[10]:


@time_usage
def normalization(adata):
    # Saving raw count data to separate layer
    adata.layers["counts"] = adata.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)

    return adata


# In[11]:


adata = normalization(adata)


# # Feature Selection

# In[12]:


@time_usage
def feature_selection(adata, n_top_genes):
    # select top 2000 most variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = n_top_genes)
    
    # plot highly variable genes
    sc.pl.highly_variable_genes(adata, show=False)
    pdf.savefig()
    plt.close()

    return adata


# In[13]:


adata = feature_selection(adata, n_top_genes)


# # Dimensionality Reduction

# In[14]:


@time_usage
def dim_reduction(adata):
    # do PCA reduction
    sc.tl.pca(adata)
    # plot PCA variance ratios
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True,show=False)
    pdf.savefig()
    plt.close()
    
    return adata


# In[15]:


adata = dim_reduction(adata)


# # NN graph construction and UMAP

# In[16]:


@time_usage
def umap(adata):
    # compute the neighborhood graph of cells using the PCA representation of the data matrix
    sc.pp.neighbors(adata)
    
    # embed graph in 2D for visualization with UMAP
    sc.tl.umap(adata)
    sc.pl.umap(adata,show=False)
    pdf.savefig()
    plt.close()
    
    return adata


# In[17]:


adata = umap(adata)


# # Leiden Clustering

# In[18]:


@time_usage
def clustering(adata):
    for res in [0.2, 0.5, 0.8]:
        sc.tl.leiden(
            adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
        )
    
    sc.pl.umap(
        adata,
        color=["leiden_res_0.20", "leiden_res_0.50", "leiden_res_0.80"],
        legend_loc="on data",
        show=False
    )
    pdf.savefig()
    plt.close()

    return adata


# In[19]:


adata = clustering(adata)


# # Cell Type Annotation

# In[20]:


@time_usage
def cell_type_annotation(adata, model):
    ## INITIAL ANNOTATION WITH CELLTYPIST
    # download appropriate cell annotation model
    ct.models.download_models(model=model, force_update=True)
    model = ct.models.Model.load(model=model)
    
    # pick preferred clustering resolution to annotate
    predictions = ct.annotate(adata, model=model, majority_voting=True, over_clustering=adata.obs['leiden_res_0.50'])
    
    # convert back to anndata
    adata = predictions.to_adata()


    ## DECOUPLER VERIFICATION
    # Query Omnipath and get PanglaoDB
    markers = dc.get_resource(name="PanglaoDB", organism="human")
    # Keep canonical cell type markers alone
    markers = markers[markers["canonical_marker"]]

    # Remove duplicated entries
    markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

    dc.run_mlm(mat=adata, net=markers, weight=None, source="cell_type", target="genesymbol", verbose=True, use_raw=False)
    acts = dc.get_acts(adata=adata, obsm_key="mlm_estimate")
    
    # plot cell annotations
    sc.pl.umap(
        acts,
        color=[
            "majority_voting",
            "B cells",
            "T cells",
            "Endothelial cells",
            "Fibroblasts",
            "Dendritic cells",
        ],
        wspace=0.5,
        ncols=3,
        show=False
    )
    pdf.savefig()
    plt.close()

    return adata


# In[21]:


adata = cell_type_annotation(adata, model=cell_annot_model)


# In[22]:


# close and save pdf file with plots
pdf.close() 

# write benchmarks to log file
with open(log_path, 'w') as f:
    for item in benchmarks:
        f.write(item)
print(f"Finished script!\nPlots saved to {plot_path} and log saved to {log_path}") 

