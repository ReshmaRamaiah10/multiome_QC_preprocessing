import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import os
import anndata
from scipy.io import mmread

RIBO_GENE_PATH = '/data/peer/lpaddock/data/utils/RB_genes_human'

def ribo_content_plots(adata, output_dir):

    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))

    axes[0].hist(adata.obs['ribo_pct'], bins=100)
    axes[0].set_xlabel('% Ribo-content', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)

    axes[1].scatter(adata.obs['log1p_total_counts'], adata.obs['ribo_pct'])
    axes[1].set_xlabel('Log library size', fontsize=14)
    axes[1].set_ylabel('% Ribo-content', fontsize=14)

    axes[2].scatter(adata.obs['log1p_n_genes_by_counts'], adata.obs['ribo_pct'])
    axes[2].set_xlabel('Log num. genes per cell', fontsize=14)
    axes[2].set_ylabel('% Ribo-content', fontsize=14)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ribo_content.png'), dpi=150)

def mt_content_plots(adata, output_dir):

    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))

    axes[0].hist(adata.obs['mito_pct'], bins=100)
    axes[0].set_xlabel('% MT-content', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)

    axes[1].scatter(adata.obs['log1p_total_counts'], adata.obs['mito_pct'])
    axes[1].set_xlabel('Log library size', fontsize=14)
    axes[1].set_ylabel('% MT-content', fontsize=14)

    axes[2].scatter(adata.obs['log1p_n_genes_by_counts'], adata.obs['mito_pct'])
    axes[2].set_xlabel('Log num. genes per cell', fontsize=14)
    axes[2].set_ylabel('% MT-content', fontsize=14)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'MT_content_pre_filt.png'), dpi=150)

def genes_metadata_plots(n_cells_by_counts, output_dir):
    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))
    
    axes[0].hist(n_cells_by_counts, bins=100)
    axes[0].set_xlabel('Number of cells a gene is expressed in', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)
    axes[0].set_title('Histogram of number of cells each gene is expressed in', fontsize=14)

    axes[1].hist(np.log(n_cells_by_counts + 1), bins=100)
    axes[1].set_xlabel('Log - Number of cells a gene is expressed in', fontsize=14)
    axes[1].set_ylabel('Frequency', fontsize=14)
    axes[1].set_title('Histogram of log of number of cells each gene is expressed in', fontsize=14)

    axes[2].hist(np.log(n_cells_by_counts + 1), bins=100)
    axes[2].set_xlabel('Log - Number of cells a gene is expressed in', fontsize=14)
    axes[2].set_ylabel('Frequency, ylim trimmed', fontsize=14)
    axes[2].set_title('Histogram of log of number of cells each gene is expressed in', fontsize=14)
    axes[2].set_ylim([0, 1000])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'gene_metadata_pre_filt.png'), dpi=150)

def cells_genes_correlation_plots(log1p_total_counts, log1p_n_genes_by_counts, output_dir):
    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))
    
    sns.histplot(log1p_total_counts, bins=100, ax=axes[0])
    axes[0].set_xlabel('Log-library size', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)
    axes[0].set_title('Histogram of log library size', fontsize=14)

    sns.histplot(log1p_n_genes_by_counts, bins=100, ax=axes[1])
    axes[1].set_xlabel('Log of num. genes per cell', fontsize=14)
    axes[1].set_ylabel('Frequency', fontsize=14)
    axes[1].set_title('Histogram of number of genes expressed in each cell', fontsize=14)

    sns.scatterplot(x=log1p_total_counts, y=log1p_n_genes_by_counts, ax=axes[2])
    axes[2].set_ylabel('Log of num. genes per cell', fontsize=14)
    axes[2].set_xlabel('Log library size', fontsize=14)
    corr_coef = np.corrcoef(log1p_total_counts, log1p_n_genes_by_counts)[0, 1]
    axes[2].set_title(f'Correlation = {corr_coef:.3f}', fontsize=14)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells_genes_correlation_pre_filt.png'), dpi=150)

def load_anndata(mtx_dir):
    # load in mtx, barcodes and features
    matrix_path = os.path.join(mtx_dir, 'matrix.mtx.gz')
    matrix = mmread(matrix_path).tocsr()
    barcodes_path = os.path.join(mtx_dir, 'barcodes.tsv.gz')
    barcodes = pd.read_csv(barcodes_path,sep='\t',compression='gzip',header=None,
                          names = ['og_barcode'])
    features_path = os.path.join(mtx_dir, 'features.tsv.gz')
    features = pd.read_csv(features_path, sep='\t',compression='gzip',header=None, 
                           names=['ensemblID','name', 'feature_type', 'chr', 'start', 'end'])
    
    # select gene features
    gene_features = np.where(features['feature_type'] == 'Gene Expression')[0]
    
    # generate anndata
    var = features.iloc[gene_features]
    var.index = var['name'].tolist()
    obs = barcodes
    obs.index = obs['og_barcode'].tolist()
    rna_adata = anndata.AnnData(matrix[gene_features,:].T,
                    var = var,
                    obs = obs)

    return rna_adata
    
def prefilt_anndata(data, outs_dir):
    for sample_name in os.listdir(data):
        sample_path = os.path.join(data, sample_name)
        if os.path.isdir(sample_path):
            print('Processing sample', sample_name)
            mtx_dir = os.path.join(sample_path, "raw_feature_bc_matrix")
            if os.path.exists(mtx_dir):
                adata = load_anndata(mtx_dir)
                adata.var_names_make_unique()
                # rename barcodes to sample_name+barcodes and add the sample data to anndata
                adata.obs_names = [f"{sample_name}#{barcode}" for barcode in adata.obs_names]
                adata.obs['batch'] = sample_name

                # calculate qc metrics
                sc.pp.calculate_qc_metrics(adata, inplace = True)
                # ribo content
                ribo_genes = pd.read_csv(RIBO_GENE_PATH, header = None)[0].tolist()
                ribo_genes = adata.var_names[adata.var_names.isin(ribo_genes)]
                adata.var['ribo'] = adata.var_names.isin(ribo_genes)
                row_sum_adata_ribo = np.sum(adata[:, ribo_genes].X.toarray(), axis=1)
                adata.obs['ribo_pct'] = row_sum_adata_ribo / adata.obs['total_counts'] * 100
                # mito content
                mt_genes = adata.var_names[adata.var_names.str.startswith('MT-')]
                adata.var['mito'] = adata.var_names.isin(mt_genes)
                row_sum_adata_mito = np.sum(adata[:, mt_genes].X.toarray(), axis=1)
                adata.obs['mito_pct'] = row_sum_adata_mito / adata.obs['total_counts'] * 100

                
                # save as a new h5ad file
                sample_outs_dir = os.path.join(outs_dir, sample_name)
                os.makedirs(sample_outs_dir, exist_ok=True)
                adata.write_h5ad(os.path.join(sample_outs_dir, f"pre_filt_rna_adata.h5ad"))
                print("Saving modified adata as pre_filt_rna_adata.h5ad")

                # create a new directory to save the plots
                qc_plots_dir = os.path.join(sample_outs_dir, "QC_plots")
                os.makedirs(qc_plots_dir, exist_ok=True)

                # generate qc plots prior to filteration
                print('Generating pre-filtering QC plots')
                cells_genes_correlation_plots(adata.obs['log1p_total_counts'], adata.obs['log1p_n_genes_by_counts'], qc_plots_dir)
                genes_metadata_plots(adata.var['n_cells_by_counts'], qc_plots_dir)
                mt_content_plots(adata, qc_plots_dir)
                ribo_content_plots(adata, qc_plots_dir)
                print('Processed sample', sample_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process data for quality control.')
    parser.add_argument('--path', '-p', type=str, help='Path to the data directory.')
    parser.add_argument('--out', '-o', type=str, help='Path to the output directory.')
    args = parser.parse_args()

    if args.path:
        sample_dir = args.path
    else:
        sample_dir = os.getcwd()
        
    if args.out:
        outs_dir = args.out
    else:
        outs_dir = os.getcwd()

    prefilt_anndata(sample_dir, outs_dir)