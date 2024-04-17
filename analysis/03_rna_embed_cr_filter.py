import anndata
import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import doubletdetection
import sys
import utils

def cellranger_filt_anndata(sample_dir,outs_dir):
    # for each sample
    for sample_name in os.listdir(sample_dir):
        sample_path = os.path.join(data, sample_name)
        if os.path.isdir(sample_path):
            print('Processing sample', sample_name)
            
            # get cellranger filtered barcodes
            filt_bc_path = os.path.join(sample_path, "filtered_feature_bc_matrix/barcodes.tsv.gz")
            cr_cells = pd.read_csv(filt_bc_path,sep='\t',compression='gzip',header=None)[0].tolist()
            cr_cells = [sample_name + '#' + x for x in cr_cells]
            
            # load unfiltered anndata and subset + embed
            prefilt_ad_path = os.path.join(outs_dir, sample_name,'pre_filt_rna_adata.h5ad')
            if os.path.exists(prefilt_ad_path):
                prefilt_ad = sc.read_h5ad(prefilt_ad_path)
                filt_ad = utils.subset_and_reprocess_rna(prefilt_ad,cr_cells)
                filt_ad.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Subset data by cellranger filtered cells and embed.')
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

    cellranger_filt_anndata(sample_dir,outs_dir)