# Multiome QC preprocessing
A collection of scripts and notebooks that take single-cell multiome data from cellranger output to a QC-filtered and preprocessed dataset.

Before running the workflow:

1. Edit line 105 in `mymodule/utils.py` to replace the path to `RB_genes_human` file.
```python
RIBO_GENE_PATH = '/data/niecr/ramaiar/multiome_lucyrepo/analysis/mymodule/RB_genes_human'
```
2. Edit line 141, 170 and 204 in `multiome_qc_workflow.sh` to your `analysis` directory. (Make sure to include all the scripts and save them in the same pattern)
```python
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/
```
or 
```python
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/souporcell
```

### Run the workflow
Use `help` function to better understand the workflow inputs:
```sh
sh multiome_qc_workflow.sh --help
```
=======
# multiome_QC_preprocessing
A collection of scripts and notebooks to take multiple samples of single-cell multiome data from cellranger output to an anndata and archR project per sample with:
- A variety of calculated QC metrics.
- A filtered (default cellranger filtering) and unfiltered version of the sample.
- Initial RNA and ATAC embeddings of the filtered data.
