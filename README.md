# Multiome QC preprocessing
A collection of scripts and notebooks that take single-cell multiome data from cellranger output to a QC-filtered and preprocessed dataset.

Before running the workflow:

1. Edit line 105 in `mymodule/utils.py`
```python
RIBO_GENE_PATH = '/data/niecr/ramaiar/multiome_lucyrepo/analysis/mymodule/RB_genes_human'
```
