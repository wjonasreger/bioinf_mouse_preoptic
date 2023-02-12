# bioinf_mouse_preoptic

Bioinformatics (STAT 530) final project to answer the question **"How does the mouse hypothalamic preoptic region work?"** using data.

## Data

* https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576
  * SOFT formatted family file(s): `raw/GSE113576_family.soft`
  * MINiML formatted family file(s): `raw/GSE113576_family.xml`
  * Series Matrix File(s): `raw/GSE113576_series_matrix.txt`
  * Barcodes: `raw/GSE113576_barcodes.tsv`
  * Genes: `raw/GSE113576_genes.tsv`
  * Matrix: `raw/GSE113576_matrix/`
* https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248
  * Data: `raw/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells/`
  
**Loading large data components stored on GitHub repository**
```r
source("lib/data_split_merge.R")

# load GSE113576 Matrix data
mat = matrixDataMerge("raw/GSE113576_matrix")

# load Merfish all cells data
df = csvDataMerge("raw/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells")
```

**Note: The above data files are raw. Any new processed data assets derived from these raw data, should be saved in the `data` directory and the script should be saved in the repository.**

