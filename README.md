
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bioinf_mouse_preoptic

<!-- badges: start -->
<!-- badges: end -->

The goal of `bioinf_mouse_preoptic` is to store analysis documents for
the Bioinformatics (UIUC STAT 530 Sp23) final project to answer the
question:

> **“How does the mouse hypothalamic preoptic region work?”** using data
> from
> [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576)
> and
> [DRYAD](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248).

## Install Project Package

You can install the development version of the associated project
package `mousePreopticR` from
[GitHub](https://github.com/wjonasreger/mousePreopticR) with:

``` r
# install.packages("devtools")
devtools::install_github("wjonasreger/mousePreopticR")
```

## Download Project Data

**Source**

- [GSE113576 Experiment
  Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576)
  - Barcodes:
    [GSE113576_barcodes.tsv](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fbarcodes%2Etsv%2Egz)
  - Genes:
    [GSE113576_genes.tsv](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fgenes%2Etsv%2Egz)
  - Matrix:
    [GSE113576_matrix.mtx](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fmatrix%2Emtx%2Egz)
- [Merfish
  Data](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248)
  - Data:
    [Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv](https://datadryad.org/stash/downloads/file_stream/67671)

**Setting up large data folder**

- Create `large_data` folder in repository directory to store large data
  that will be ignored by GitHub. Store the above source data in
  `large_data/GSE113576` and `large_data/merfish`, respectively. Rename
  data files as shown below:
  - `/GSE113576_barcodes.tsv` –\> `/barcodes.tsv`
  - `/GSE113576_genes.tsv` –\> `/genes.tsv`
  - `/GSE113576_matrix.mtx` –\> `/matrix.mtx`
  - `/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv` –\>
    `/data.csv`

## Setup Analysis Workspace

### Import packages

``` r
## import packages
library(mousePreopticR)
library(Seurat)
library(monocle3)
library(SpatialExperiment)
```

### Load Large Data

``` r
## load local data as expression matrix
#     - requires "matrix.mtx", "genes.tsv", and "barcodes.tsv" in data directory
mhpr_matrix = readMPR(data_dir = "large_data/GSE113576", cell_column = 1, gene_column = 1)
mhpr_genes = read.table("large_data/GSE113576/genes.tsv", sep='\t', col.names = c("ensembl_id", "gene_short_name"), row.names = 1)
mhpr_barcodes = read.table("large_data/GSE113576/barcodes.tsv", sep='\t', col.names = c("barcode"), row.names = 1)
mhpr_merfish = read.csv("large_data/merfish/data.csv", header = TRUE)
```

### Gene Expression Matrix Data Objects

``` r
## Seurat object with mhpr_matrix (Seurat)
mhpr_seurat = CreateSeuratObject(counts = mhpr_matrix)
mhpr_seurat
#> An object of class Seurat 
#> 27998 features across 31299 samples within 1 assay 
#> Active assay: RNA (27998 features, 0 variable features)
```

``` r
## CellDataSet object with mhpr_matrix, mhpr_barcodes, and mhpr_genes (monocle3)
mhpr_cds = new_cell_data_set(expression_data = mhpr_matrix,
                             cell_metadata = mhpr_barcodes,
                             gene_metadata = mhpr_genes)
mhpr_cds
#> class: cell_data_set 
#> dim: 27998 31299 
#> metadata(1): cds_version
#> assays(1): counts
#> rownames(27998): ENSMUSG00000051951 ENSMUSG00000089699 ...
#>   ENSMUSG00000096730 ENSMUSG00000095742
#> rowData names(1): gene_short_name
#> colnames(31299): AAACCTGAGATGTGGC-1 AAACCTGCACACAGAG-1 ...
#>   TTTGTCATCGTGGGAA-6 TTTGTCATCTTTACAC-6
#> colData names(1): Size_Factor
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```

``` r
## SpatialExperiment object with mhpr_merfish (SpatialExperiment)
names = list(
  col = colnames(mhpr_merfish),
  row = mhpr_merfish$Cell_ID,
  meta = c("Cell_ID", "Centroid_X", "Centroid_Y")
)

subs = list(
  col = setdiff(names[["col"]][10:170], "Fos"), # ignore Fos column due to NAs
  row = names[["row"]]
)

mhpr_mf_mat = matrixify(data = mhpr_merfish, names, subs, transpose = TRUE)

# define features for SpatialExperiment object
# expression matrix
counts = mhpr_mf_mat[["matrix"]]
# row data
row_data = data.frame(gene_short_name = counts@Dimnames[[1]])
rownames(row_data) = row_data$gene_short_name
# column data
col_data = data.frame(barcode = counts@Dimnames[[2]])
rownames(col_data) = col_data$barcode
# spatial coordinates
spatial_coords = as.matrix(apply(as.data.frame(mhpr_mf_mat[["metadata"]][3:4]), 2, as.numeric))

# create SpatialExperiment
mhpr_spe = SpatialExperiment(
  assays = list(counts = counts),
  rowData = row_data,
  colData = col_data,
  spatialCoords = spatial_coords
)
mhpr_spe
#> class: SpatialExperiment 
#> dim: 160 1027848 
#> metadata(0):
#> assays(1): counts
#> rownames(160): Ace2 Adora2a ... Ucn3 Vgf
#> rowData names(1): gene_short_name
#> colnames(1027848): 6749ccb4-2ed1-4029-968f-820a287f43c8
#>   6cac74bd-4ea7-4701-8701-42563cc65eb8 ...
#>   180ae0ff-9817-48b9-8d34-b374bda6e316
#>   83463d3a-29c5-40c3-b762-ffc3b7a11ad3
#> colData names(2): barcode sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : Centroid_X Centroid_Y
#> imgData names(0):
```
