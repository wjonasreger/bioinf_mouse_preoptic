# bioinf_mouse_preoptic

Bioinformatics (UIUC STAT 530) final project to answer the question **"How does the mouse hypothalamic preoptic region work?"** using data.

## Data

**Source**

* [GSE113576 Experiment Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576)
  * Barcodes: [GSE113576_barcodes.tsv](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fbarcodes%2Etsv%2Egz)
  * Genes: [GSE113576_genes.tsv](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fgenes%2Etsv%2Egz)
  * Matrix: [GSE113576_matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113576&format=file&file=GSE113576%5Fmatrix%2Emtx%2Egz)
* [Merfish Data](https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248)
  * Data: [Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv](https://datadryad.org/stash/downloads/file_stream/67671)
  
**Setting up large data folder**

* Create `large_data` folder in repository directory to store large data that will be ignored by GitHub. Store the above source data in `large_data/GSE113576` and `large_data/Merfish`, respectively.
