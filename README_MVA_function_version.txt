#***contents in README
## This is to calculate moving averages from single cell RNAseq data, basically application is to tell apart aneuploidy tumor cells from diploid normal cells
## Tt won't serve the purpose if a tumor is diploid
## by Ruli Gao, updated on Aug 22nd, 2018

#***input:
## tumor_mat, is the UMI count matrix from test sample, gene by cells
## normal_mat, is the UMI count matrix from normal control sample, gene by cells
## rownames, rownames of the input matrix, default "ENSEMBLE_id" or "GENE_SYMBOL"
## plot=TRUE, default, will generate a heatmap

#***output a list of two matrix objects
## "MVA_results", the MVA results for all cells 
## "defined_celltype", defined tumor or normal cells by clustering. 
## heatmap plot in PDF in working directory. I forced to have 2 clusters, always double check with the plot to make sure tumor cells have CNAs

#***algorithm
## genes applied are top genes that are detected in 40% of cells in either tumor or normal tissue, can tune down to 30% if you have much less gene coverage
## gene annotations is done with hg20, Biomart(hsapients, ensembl 93), to install Biomart, BiocInstaller::biocLite('grimbough/biomaRt')
## moving windows of 50 genes, ordered by chromosomal positions

#***To use, example run directly from 10X output
#####prepare count matrix
## > source("/volumes/lab/users/ruligao/code/update_MVA_function_version.R")
## > library(cellrangerRkit)
## > path1 <- "/Volumes/seq/projects/ATC_thyroid/LAI_p259T"   ##path one layer above out/ from the 10x output
## > tumor_obj <- load_cellranger_matrix(path1);
## > tumor_mat <- exprs(tumor_obj)

## > path2 <- "/Volumes/seq/projects/ATC_thyroid/LAI_p259N"
## > normal_obj <- load_cellranger_matrix(path2);
## > normal_mat <- exprs(normal_obj)

##run the method
###> CNA_result <- cal_CNAs(tumor_mat, normal_mat, ROW.name="ENSEMBLE_id", plot=TRUE)

