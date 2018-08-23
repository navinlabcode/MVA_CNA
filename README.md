# MVA_from_10XRNAseq


Description:
This is a straightforward method written in R to calculate copy numbers from 10X RNAseq data using moving averages of gene expression. Genome cordinates are annotated with hg20 and ensembl 93 through BioMart. Top genes that are expressed in 40% of either testing group of control group are selected to smoothing. The sliding windows are 50 genes.  This method is designed to separate aneuploid tumor cells from diploid normal cells within tumor tissues. By default it will not help if you meet a diploid tumor.

Get ready:
install {BioMart} using BiocInstaller::biocLite('grimbough/biomaRt')

To run:
> source("/volumes/lab/users/ruligao/code/update_MVA_function_version.R")

> CNA_result <- cal_CNAs(tumor_mat, normal_mat, ROW.name="ENSEMBLE_id", plot=TRUE)

To prepare input:
1) tumor_mat, is the UMI count matrix from testing sample, genes in the rows and cells in the columns
2) normal_mat, is the UMI count matrix from normal control sample, genes in the rows and cells in the columns
3) ROW.name, rownames of the input matrix, default "ENSEMBLE_id" or "GENE_SYMBOL"
4) plot=TRUE, default, will generate a heatmap of the copy numbers of all cells

Description of output:
1) it outputs a list of two matrix objects
2) "MVA_results", the MVA results for all cells 
3) "defined_celltype", defined tumor or normal cells by clustering
4) heatmap plot in PDF in working directory. always double check with the plot to make sure tumor cells have CNAs

An example run directly from 10X output:
> source("/volumes/lab/users/ruligao/code/update_MVA_function_version.R")

> library(cellrangerRkit)	

#prepare count matrix of testing sample

> path1 <- "/Volumes/seq/projects/ATC_thyroid/LAI_p259T"   ##path that is one layer above outs/ from the 10x output	

> tumor_obj <- load_cellranger_matrix(path1)	

> tumor_mat <- exprs(tumor_obj)	

#prepare count matrix of normal control sample	

> path2 <- "/Volumes/seq/projects/ATC_thyroid/LAI_p259N"

> normal_obj <- load_cellranger_matrix(path2)	

> normal_mat <- exprs(normal_obj)	

> CNA_result <- cal_CNAs(tumor_mat, normal_mat, ROW.name="ENSEMBLE_id", plot=TRUE)



An example run directly from Seurate object:
> source("/volumes/lab/users/ruligao/code/update_MVA_function_version.R")

> tumor_mat <- as.matrix(tumor_obj@raw.data)

> normal_mat <- as.matrix(normal_obj@raw.data)

> CNA_result <- cal_CNAs(tumor_mat, normal_mat, ROW.name="GENE_SYMBOL", plot=TRUE)

An example run directly from count matrix:
> source("/volumes/lab/users/ruligao/code/update_MVA_function_version.R")

> CNA_result <- cal_CNAs(tumor_mat, normal_mat, ROW.name="GENE_SYMBOL", plot=TRUE)
