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

#############begin MVA method
cal_CNAs <- function(tumor_mat, normal_mat, ROW.name="ENSEMBLE_id", plot=TRUE){
library(cellrangerRkit)
library(biomaRt)   
matT <- tumor_mat
matN <- normal_mat
####combine matrix
share.rows <- intersect(rownames(matT), rownames(matN))
matT <- matT[which(rownames(matT) %in% share.rows),]
matN <- matN[which(rownames(matN) %in% share.rows),]
matN <- matN[order(match(rownames(matN), rownames(matT))),]

colnames(matT)  <- paste("T_", colnames(matT), sep="")
colnames(matN)  <- paste("N_", colnames(matN), sep="")
mat <- cbind(data.frame(as.matrix(matT)),data.frame(as.matrix(matN)))

## clean matrix to select genes for moving average
Z1 <- apply(mat[, which(substring(colnames(mat),1,1)=="T")],1, function(x)(sum(x>0)))/length(which(substring(colnames(mat),1,1)=="T"))
Z2 <- apply(mat[, which(substring(colnames(mat),1,1)=="N")],1, function(x)(sum(x>0)))/length(which(substring(colnames(mat),1,1)=="N"))
matx <- mat[Z1 > 0.4 | Z2 > 0.4, ]; dim(matx)

#annotation for gene coordinates, ###other data matrix can start from here
##if met errors, try to reinstall Biomart
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
if(ROW.name=="ENSEMBLE_id"){
anno <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name','start_position', 'end_position', 'band'),
      filters = 'ensembl_gene_id', values =rownames(matx), mart = ensembl)
anno <- anno[order(match(anno$ensembl_gene_id, rownames(matx))),]
data <- cbind(anno, matx)
}else if(ROW.name=="GENE_SYMBOL") {
  anno <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name','start_position', 'end_position', 'band'),
                filters = 'hgnc_symbol', values =rownames(matx), mart = ensembl) 
  anno <- anno[!duplicated(anno$hgnc_symbol),]
  matx <- matx[which(rownames(matx) %in% anno$hgnc_symbol),]
  anno <- anno[order(match(anno$hgnc_symbol, rownames(matx))),]
  data <- cbind(anno, matx)
}
##only keep genes in chr 1:Y
chrom <- c(1:22, "X", "Y")
datax <- data[which(data$chromosome_name %in% chrom),]
datax$chromosome_name[which(datax$chromosome_name=="X")] <- 23
datax$chromosome_name[which(datax$chromosome_name=="Y")] <- 24

#here add absolute chromosome length from hg20, ensemble 93
chr.L <- c(0, 248956422, 242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
           135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,
           58617616,64444167,46709983,50818468,156040895,57227415)

abspos <- NULL
for(i in 1:nrow(datax)){
   chr <- datax$chromosome_name[i]
   mid <- as.numeric(datax$start_position[i]) 
   if(chr==1){
     X <- mid
   } else {
     X <- mid+sum(chr.L[1:chr])
     }
  abspos <- c(abspos, X)
  i<- i+1
}

datax1 <- cbind(abspos, datax)
sam <- datax1[order(datax1$abspos),] ##order genes by abspos, use mat!
##FTT (Freeman and Tukey 1950) is found to be a better transformation to stabilize variance
norm.mat <- log(sqrt(sam[, 8:ncol(sam)])+ sqrt(sam[, 8:ncol(sam)]+1))
m <- apply(norm.mat, 2, mean)
norm.mat <- norm.mat-m

##calculate relative expression
mean.N <- apply(norm.mat[, which(substring(colnames(norm.mat), 1,1)=="N")], 1, median)
expr.relat <- as.matrix(norm.mat -mean.N); dim(expr.relat)

### adjust extreme values
expr.relat[expr.relat >=3]  <- 3 
expr.relat[expr.relat<= -3]  <- -3

#####################begin MVA and centered#####################
 MVA <- function(arr, s) {
	CNV <- arr
	n <- length(arr)
	for (i in 1:n){
		if(i<=s/2){
		sample <- arr[1:s]
		CNV[i] <- mean(sample)
		}
		else if(i >s/2 & i <=(n-s/2)){
				sample2 <- arr[(i-s/2):(i+s/2)]
				CNV[i] <- mean(sample2)
					}
					else {
						sample3 <- arr[(n-s+1):n]
						CNV[i] <- mean(sample3)
						}
						i <- i+1
			}
			
			return(CNV)		
	}
	
########################sample analysis#####################
print("start from chromosome ")
CNA.sam <- NULL
for (x in 1:max(as.numeric(datax$chromosome_name))){
	expr <- expr.relat[datax1$chromosome_name==x, ]
	print(x)
	print(nrow(expr))
    CNV.sam <- NULL
		i<-1
		for(i in 1:ncol(expr)) {
			sample <- expr[,i]
			n <- length(sample)
			step <- min(n,50)
		cnv <-MVA(sample,step)
		CNV.sam <- cbind(CNV.sam, cnv)
		i <- i+1
		}
			x <- x+1
	CNA.sam <- rbind(CNA.sam, CNV.sam)
}

count <- CNA.sam
#####normalize to genome wise express to get ratio values
mean.c <- apply(count, 2, mean)
log2ratio <- t(t(count)-mean.c)
colnames(log2ratio) <- colnames(expr.relat)
mean.cn <- apply(log2ratio, 2, mean)

final <- cbind(sam[, 1:7], log2ratio) ##remember abspos is sorted in mat, not in datax1
##assign a good rowname
final <- final[!duplicated(final$hgnc_symbol),]; dim(final)
rownames(final) <- final$hgnc_symbol  ##ready to output


log2ratio <- final[, 8:ncol(final)]

d <- dist(as.matrix(t(log2ratio)),method="euclidean")
fit <- hclust(d, method="ward.D")
group <- cutree(fit, k=2);CL <- data.frame(group);

###define tumor cells, normal cells
sdgroup1 <- sd(apply(log2ratio[, which(CL$group==1)],1,mean))
sdgroup2 <- sd(apply(log2ratio[, which(CL$group==2)],1,mean))
tumorid <- ifelse(sdgroup1>sdgroup2,1,2)
normalid <- ifelse(sdgroup1<sdgroup2,1,2)

class<- rownames(CL)  ##coloring clustering results
class[CL$group==tumorid] <- "defined_tumor"
class[CL$group==normalid] <- "defined_normal"

tissue<- rep("NA", length(rownames(CL))) ## tissue source
tissue[which(substring(rownames(CL),1,1) =="T")] <- "TumorSample"
tissue[which(substring(rownames(CL),1,1) =="N")] <- "NormalSample"

defMeta <- cbind(CL,tissue,class)
colnames(defMeta) <- c("cluster_id","defined_celltype","tissue")
results <- list(final, defMeta)
names(results) <- c("MVA_results","defined_celltype")

if(plot==TRUE){
  print('start ploting')
  library("gplots")
  library("devtools")
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  #source('/volumes/lab/ruligao/code/heatmap.3.R', chdir = TRUE) #alternative 
  
  library(RColorBrewer)
  my_palette <- colorRampPalette(rev(brewer.pal(n = 3, name = "RdBu")))(n = 999)
  col_breaks = c(seq(-1,-0.5001,length=50),seq(-0.5,-0.3,length=320),seq(-0.3,0.3,length=200),seq(0.3,0.5,length=380),seq(0.5001,1,length=50))

  rbPal <- colorRampPalette(brewer.pal(n = 6, name = "Set2")[1:2])
  hclust<- rbPal(2)[as.numeric(CL$group)]

  rbPal1 <- colorRampPalette(brewer.pal(n = 8, name = "Dark2")[1:2])
  tissuex<- rbPal1(2)[as.numeric(factor(tissue))]
  cells <- rbind(hclust,tissuex)

  chr <- as.numeric(final$chromosome_name) %% 2+1
  rbPal6 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal6(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)

  rang <- 0.5*(max(as.matrix(t(log2ratio)))-min(as.matrix(t(log2ratio))))
  
  if (ncol(norm.mat)< 3000){
    h <- 25
  } else {
    h <- 50
  }
  pdf("MVA50_CNA_heatmap.pdf", height=h, width=35)
  heatmap.3(as.matrix(t(log2ratio))/rang, dendrogram="r", distfun = dist, hclustfun = function(x) hclust(x, method="ward.D"),RowSideColors=cells, ColSideColors=chr1,Colv=NA, Rowv=T,notecol="black",col=my_palette,breaks=col_breaks,key=TRUE, keysize=1, density.info="none", trace="none", cexRow=0.4,cexCol=0.1,cex.main=1,cex.lab=0.1, symm=F,symkey=F,symbreaks=T,cex=1, margins=c(10,10))
  legend("topright", c("Normal","Tumor","def_norm","def_tumor"),col=c(brewer.pal(n = 8, name = "Dark2")[1:2],brewer.pal(n = 8, name = "Set2")[1:2]), pch=15, cex=2, bty='n')
  dev.off()
}

return(results)
}


