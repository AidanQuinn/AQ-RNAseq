# Analysis of Biomek RNAseq QC test data
library(tibble)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(gglabeller)

source("https://raw.githubusercontent.com/AidanQuinn/codebits/master/read_cols_counts.R")

setwd("/n/data1/dfci/pedonc/kadoch/Aidan/20230810_Biomek_RNAseq_QC/full_analysis/aligned/")
################################################################################
# Read in the data
count_files <- list.files(pattern = "ReadsPerGene.out.tab")
counts <- read_cols_counts(
  file_list = count_files, col_n = 2, row_names = 1, skip_rows = 4)

# clean up column names and generate column data
clean_col_name <- function(column_name){
  a = strsplit(column_name, split = "\\/")[[1]][1]
  return(strsplit(a, split = "\\_")[[1]][1])
}
colnames(counts) <- sapply(colnames(counts), clean_col_name)
colnames(counts)[1] <- "GeneID"

################################################################################
# Generate column data
# pull treatment condition from filename
get_condition <- function(column_name){
  a = strsplit(column_name, split = "\\/")[[1]][1]
  return(strsplit(a, split = "\\_")[[1]][4])
}
condition <- sapply(count_files, get_condition)
condition[1:24] <- condition[25:48]

# pull rep from filename
get_rep <- function(column_name){
  a = strsplit(column_name, split = "\\/")[[1]][1]
  return(strsplit(a, split = "\\_")[[1]][5])
}
repp <- sapply(count_files, get_rep)
repp[1:24] <- repp[25:48]

# pull cell line from filename
get_cl <- function(column_name){
  a = strsplit(column_name, split = "\\/")[[1]][1]
  return(strsplit(a, split = "\\_")[[1]][3])
}
cell_line <- sapply(count_files, get_cl)
cell_line[1:24] <- cell_line[25:48]

# build the dataframe
col_data <- data.frame(
  "ID" = colnames(counts[,-1]),
  "batch" = c(rep(2,24), rep(1,24)),
  "condition" = condition,
  "cell_line" = cell_line,
  "rep" = repp
)

# clean up rownames
rownames(col_data) <- col_data$ID

################################################################################
# Set up gene/row data
# remove non-detected genes
cutoff <- 10^1.5
row_sums <- apply(counts[,-1], 1, sum)
hist(log10(row_sums[row_sums>0]), breaks = 50)
abline(v = log10(cutoff), lty = 3, col = "red2", lwd = 2)

keepers <- row_sums > cutoff
counts <- counts[keepers,]

# OPTIONAL STEP
# remove ERCC spike-in genes
keepers <- counts[,1]

# remove version from ENS GeneID
remove_gene_ID_vers <- function(gid_vers){
  return(
    strsplit(gid_vers, split = "\\.")[[1]][1]
  )}

gids <- sapply(counts$GeneID, remove_gene_ID_vers)
counts$GeneID <- gids 

# annotate gene IDs
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
#filters[53,]
attributes = listAttributes(ensembl)
attributes[grep("ID", attributes[,2]),]
#attributes[c(24,8,1,2),]

bm <- getBM(
  attributes = c('external_gene_name', 'description', 'ensembl_gene_id'),
  filters = 'ensembl_gene_id',
  values = counts[,1], 
  mart = ensembl)

i <- match(counts[,1], bm[,3], nomatch = NA)
all(bm[i,3] == counts[,1], na.rm = T) # sanity check

row_data <- tibble(bm[i,], col_row=1:nrow(counts))

# replace genes without names to ENSGID
i <- which(row_data$external_gene_name == "")
row_data$external_gene_name[i] <- counts[i,1]

# replace genes not found in biomart with ENSGID
i <- which(row_data$external_gene_name %in% NA)
row_data$ensembl_gene_id[i] <- counts[i,1]
row_data$external_gene_name[i] <- counts[i,1]

# replace non unique names
i <- which(duplicated(row_data$external_gene_name))
row_data$external_gene_name[i] <- counts[i,1]

counts <- counts[,-1]

#rownames(counts) <- row_data$ensembl_gene_id_version
rownames(counts) <- row_data$external_gene_name

################################################################################
# build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = col_data, 
  rowData = row_data, 
  design = ~ cell_line + condition)

dds <- DESeq(dds)
#vsd <- vst(dds, blind = T)
vsd <- vst(dds)

# compute Euclidean distances among samples
sampleDists <- dist(t(assay(vsd)))

# plot distance heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$batch, vsd$cell_line, vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# plot PCA of the data
pcaData <- plotPCA(vsd, 
                   intgroup=c("batch", "cell_line", "condition"), 
                   returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$batch <- factor(pcaData$batch, labels = c("KS", "AQ"))
pca <- ggplot(pcaData,
              aes(PC1, PC2, color=batch, shape = "cell_line",
                  label=paste(batch, cell_line, condition, 
                              sep = ":"))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  theme_light()
pca

gglabeller(pca)

resultsNames(dds)







