#!/usr/bin/env Rscript
##
##R Script to make GSEA analysis
##The input must be a tab-separated file with 2 columns: first column must contain the gene IDs
##and the second one the log2FC values.
##
##THE ORDER OF THE COLUMNS IS VERY IMPORTANT!!
## 
## Eva Sacristán and Sandra González (GENGS CBMSO)
##################################################################

#
#######################

##CLEAN AND LIBRARIES
#####################

cat("\n Checking if all libraries are installed \n")
#rm(list=ls())

repos = "http://cran.us.r-project.org"
if ("UpSetR" %in% row.names(installed.packages())  == FALSE) install.packages("UpSetR", repos = repos)
if ("ggplot2" %in% row.names(installed.packages())  == FALSE) install.packages("ggplot2", repos = repos)
if ("RColorBrewer" %in% row.names(installed.packages())  == FALSE) install.packages("RColorBrewer", repos = repos)
if ("optparse" %in% row.names(installed.packages())  == FALSE) install.packages("optparse", repos = repos)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if ("clusterProfiler" %in% row.names(installed.packages()) == FALSE) BiocManager::install("clusterProfiler")
if ("enrichplot" %in% row.names(installed.packages()) == FALSE) BiocManager::install("enrichplot")
if ("org.Hs.eg.db" %in% row.names(installed.packages()) == FALSE) BiocManager::install("org.Hs.eg.db")

suppressPackageStartupMessages({
  library(BiocManager, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(UpSetR, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(optparse, quietly = TRUE)
  orgdb <- "org.Hs.eg.db"
  library(orgdb, quietly = TRUE, character.only = TRUE)
})


##GET PARAMETERS AND DATA
##########################

input <- "../BWH_counts/definitive_results/tsv/0.05_sig_padj_Affected_vs_Unaffected.tsv"

data <- read.table(input, sep= "\t", quote = "", header=T)
#Solo carga 21948 observaciones en vez de las 28526 que tiene el excel ?????????
#Ya carga todos los datos, problema era que faltaba la opción quote = ""
prefix <- "ORA_results"

##GENERATE PLOT
#####################

cat("\n Calculation ORA and making some nice plots \n")
    
##Calculate GSEA

keep_FC <- abs(data$log2FoldChange) <= 4.5
data_filter <- data[keep_FC,]
dat_names <- data_filter$X

#Calculate ORA and write tables of results

ego2 <- enrichGO(gene         = dat_names,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

write.table(ego2, file = paste("ORA/tableGO_0.05_",prefix,".txt", sep =""), sep= "\t", quote = F)



#############

set.seed(123)
#cada vez sale un resultado diferente -> utilizar seed? -> SI
egs <- enrichGO(gene = dat_sort, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont= "ALL", pvalueCutoff=0.05)
egs_genename <- setReadable(egs, OrgDb = "org.Hs.eg.db")

write.table(egs, file = paste("ORA/tableGO_0.05_",prefix,".txt", sep =""), sep= "\t", quote = F)
write.table(egs_genename, file = paste("ORA/tableGO_0.05_",prefix,"_genename.txt", sep =""), sep= "\t", quote = F)

##Dotplot / BarPlot

jpeg(file = paste("ORA/", prefix, "_dotplot.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    dotplot(egs, x = "GeneRatio", color = "pvalue", showCategory = 20, font.size = 15)
invisible(dev.off())

##Gene-concept network

jpeg(file = paste("ORA/", prefix, "_gene_concept_net.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    cnetplot(egs_genename, categorySize="pvalue", foldChange=dat_sort, font.size = 15, colorEdge = T)
invisible(dev.off())

##Ridgeline plot

jpeg(file = paste("ORA/", prefix, "_ridge.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    ridgeplot(egs, fill="pvalue")
invisible(dev.off())

##Heatplot

jpeg(file = paste("ORA/", prefix, "_heatplot.jpeg", sep =""), units = 'in', width = 20, height = 10, res = 300)
    heatplot(egs, foldChange=dat_sort)
invisible(dev.off())

##GSEAplot
jpeg(file = paste("ORA/", prefix, "_heatplot.jpeg", sep =""), units = 'in', width = 20, height = 10, res = 300)
    gseaplot(egs, geneSetID = 1)
invisible(dev.off())

##Upset plot (of the 20 first terms)

genes_20first <- as.data.frame(as.factor(head(egs@result$core_enrichment, 20)))
lista_20first <- list()
for (i in 1:nrow(genes_20first)){
    lista_20first[[i]] <- unlist(strsplit(as.character(genes_20first[i,1]),split="/"))   
}
uniq_genes <- as.character(unique(names(dat_sort)))

func_20first <- egs$Description[1:20]
mat <- matrix(0L, nrow = length(uniq_genes), ncol = length(func_20first)) 

for (i in 1:length(uniq_genes)) {
for (j in 1:length(func_20first)) {
gen <- uniq_genes[i]
if (gen %in% lista_20first[[j]]) {
mat[i,j] =  1
}}} 
 
mat_20first <- as.data.frame(mat)
colnames(mat_20first) <- func_20first
row.names(mat_20first) <- uniq_genes

jpeg(file = paste("ORA/", prefix, "_upset_20first.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    upset(mat_20first, nsets=10, order.by="freq", sets.bar.color="skyblue")
invisible(dev.off())

##################################

###Plot the GSEA if terms are provided and if not, plot the first 5 more abundant terms

for (j in 1:5){
  pl <- gseaplot2(egs, geneSetID=j, title = egs$Description[j], base_size=40, color="red")
  desc <- gsub(" ", "_", egs$Description[j], fixed = TRUE) 
  filename <- paste("ORA/", desc, "_", prefix, ".jpeg", sep ="")
  ggsave(pl, file=filename, device = "jpeg", units= "in", height = 15, width = 20)
}

##QUIT
####################

cat("\n ALL DONE!!! ^w^ \n")

