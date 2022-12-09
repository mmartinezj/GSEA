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
if ("ComplexHeatmap" %in% row.names(installed.packages()) == FALSE) BiocManager::install("ComplexHeatmap")
if ("msigdbr" %in% row.names(installed.packages()) == FALSE) BiocManager::install("msigdbr")
if ("stringr" %in% row.names(installed.packages()) == FALSE) BiocManager::install("stringr")


suppressPackageStartupMessages({
  library(BiocManager, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(UpSetR, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(ggupset, quietly = TRUE)
  library(RColorBrewer, quietly = TRUE)
  library(optparse, quietly = TRUE)
  orgdb <- "org.Hs.eg.db"
  library(orgdb, quietly = TRUE, character.only = TRUE)
  library(dplyr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  library(fgsea, quietly = TRUE)
  library(reshape2, quietly = TRUE)
  library(ComplexHeatmap, quietly = TRUE)
  library(circlize, quietly = TRUE)
  library(msigdbr, quietly = TRUE)
  library(data.table, quietly = TRUE)
  library(DT, quietly = TRUE)
  library(stringr, quietly = TRUE)
})


##GET PARAMETERS AND DATA
##########################
#Paths
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
input <- "C:/Users/CBM/Desktop/BWH_counts/definitive_results/tsv/deseq2_noNA_padj_Affected_vs_Unaffected.tsv"
# file no NA -> deseq2_noNA_padj_Affected_vs_Unaffected
# file all -> all_genes_Affected_vs_Unaffected.tsv

data <- read.table(input, sep= "\t", quote = "", header=T, row.names = 1)

prefix <- "GSEA_results"

category <- "H"

subcategory <- NULL

genes <- "noNA_genes"

statistic <- "log+pvalue" #select from stat, log2foldchange, log+pvalue

#Outputs
resD0 <- paste0("results/GSEA_", statistic, "_")
resD <- gsub(':','_', paste0(resD0, category,'_', subcategory, genes, '/'))

if (!file.exists(resD)){
  dir.create(file.path(resD))
}
    
#Rank data using the desired metric
if (statistic == "log+pvalue"){
  data$fcsign <- sign(data$log2FoldChange)
  data$logP = -log10(data$pvalue)
  data$metric = data$logP/data$fcsign
  dat <- data$metric
  names(dat) <- as.character(rownames(data))
  dat_filtered <- dat[!duplicated(names(dat))] #remove rows with duplicate names = removes 7650 entries with NA as gene symbol
  dat_sort <- sort(dat_filtered, decreasing=TRUE)
  cat("Using log+pvalue metric")
} else{ 
    dat <- data$stat 
    names(dat) <- as.character(rownames(data))
    dat_filtered <- dat[!duplicated(names(dat))] #remove rows with duplicate names = removes 7650 entries with NA as gene symbol
    dat_sort <- sort(dat_filtered, decreasing=TRUE)
    cat("Using stat metric")
}


#Obtain hallmark gene sets relevant to Homo sapiens
mm_hallmark_sets <- msigdbr(species = "Homo sapiens", category = category) %>% 
  dplyr::select(gs_name, ensembl_gene)
head(mm_hallmark_sets) #Check using sets


#Calculate GSEA and write tables of results
set.seed(123)
egs <- GSEA(geneList = y_sort, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", seed = T, TERM2GENE = mm_hallmark_sets)
head(egs@result)
egs_df <- data.frame(egs@result)
egs_df_excel2 <- egs_df[, 2:length(egs_df)]
#egs_df_excel <- egs_df
#names(egs_df_excel)[names(egs_df_excel) == 'ID'] <- 'Identifier'
 
write.table(egs_df_excel2, file = paste(resD, "tableGSEA_0.05_", statistic, "_ENSEMBL_", category, "_", subcategory, genes, "_",prefix,".txt", sep =""), sep= "\t", quote = F, row.names = F)

if (category == "C5"){
  #Write table for GOBP
  C5_GOBP <- egs_df_excel2[grep("^GOBP", egs_df_excel2$Description),]
  #C5_GOBP <- egs_df[str_detect(egs_df$Description, "GOBP"),] #Alternative way to extract subsets
  write.table(C5_GOBP, file = paste(resD, "tableGSEA_0.05_", statistic, "_ENSEMBL_", category, "_", subcategory, genes, "_C5_GOBP_",prefix,".txt", sep =""), sep= "\t", quote = F, row.names = F)

  #Write table for GOCC
  C5_GOCC <- egs_df_excel2[grep("^GOCC", egs_df_excel2$Description),]
  write.table(C5_GOCC, file = paste(resD, "tableGSEA_0.05_", statistic, "_ENSEMBL_", category, "_", subcategory, genes, "_C5_GOCC_",prefix,".txt", sep =""), sep= "\t", quote = F, row.names = F)

  #Write table for GOMF
  C5_GOMF <- egs_df_excel2[grep("^GOMF", egs_df_excel2$Description),]
  write.table(C5_GOMF, file = paste(resD, "tableGSEA_0.05_", statistic, "_ENSEMBL_", category, "_", subcategory, genes, "_C5_GOMF_",prefix,".txt", sep =""), sep= "\t", quote = F, row.names = F)

  #Write table for HP
  C5_HP <- egs_df_excel2[grep("^HP", egs_df_excel2$Description),]
  write.table(C5_HP, file = paste(resD, "tableGSEA_0.05_", statistic, "_ENSEMBL_", category, "_", subcategory, genes, "_HP_",prefix,".txt", sep =""), sep= "\t", quote = F, row.names = F)
}

#PLOTS
##Dotplot / BarPlot

jpeg(file = paste(resD, prefix, "_dotplot.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
  par(mar = c(2, 2, 2, 5)) 
  dotplot(egs, x = "GeneRatio", color = "pvalue", showCategory = 20, font.size = 15)
invisible(dev.off())

##Gene-concept network

jpeg(file = paste(resD, prefix, "_gene_concept_net_", statistic, "_.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    cnetplot(egs, categorySize="pvalue", foldChange=NULL, font.size = 15, colorEdge = T)
invisible(dev.off())

##Ridgeline plot

jpeg(file = paste(resD, prefix, "_ridge_", statistic, "_.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    par(mar = c(2, 2, 2, 5)) 
    ridgeplot(egs, fill="p.adjust", core_enrichment = TRUE, orderBy = "NES")
invisible(dev.off())

##Heatplot

jpeg(file = paste(resD, prefix, "_heatplot_", statistic, "_.jpeg", sep =""), units = 'in', width = 20, height = 10, res = 300)
    heatplot(egs, foldChange=NULL)
invisible(dev.off())


##GSEAplot
jpeg(file = paste(resD, prefix, "_gsea_plot_", statistic, "_.jpeg", sep =""), units = 'in', width = 20, height = 10, res = 300)
    gseaplot(egs, geneSetID = 1)
invisible(dev.off())

##Upset plot (of the 20 first terms)

genes_20first <- as.data.frame(as.factor(head(egs@result$core_enrichment, 20)))
lista_20first <- list()
for (i in 1:nrow(genes_20first)){
    lista_20first[[i]] <- unlist(strsplit(as.character(genes_20first[i,1]),split="/"))   
}
uniq_genes <- as.character(unique(names(dat_sort)))
keep <- !is.na(uniq_genes)
uniq_genes <- uniq_genes[keep] #if there is a NA, row.names(mat_20first) <- uniq_genes, wont work

func_20first <- egs$Description[1:10]
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

jpeg(file = paste(resD, prefix, "_upset_20first_", statistic, "_.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    upset(mat_20first, nsets=10, order.by="freq", sets.bar.color="skyblue")
invisible(dev.off())

jpeg(file = paste(resD, prefix, "_upset.jpeg", sep =""), units = 'in', width = 15, height = 10, res = 300)
    enrichplot::upsetplot(egs)
invisible(dev.off())

p2 <- emapplot(egs, showCategory = 10)
cowplot::plot_grid(p2, ncol = 1, lables = LETTERS[1])
##################################

###Plot the GSEA if terms are provided and if not, plot the first 5 more abundant terms
sig_categories <- nrow(egs_df)
if (sig_categories <= 20){
  for (j in 1:sig_categories){
    pl <- gseaplot2(egs, geneSetID=j, title = egs$Description[j], base_size=40, color="red")
    desc <- gsub(" ", "_", egs$Description[j], fixed = TRUE) 
    filename <- paste(resD, desc, "_", prefix, ".jpeg", sep ="")
    ggsave(pl, file=filename, device = "jpeg", units= "in", height = 15, width = 20)
  }
  gseap <- gseaplot2(egs, geneSetID = 1:sig_categories, pvalue_table = TRUE)
  ggsave(gseap, file= paste(resD, prefix, "_gseaplot.jpeg", sep =""), device = "jpeg", units= "in", height = 20, width = 25)
  
} else{
    for (j in 1:19){
      pl <- gseaplot2(egs, geneSetID=j, title = egs$Description[j], base_size=40, color="red")
      desc <- gsub(" ", "_", egs$Description[j], fixed = TRUE) 
      filename <- paste(resD, desc, "_", prefix, ".jpeg", sep ="")
      ggsave(pl, file=filename, device = "jpeg", units= "in", height = 15, width = 20)
    }
    gseap <- gseaplot2(egs, geneSetID = 1:19, pvalue_table = TRUE)
    ggsave(gseap, file= paste(resD, prefix, "_gseaplot.jpeg", sep =""), device = "jpeg", units= "in", height = 20, width = 25)
}


######TESTING

