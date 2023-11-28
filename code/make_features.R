#!/usr/bin/env Rscript
#------------------------------------------------------SETUP-----------------------------------------------------#

# Load libraries
library(tidyverse)
library(data.table)
library(BuenColors)
library(Seurat)
library(irlba)
library(Matrix)
library(future)
library(reticulate)
library(ggrastr)
library(tidytext)
library(matrixTests)
source("utils.R")

# Read commandline args
args <- commandArgs(trailingOnly = TRUE)

# Parameters
# Assumes all data is present in ../data
# name <- args[1L]
# cores <- as.integer(args[2L])
# number_pcs <- as.integer(args[3L])
# vargenes <- as.integer(args[4L])
# clus_res <- as.numeric(args[5L])
# inputData <- as.numeric(args[6L])
# inputAnnot <- as.numeric(args[7L])
# isProcessed <- as.numeric(args[8L])
# TODO: find a way to do something like argparse with list/nullable input so parameters can go to default like annotations that could be provided
name <- "human_heme_test2"
cores <- 6
number_pcs <- NULL #35
vargenes <- NULL #1500
clus_res <- NULL #0.6
inputData <- "/home/louwenjjr/Documents/opentargets-bi/pops_gene_features/data/human_heme/16populations_RNAcounts.txt"
inputAnnot <- NULL
isProcessed <- TRUE

# todo: write commands used to log file in code dir


# Set up parallelization
# Remember to use htop to delete forgotten forks
Sys.setenv(R_FUTURE_FORK_ENABLE = T)
options(future.globals.maxSize = 2048 * 1024^2)
plan(strategy = "multicore", workers = cores)

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
# file <- paste0("../data/", name, "/16populations_RNAcounts.txt")
if (dir.exists(inputData)) {
    filelist = Sys.glob(paste0(inputData, "/*gz"))  # make sure inputData dir only contains expression data files
    # Read all matrices into a list
    # assumes rownames to be present as first column in each file
    datalist <- lapply(filelist, read_sparse_mat(rowIdType = rowIdType))

    # Bind assuming same row order
    mat <- do.call("cbind", datalist)
} else {
    mat <- data.frame(fread(inputData, sep = "\t"))[,-1] %>%
        data.matrix() %>%
        Matrix(sparse = TRUE)
    rownames(mat) <- data.frame(fread(inputData), select = 1, skip = 1, sep = "\t")[,1]

    # Convert to ENSG, drop duplicates, and fill in missing genes
    mat <- ConvertToENSGAndProcessMatrix(mat, "human_symbol")
}

# read annotations, if any
if (inputAnnot != NULL) {
    mat.annot <- data.frame(fread(inputAnnot), row.names=1, header=T)
}

# colnames(mat) <-  gsub("[.]", "-", colnames(mat)) used for human_airway

# Convert to ENSG, drop duplicates, and fill in missing genes
mat <- ConvertToENSGAndProcessMatrix(mat, "human_symbol")

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

minFeatures <- 200  # todo: why is this sometimes used but not always?
so <- CreateSeuratObject(counts = mat, project = name, meta.data = mat.annot)

# Clean up
rm(mat)

# QC
if (isProcessed == FALSE){
    so <- subset(so,
                 subset = nFeature_RNA > quantile(so$nFeature_RNA, 0.05) &
                     nFeature_RNA < quantile(so$nFeature_RNA, 0.95))
}
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1000000)
so <- ScaleData(so)

# Identify variable genes
# todo: probably only when data is raw
so <- FindVariableFeatures(so, nfeatures = vargenes)
# Plot variable genes with and without labels
PlotAndSaveHVG(so, name)

# Run PCA todo: see if features param is needed
so <- RunPCA(so, npcs = 15, features = row.names(so))
# Project PCA to all genes
so <- ProjectDim(so, do.center = T)
# Plot Elbow
# todo: automate finding number of PCs based on elbow plot
PlotAndSavePCAElbow(so, 15, name)

# Run ICA
# todo: Wrap ICA in a try since it seems to fail for some studies (see human_heme where only projected_pcaloadings, average_expression and diffexprs_tstat_clusters_pre_def is outputted)
so <- RunICA(so, nics = number_pcs)
# Project ICA to all genes
so <- ProjectDim(so, reduction = "ica", do.center = T)

# Cluster cells
so <- FindNeighbors(so, dims = 1:number_pcs, nn.eps = 0)
so <- FindClusters(so, resolution = clus_res, n.start = 100)

# UMAP dim reduction
so <- RunUMAP(so, dims = 1:number_pcs, min.dist = 0.4, n.epochs = 500,
              n.neighbors = 10, learning.rate = 0.1, spread = 2)

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)
# Plot known clusters on UMAP (if applicable)
PlotAndSaveUMAPClusters(so, so@meta.data$CellType, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP
marker_genes <- c("CCDC67", "DEUP1", "FOXN4", "CDC20B", "RERGL", "MCAM", "PDGFRB", "ACTA2", "MYL9", "ASCL3", "CFTR", "FOXJ1", "MUC5AC", "SFTPA2", "CA2", "CAV1", "ANXA3", "CAV2")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "CellType"
clus <- unique(so@meta.data$CellType)
demarkers_pre_def <- WithinClusterFeatures(so, "CellType", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 30, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

