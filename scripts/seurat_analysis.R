
## ------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(pheatmap)
library(SingleR)
library(viridis)
library(DoubletFinder)

## ------------------------------------------------------------------------------------------------------------------------------
## CLI arguments:
args <- commandArgs(trailingOnly = TRUE)

check_args <- function(args) {
  if (length(args) != 2) {
    stop("Error: Exactly two arguments are required - the input .h5 path, followed by the sample name.")
  }
}

# Check the number of arguments
check_args(args)

h5_path <- args[1]
samplename <- args[2] 

## ------------------------------------------------------------------------------------------------------------------------------
# scRNA analysis parameters:
plots_path <- paste0("seurat_", samplename, "_plots.pdf")
log_path <- paste0("seurat_", samplename, "_benchmarks.log")

# create pdf file to hold plots
pdf(plots_path)

# create list to hold benchmarks
benchmarks <- list("Seurat Benchmarks:")

min.cells <- 3  # minimum number of unique cells that a gene should be found in to be kept for analysis
nFeature_RNA_lower_lim <- 200
nCount_RNA_lower_lim <- 800
percent_MT_upper_lim <- 5

nfeatures <- 2000 # number of top highly variable genes to analyze
num_pcs <- 15 # number of Principal Components to keep for analysis, depends on data
doublet_formation_rate <- 0.075  # for doublet prediction purposes

# ref for cell annotation
#ref <- celldex::MonacoImmuneData()  
ref <- celldex::HumanPrimaryCellAtlasData()

## ------------------------------------------------------------------------------------------------------------------------------
# Define timing decorator function

timing_decorator <- function(func) {
  func_name <- deparse(substitute(func)) # Extract function name
  function(...) {
    start_time <- Sys.time()
    result <- func(...)
    
    time_diff <- Sys.time() - start_time
    time_diff_seconds <- round(as.numeric(time_diff, units = "secs"), 2)
    benchmark_str <- paste0("Function '", func_name, "' took ", time_diff_seconds, ' seconds')
    benchmarks <<- append(benchmarks, benchmark_str)
    return (result)
  }
}
## ------------------------------------------------------------------------------------------------------------------------------

read_input_data <- function(h5_path, samplename) {
#now <- Sys.time()
  dgCMatrix <- Read10X_h5(h5_path)
  seurat_obj <- CreateSeuratObject(counts=dgCMatrix, min.cells = min.cells, project=samplename)
  return (seurat_obj)
}

## ------------------------------------------------------------------------------------------------------------------------------
timed_read_input_data <- timing_decorator(read_input_data)
seurat_obj <- timed_read_input_data(h5_path, samplename)

## ------------------------------------------------------------------------------------------------------------------------------
quality_control <- function(seurat_obj, nFeature_RNA_lower_lim, nCount_RNA_lower_lim, percent_MT_upper_lim) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  print(VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3))
  print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm'))
  
  seurat_obj<- subset(seurat_obj, subset = nFeature_RNA > nFeature_RNA_lower_lim & nCount_RNA > nCount_RNA_lower_lim & percent.mt < percent_MT_upper_lim)
  
  return (seurat_obj)
}

## ------------------------------------------------------------------------------------------------------------------------------
timed_quality_control <- timing_decorator(quality_control)
seurat_obj <- timed_quality_control(seurat_obj, nFeature_RNA_lower_lim, nCount_RNA_lower_lim, percent_MT_upper_lim)

## ------------------------------------------------------------------------------------------------------------------------------
normalization <- function(seurat_obj) {
  
  seurat_obj <- NormalizeData(seurat_obj) # default: LogNormalize, scale factor of 10000
  return (seurat_obj)
}

## ------------------------------------------------------------------------------------------------------------------------------
timed_normalization <- timing_decorator(normalization)
seurat_obj <- timed_normalization(seurat_obj)

## ------------------------------------------------------------------------------------------------------------------------------
feature_selection <- function(seurat_obj, nfeatures) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  
  # get top 10 variable features for plot label
  top10 <- head(VariableFeatures(seurat_obj), 10)

  # plot variable features
  plot1 <- VariableFeaturePlot(seurat_obj)
  print(LabelPoints(plot = plot1, points = top10, repel = TRUE))
  
  return (seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_feature_selection <- timing_decorator(feature_selection)
seurat_obj <- timed_feature_selection(seurat_obj, nfeatures)


## ------------------------------------------------------------------------------------------------------------------------------
scale_data <- function(seurat_obj) {
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  return (seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_scale_data <- timing_decorator(scale_data)
seurat_obj <- timed_scale_data(seurat_obj)


## ------------------------------------------------------------------------------------------------------------------------------
dim_reduction <- function(seurat_obj) {
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  
  print(DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE))
  
  print(ElbowPlot(seurat_obj)) # show PC variance ratios
  return (seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_dim_reduction <- timing_decorator(dim_reduction)
seurat_obj <- timed_dim_reduction(seurat_obj)


## ------------------------------------------------------------------------------------------------------------------------------
clustering <- function(seurat_obj, num_pcs) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_pcs)
  seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1,0.3, 0.5, 0.7, 1)) 
  # higher resolutions -> more clusters

  Idents(seurat_obj) <- "RNA_snn_res.0.7" # pick proper resolution for given data
  
  return (seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_clustering <- timing_decorator(clustering)
seurat_obj <- timed_clustering(seurat_obj, num_pcs)


## ------------------------------------------------------------------------------------------------------------------------------
umap <- function(seurat_obj, num_pcs) {
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_pcs)
  print(DimPlot(seurat_obj, reduction = "umap" ,label=TRUE))
  
  return (seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_umap <- timing_decorator(umap)
seurat_obj <- timed_umap(seurat_obj, num_pcs)


## ------------------------------------------------------------------------------------------------------------------------------
doublet_detection <- function(seurat_obj, num_pcs, doublet_formation_rate) {

  ## pK Identification (no ground-truth) 
  sweep.res.list_seurat <- paramSweep(seurat_obj, PCs = 1:num_pcs, sct = FALSE, num.cores=1)
  sweep.stats_seurat <- summarizeSweep(sweep.res.list_seurat, GT = FALSE)
  bcmvn_seurat <- find.pK(sweep.stats_seurat)

  print(ggplot(bcmvn_seurat, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line())
  
  pK <- bcmvn_seurat %>%  # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  ## Homotypic Doublet Proportion Estimate 
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # run doubletFinder 
  seurat_obj <- doubletFinder(seurat_obj, 
                                 PCs = 1:num_pcs, 
                                 pN = 0.25, # default recommended
                                 pK = pK, 
                                 nExp = nExp_poi.adj,
                                 reuse.pANN = FALSE, sct = FALSE)
  
  df_classifications <- paste0("DF.classifications_0.25_", pK, '_', nExp_poi.adj)
  
  # visualize doublets
  print(DimPlot(seurat_obj, reduction = 'umap', group.by = df_classifications))
  
  # number of singlets and doublets
  table(seurat_obj@meta.data[[df_classifications]])

  # filter out doublets  
  filtered_metadata <- seurat_obj@meta.data[seurat_obj@meta.data[[df_classifications]] != "Doublet", ]
  seurat_obj <- subset(seurat_obj, cells = rownames(filtered_metadata))
  
  return(seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_doublet_detection <- timing_decorator(doublet_detection)
seurat_obj <- timed_doublet_detection(seurat_obj, num_pcs, doublet_formation_rate)


## ------------------------------------------------------------------------------------------------------------------------------
cell_type_annotation <- function(seurat_obj, ref) {
  # get counts data from Seurat object
  seurat_counts <- GetAssayData(seurat_obj, layer = 'counts')
  
  # run SingleR automatic cell annotation
  pred <- SingleR(test = seurat_counts, ref = ref, labels = ref$label.main)
  
  # save cell type label predictions
  seurat_obj$singleR.labels <- pred$labels[match(rownames(seurat_obj@meta.data), rownames(pred))]
  
  # Plot cell type label predictions
  print(DimPlot(seurat_obj, reduction = 'umap', group.by = 'singleR.labels'))
  
  ## Annotation Diagnostics 
  print(plotDeltaDistribution(pred)) # Based on deltas across cells
  
  # comparing to unsupervised clustering
  tab <- table(Assigned=pred$labels, Clusters=seurat_obj$seurat_clusters)
  print(pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10)))

  return (seurat_obj)
}


## ------------------------------------------------------------------------------------------------------------------------------
timed_cell_type_annotation <- timing_decorator(cell_type_annotation)
seurat_obj <- timed_cell_type_annotation(seurat_obj, ref)


## ------------------------------------------------------------------------------------------------------------------------------
# to save Seurat object as RDS
#saveRDS(seurat_obj, file = "./seurat_seurat.rds")

writeLines(unlist(benchmarks), log_path)

dev.off()

print(paste("Finished script! Plots saved to", plots_path, "and log saved to", log_path))
