#!/usr/bin/R

#############
# Monocle 3 #
#############

# This script performs trajectory analysis of single-cell data using Monocle 3
# If it's a Seurat object
# If you wonder why DDRTree is no longer an option, it is not as fast nor robus;
# check: https://github.com/cole-trapnell-lab/monocle3/issues/187

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
library(monocle3)
library(Seurat)

### Functions
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(
  cds,
  time_bin = 1
){
  time_bin_i <- names(which.max(table(cds[[time_bin]])))
  cell_ids <- which(colData(cds)[, time_bin] == time_bin_i)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  tvar <- as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
  igraph::V(principal_graph(cds)[["UMAP"]])$name[tvar]
}

### Input
redu = c("UMAP", "tSNE", "PCA", "LSI", "Aligned")[1]
outdir = NULL
suffix = NULL
use_seurat = TRUE

# CD4T6
edataf = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2n3_sng_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
npcs = 38
cclust = "RNA_snn_res.0.6"
selectss = list(c("orig.virus2", "CV")); suffix = "_cv"
hvgf = '/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng_cv_25p/filtering/hvGenes.csv'

### Working folder
if(is.null(outdir)) outdir <- dirnamen(edataf, 3)
setwdc(paste0(sub('clustering', 'trajectory', outdir), "_", npcs, "PCs", suffix))
cat("Working at:", getwd(), "\n")

### Reading
edata <- readfile(edataf)

### Operations
scells = getsubset(selectss, edata@meta.data, v = TRUE)
edata <- edata[, scells]
myuse_genes <- if(!is.null(hvgf)){
  hvgdat <- read.csv(hvgf, stringsAsFactors = FALSE, row.names = 1)
  sum(hvgdat[, 'variable']); tvar <- rownames(hvgdat[hvgdat[, 'variable'], ])
  cat(commas(tvar), '\n')
  tvar
}else if(isTRUE(use_seurat) && casefold(class(edata)) == "seurat"){
  VariableFeatures(edata)
}else if(is.numeric(use_seurat) && casefold(class(edata)) == "seurat"){
  edata <- FindVariableFeatures(object = edata, selection.method = "vst", nfeatures = use_seurat, verbose = TRUE)
  VariableFeatures(edata)
}else{
  NULL
}

## The first step in working with Monocle 3 is to load up your data into Monocle 3's main class, cell_data_set:
gene_annotation <- edata@assays$RNA@meta.features
gene_annotation <- cbind(gene_short_name = rownames(gene_annotation), gene_annotation)
cds <- new_cell_data_set(
  expression_data = edata@assays$RNA@counts,
  cell_metadata = edata@meta.data,
  gene_metadata = gene_annotation
)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(
  cds = cds,
  method = "PCA",
  num_dim = npcs,
  norm_method = "log",
  use_genes = myuse_genes,
  residual_model_formula_str = "~nCount_RNA+percent.mt",
  alignment_group = NULL,
  pseudo_count = NULL,
  scaling = TRUE,
  verbose = TRUE
)
cds[[cclust[1]]] <- factormix(cds[[cclust[1]]])

if(1){
  cat("Total number of genes:", nrow(cds), "\n");
  cat("Highly variable genes:\n"); str(cds@preprocess_aux$gene_loadings)
}

p <- plot_pc_variance_explained(cds)
pdf("1_variance_explained_elbow.pdf")
print(p)
dev.off()
if(isTRUE(use_seurat)){
  cds@reducedDims@listData[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings
  p <- plot_pc_variance_explained(cds)
  pdf("1_variance_explained_elbow_seurat.pdf")
  print(p)
  dev.off()
  cds@reducedDims@listData[["PCA"]]  <- edata@reductions[["pca"]]@cell.embeddings[, 1:npcs]
}

## Step 2: Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "batch")
cds <- align_cds(
  cds = cds,
  preprocess_method = "PCA",
  alignment_group = NULL,
  alignment_k = 20,
  residual_model_formula_str = "~nCount_RNA+percent.mt",
  verbose = TRUE
)

## Step 3: Reduce the dimensions using UMAP
if(isTRUE(use_seurat)){
  cds@reducedDims@listData[["UMAP"]] <- edata@reductions[["umap"]]@cell.embeddings
}else{
  cds <- reduce_dimension(
    cds = cds,
    max_components = 2,
    reduction_method = redu,
    preprocess_method = NULL,
    umap.fast_sgd = TRUE,
    verbose = TRUE
  )
  if(redu != "UMAP"){
    cds@reducedDims@listData$UMAP <- cds@reducedDims@listData$tSNE
  }
}

## Step 4: Cluster the cells
# Monocle is able to learn when cells should be placed in the same trajectory as
# opposed to separate trajectories through its clustering procedure. Recall that
# we run cluster_cells(), each cell is assigned not only to a cluster but also to
# a partition. When you are learning trajectories, each partition will eventually
# become a separate trajectory
cds <- cluster_cells(
  cds = cds,
  reduction_method = "UMAP", # sounds like it only accepts UMAP
  k = 20,
  cluster_method = c("leiden", "louvain")[1],
  num_iter = 2,
  partition_qval = 0.05,
  weight = FALSE,
  resolution = NULL,
  random_seed = NULL,
  verbose = TRUE
)

## Step 5: Learn a graph
cds <- learn_graph(
  cds = cds,
  use_partition = TRUE,
  close_loop = TRUE,
  learn_graph_control = NULL,
  verbose = TRUE
)

p <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
pdf("2_partitions.pdf")
print(p)
dev.off()
table(cds@clusters@listData[["UMAP"]]$partitions)
if(!'partitions_bk' %in%  names(cds@clusters@listData$UMAP))
  cds@clusters@listData$UMAP$partitions_bk <- cds@clusters@listData$UMAP$partitions
tvar <- factor(rep('1', ncol(cds))); names(tvar) <- names(cds@clusters@listData$UMAP$partitions)
cds@clusters@listData$UMAP$partitions <- tvar

## Step 5: Learn a graph - one partition
cds <- learn_graph(
  cds = cds,
  use_partition = TRUE,
  close_loop = TRUE,
  learn_graph_control = NULL,
  verbose = TRUE
)
p <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
pdf("2_partitions2.pdf")
print(p)
dev.off()
table(cds@clusters@listData$UMAP$partitions)

pdf("3_select_root.pdf", width = 9, height = 7)
plot_cells(
  cds = cds,
  reduction_method = redu,
  color_cells_by = cclust[1],
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5
)
dev.off()

## Step 6: Order cells
# root_clust <- names(which.max(table(cds[[cclust[1]]])))
# root_cell <- rownames(colData(cds)[as.character(colData(cds)[[cclust[1]]]) == root_clust, ])[1]
cds <- order_cells(
  cds = cds,
  reduction_method = "UMAP",
  root_pr_nodes = get_earliest_principal_node(cds),
  root_cells = NULL,#root_cell,
  verbose = TRUE
)

p <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5
)
pdf("4_pseudotime.pdf")
print(p)
dev.off()

for(cname in getfound(cclust, colnames(colData(cds)), v = TRUE)){
  p <- plot_cells(
    cds = cds,
    reduction_method = "UMAP",
    color_cells_by = cname,
    group_cells_by = 'cluster',
    label_cell_groups = !TRUE,
    label_groups_by_cluster = !TRUE,
    # group_label_size = 6,
    # labels_per_group = 1,
    label_branch_points = TRUE, # black circles
    label_roots = TRUE, # white circles
    label_leaves = !TRUE, # gray circles
    graph_label_size = 3
  )
  pdf(paste0("5_", casefold(redu), "_", cname, ".pdf"), width = 7.5)
  print(p)
  dev.off()
}
graphics.off()

save(cds, file = "object.rdata")

## 7. Finding modules of co-regulated genes
# The function graph_test() uses a statistic from spatial autocorrelation
# analysis called Moran's I, which Cao & Spielmann et al showed to be effective
# in finding genes that vary in single-cell RNA-seq datasets.
# 'principal_graph' is recommended for trajectory analysis
pr_graph_test_res <- graph_test(
  cds = cds,
  neighbor_graph = c("knn", "principal_graph")[2],
  reduction_method = "UMAP",
  k = 25,
  method = c("Moran_I"),
  alternative = "greater",
  expression_family = "quasipoisson",
  cores = 1,
  verbose = TRUE
)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(
  cds = cds[pr_deg_ids,],
  resolution = 1e-2
)
cell_group_df <- tibble::tibble(
  cell = row.names(colData(neurons_cds)),
  cell_group = partitions(cds)[colnames(cds)]
)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
pheatmap::pheatmap(
  mat = agg_mat,
  cluster_rows = TRUE, cluster_cols = TRUE,
  scale = "column", clustering_method = "ward.D2",
  fontsize = 6
)
