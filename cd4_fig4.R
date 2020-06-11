#!/usr/bin/R

############
# Figure 4 #
############

# This script will create plots for figure 4 from cd4 data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')
source('/mnt/BioHome/ciro/scripts/seurat/utilities.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/F4'
setwd(dirfig)

# Global variables
colsname <- "/home/ciro/covid19/info/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

### 4.C Trajectory ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(monocle3)
## Input
clustname = "/home/ciro/large/covid19/results/trajectory/CD4T6_R1n2_sng2_25p_cv/object.rdata"

## Reading
cds <- theObjectSavedIn(clustname)

## Operations
p <- plot_cells(
  cds = cds,
  reduction_method = "UMAP",
  color_cells_by = "RNA_snn_res.0.6",
  group_cells_by = 'cluster',
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_branch_points = TRUE, # black circles
  label_roots = TRUE, # white circles
  label_leaves = TRUE, # gray circles
  graph_label_size = 4
) + scale_y_reverse() + theme_cowplot() + theme(legend.title = element_blank())
pdf(paste0("f4b_trajectory_umap.pdf"), width = 8, height = 8)
print(p)
dev.off()
p2 <- shut_it(p, "text")
p2$layers <- p2$layers[1:2]
pdf(paste0("f4b_trajectory_umap_blank.pdf"), width = 8, height = 8)
print(p2)
dev.off()

#### 4.D Principal UMAP ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.2'

## Reading
mycells <- theObjectSavedIn(clustname)

## Operations
mycells@meta.data$Cluster <- factormix(mycells@meta.data[, nres])
p <- DimPlot(
  object = mycells,
  reduction = 'umap',
  group.by = "Cluster",
  label = TRUE
)
pdf('f4d_umap.pdf')
print(p)
dev.off()
pdf('f4d_umap_blank.pdf')
print(shut_it(p, lays = "ext"))
dev.off()

#### 2.E Cluster markers ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData"
marknames = "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/markers/16PCs_RNA_snn_res.0.2_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
nres = 'RNA_snn_res.0.2'
suffix = "markers_shared"
ntop = 200
identnames = c(
  "0" = "Treg", "1" = "CTL", "4" = "CTL", "2" = "TCM", "3" = "TFH", "5" = "TH17"
)

## Reading
mycells <- theObjectSavedIn(clustname)
mygenes <- readfile(marknames, stringsAsFactors = FALSE)

## Operations
mycells$Name = factor(identnames[as.character(mycells@meta.data[, nres])], levels = unique(identnames))
mycells$Cluster = factor(mycells@meta.data[, nres], levels = names(identnames))
table(mycells@meta.data[, nres], mycells$Name, useNA = 'always')
table(mycells@meta.data[, nres], mycells$Cluster, useNA = 'always')
mygenes$gene <- gsub("'", "", mygenes$gene_name)

tvar <- mygenes$Dpct > .1; table(tvar)
mygenes <- mygenes[tvar, ]
tvar <- mygenes$avg_logFC > .25; table(tvar)
mygenes <- mygenes[tvar, ]
mygenes$cluster <- factor(mygenes$cluster, names(identnames))
mygenes <- mygenes[order(mygenes$cluster), ]
topgenes <- get_top_n(x = mygenes, n = ntop)
suffy <- paste0(suffix, ifelse(nrow(topgenes) != nrow(mygenes), paste0("_top", ntop), ""))
genes <- gsub("'", "", topgenes$gene_name)
genes <- getfound(genes, rownames(mycells), v = TRUE)
genesl <- make_list(x = topgenes, colname = "cluster", col_objects = "gene_name")

### Heatmap
fname <- paste0("f4e_heatmap_mean_", suffy)
source('/home/ciro/scripts/functions/handy_functions.R')
source('/home/ciro/scripts/functions/pheatmapCorrection.R')
pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
x <- custom_heatmap(
  object = mycells,
  rnames = genes,
  orderby = "Cluster",
  use_mean = "Cluster",
  sample_it = c("Cluster", "-1000"),
  scale_row = TRUE,
  categorical_col = c("Cluster", "Name"),
  feature_order = TRUE,
  couls = NULL,
  hcouls = c('yellow', 'black', 'blue'),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE
)
graphics.off()
# Check if it's consistent
genesl[["3"]]

dgenes = c(
  'FOXP3', 'CCR8', 'TIGIT', 'IL1R2', 'NKG7', 'HOPX', 'GZMB', 'GNLY', 'CCR7', 'IL7R',
  'SESN3', 'TRAT1', 'CD200', 'GNG4', 'POU2AF1', 'IGFBP4', 'CCR6', 'CTSH', 'IL4I1', 'LGALS3'
)
# dgenes = c(
#   'FOXP3', 'CCR8', 'TIGIT', 'CTLA4', 'CCR6', 'CTSH', 'IL4I1', 'LGALS3',
#   'NKG7', 'HOPX', 'GZMB', 'GNLY', 'CCR7', 'IL7R',
#   'SESN3', 'TRAT1', 'CD200', 'GNG4', 'POU2AF1', 'IGFBP4'
# )

### Dot-plots
p <- DotPlot(
  object = mycells,
  features = dgenes,
  group.by = "Cluster",
  cols = c('#fff4ba', '#ff0000'),
  col.min = -1.5, col.max = 1.5
) + coord_flip() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold.italic"),
    # axis.text.x = element_text(face = "bold", hjust = 1, angle = 45),
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)
fname <- "f4f_curtain"
pdf(paste0(fname, "_alternative.pdf"), width = 10, height = 13)
print(p)
graphics.off()
pdf(paste0(fname, "_alternative_blank.pdf"), width = 10, height = 13)
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()


### 4.X Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData"
allconfigs <- list(
  vconfig = list(
    name = "fig4_treg",
    genes = c('FOXP3', 'CCR8', 'CTLA4'),
    sfilter = list(c('RNA_snn_res.0.2', 0:5)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 1, size = c(4, 8),
    order = c('0', '5', 'REST')
  )
)
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
allconfigs <- list(
  vconfig = list(
    name = "fig4_ctl",
    genes = c('CCL3', 'CCL4', 'CCL5', 'XCL1', 'XCL2'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 1, size = c(4, 8),
    order = c('4', '8', 'REST')
  )
)

### 4.X Expression plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# CD4T6
myset = "CD4T6"
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
genes <- c('PRF1', 'GZMB', 'CCL3', 'CCL4', 'CCL5', 'XCL1', 'XCL2')
# CD4T24
myset = "CD4T24"
clustname = "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData"
genes <- c("FOXP3")

## Reading
mycells <- theObjectSavedIn(clustname)
genes <- getfound(genes, rownames(mycells), v = TRUE)

## Operations
for(gg in genes){
  cat(gg, "\n")
  fname <- paste0("marker_", myset, "_", gg, ".pdf")
  p <- FeaturePlot(
    object = mycells,
    reduction = 'umap',
    features = gg
  ) + scale_colour_gradientn(colours = couls)
  pdf(fname)
  print(p)
  dev.off()
}

#### 2.X Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.2'
norm_type = "RNA_snn_res."
couls = c(Hospital = "#fb8071", "non-Hospital" = "#b3de69")
identnames <- LETTERS[1:6]
names(identnames) <- 0:5

## Reading
mycells <- theObjectSavedIn(clustname)

## Operations
dfplot <- FetchData(
  object = mycells,
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster = factormix(identnames[as.character(dfplot[, nres])])
dfplot$Cluster <- factormix(dfplot[, nres])
dfplot$hospital <- ifelse(dfplot$orig.hospital == "Yes", "Hospital", "non-Hospital")
dfplot$hospital <- factor(dfplot$hospital, c("Hospital", "non-Hospital"))
table(dfplot$hospital, useNA = 'always')
p <- plot_pct(x = dfplot, groups = c("hospital", "Cluster"), normalise = 1) +
  coord_flip() + scale_fill_manual(values = couls) + labs(fill = NULL)

pdf("f4h_pie_hospital.pdf", width = 5, height = 10)
print(p)
dev.off()
pdf("f4h_pie_hospital_blank.pdf", width = 5, height = 10)
print(shut_it(p + mytheme))
dev.off()

p <- plot_pct(x = dfplot, groups = c("Cluster", "hospital"), type = "donut", print_ptables = TRUE) +
  labs(fill = NULL) +
  blank_theme + theme(
    legend.position = "right",
    strip.text = element_text(face = 'bold', size = 10),
    axis.text.x = element_blank()
  )

pdf("f4h_donut_hospital.pdf", width = 5, height = 10)
print(p)
dev.off()
pdf("f4h_donut_hospital_blank.pdf", width = 5, height = 10)
print(shut_it(p + mytheme))
dev.off()

table(dfplot$orig.run)
scells <- sample_grp(annot = dfplot, cname = 'orig.run', v = TRUE)
p <- ggplot(dfplot[scells, ], aes(x = UMAP_1, y =  UMAP_2, color = Cluster)) +
  geom_point() +
  facet_wrap(facets = ~orig.run) +
  labs(color = NULL) +
  mytheme +
  SetLegendPointsGG()

pdf("f4sup_umap_runs.pdf", width = 15, height = 10)
print(p)
dev.off()

#### S4.B Volcano ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
padjthr = 0.05
fcthr = .25
showgenes = c('TIGIT', 'LAG3', 'HAVCR2', 'PDCD1', 'DUSP4', 'CD70', 'DOK5', 'PRF1', 'GZMB', 'ZBED2', 'ZBTB32')
resname = "/home/ciro/large/covid19/results/dgea/CD4T6_R1n2_sng2_25p/comprs/cluster6hospital/YesvsNo/results_YesvsNo_mastlog2cpm.csv"
prefix = "sf4b_"
trimmer = "c170"
