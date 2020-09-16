#!/usr/bin/R

############
# Figure 5 #
############

# This script will create plots for figure 5 from cd4 data

source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_plotting.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_utilities.R')
library(Seurat)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/Figure_5'
setwdc(dirfig)
clustnamebk = "none"

# Global variables
colsname <- "../data/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

### 5.X IL2 comparison heatmap ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T24_seurat.rdata"
resnamef = "../data/IL1R2pvsIL1R2n.csv"
sselect = list(c('RNA_snn_res.0.6', '0'))
prefix = "f5_"

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
res <- readfile(resnamef, row.names = 1, stringsAsFactors = FALSE)

## Operations
res <- setrows(res)
genes <- sapply(
  X = c(TRUE, FALSE),
  FUN = function(u){
    getDEGenes(
      x = res,
      pv = 0.05,
      fc = .25,
      upreg = u,
      v = TRUE
    )
  }
)
genes <- unlist(lapply(genes, head, 100), use.names = FALSE)
res[headtail(genes), ]

suffy <- paste0(basename(dirname(resnamef)), "_", summary_subset(sselect))

scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
void <- add_gene_tag(
  lgenes = "IL1R2",
  annot = mycells@meta.data,
  mat = mycells@assays$RNA@data,
  v = TRUE
)
mycells@meta.data <- joindf(mycells@meta.data, void)
table(mycells@meta.data[scells, c('RNA_snn_res.0.6', "tag_IL1R2")])

### Heatmap
fname <- paste0(prefix, "heatmap_", suffy)
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R')
pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
custom_heatmap(
  object = mycells,
  cnames = scells,
  rnames = genes,
  orderby = c("tag_IL1R2", "pca"),
  # use_mean = "orig.donor",
  sample_it = c(cname = "tag_IL1R2", maxln = "-2000"),
  scale_row = TRUE,
  categorical_col = c('tag_IL1R2', "orig.hospital"),
  feature_order = FALSE,
  couls = NULL,
  hcouls = c('yellow', 'black', 'blue'),
  # hcouls = rev(c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c")),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 2,
  verbose = TRUE,
  type = "pheat",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 3,
)
graphics.off()
