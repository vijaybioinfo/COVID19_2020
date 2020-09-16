#!/usr/bin/R

############
# Figure 4 #
############

# This script will create plots for figure 4 from cd4 data

source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_plotting.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_utilities.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/Figure_4'
setwdc(dirfig)
clustnamebk = "none"

# Global variables
colsname <- "../data/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

### 4.C Trajectory ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(monocle3)
## Input
clustname = "../data/CD4T6_monocle.rdata"

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
) + theme_cowplot() + theme(legend.title = element_blank())
pdf(paste0("f4c_trajectory_umap.pdf"), width = 8, height = 8)
print(p)
dev.off()
p2 <- shut_it(p, "text")
p2$layers <- p2$layers[1:2]
pdf(paste0("f4c_trajectory_umap_blank.pdf"), width = 8, height = 8)
print(p2)
dev.off()

#### 4.D Principal UMAP ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T24_seurat.rdata"
nres = 'RNA_snn_res.0.2'

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
str(mycells@meta.data[, nres])
mycells@meta.data[as.character(mycells@meta.data[, nres]) == "6", nres] <- "0"
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

clustnamebk <- clustname

#### 2.E Cluster markers ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T24_seurat.rdata"
marknames = "../data/CD4T24_markers.csv"
nres = 'RNA_snn_res.0.2'
suffix = "markers_shared"
ntop = 200
identnames = c(
  "0" = "Treg", "6" = "Treg", "1" = "CTL", "5" = "CTL", "2" = "TCM", "3" = "TFH", "4" = "TH17"
)

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
mygenes <- readfile(marknames, stringsAsFactors = FALSE)

## Operations
mycells@meta.data[as.character(mycells@meta.data[, nres]) == "6", nres] <- "0"
identnames <- identnames[names(identnames) %in% as.character(mycells@meta.data[, nres])]
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
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R')
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

### Dot-plots
p <- DotPlot(
  object = mycells,
  features = dgenes,
  group.by = "Cluster",
  cols = c('#fff4ba', '#ff0000'),
  col.min = -1.5, col.max = 1.5
)  + coord_flip() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold.italic"),
    # axis.text.x = element_text(face = "bold", hjust = 1, angle = 45),
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)
fname <- "f4f_curtain"
pdf(paste0(fname, ".pdf"), width = 10, height = 13)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 10, height = 13)
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()

clustnamebk <- clustname

### 4.X Expression plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# CD4T24
myset = "nf4_umapCD4T24"
clustname = "../data/CD4T24_seurat.rdata"
genes <- c("FOXP3", "IL1R2")

#### 2.X Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T24_seurat.rdata"
nres = 'RNA_snn_res.0.2'
norm_type = "RNA_snn_res."
couls = c(Hospital = "#fb8071", "non-Hospital" = "#b3de69")
identnames <- LETTERS[1:7]
names(identnames) <- 0:6
sselect = list(c('orig.hospital', '-void'))

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
mycells@meta.data[as.character(mycells@meta.data[, nres]) == "6", nres] <- "0"
dfplot <- FetchData(
  object = mycells,
  cells = getsubset(sselect, mycells@meta.data, v = TRUE),
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster = factormix(identnames[as.character(dfplot[, nres])])
dfplot$Cluster <- factormix(dfplot[, nres])
table(dfplot$orig.hospital)
dfplot$hospital <- ifelse(dfplot$orig.hospital == "Yes", "Hospital", "non-Hospital")
dfplot$hospital <- factor(dfplot$hospital, c("Hospital", "non-Hospital"))
table(dfplot$hospital, useNA = 'always')

pp <- plot_pct(x = dfplot, groups = c("hospital", "Cluster"), orderby = "Hospital", normalise = TRUE, return_table = TRUE)
propdf <- pp$table
write.csv(propdf, file = "f4h_proportions.csv")
p <- pp$plot + coord_flip() + scale_fill_manual(values = couls) + labs(fill = NULL)

pdf("f4h_hospital.pdf", width = 5, height = 10)
print(p)
dev.off()
pdf("f4h_hospital_blank.pdf", width = 5, height = 10)
print(shut_it(p + mytheme))
dev.off()

p <- plot_pct(x = dfplot, groups = c("Cluster", "hospital"), type = "donut", print_ptables = TRUE) +
  labs(fill = NULL) +
  blank_theme + theme(
    legend.position = "right",
    strip.text = element_text(face = 'bold', size = 10),
    axis.text.x = element_blank()
  )

pdf("f4i_donut_hospital.pdf", width = 5, height = 10)
print(p)
dev.off()
pdf("f4i_donut_hospital_blank.pdf", width = 5, height = 10)
print(shut_it(p + mytheme))
dev.off()

### X.X Signature calculation ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input - signature
clustname = "../data/CD4T24_seurat.rdata"
cpatterns = c("RNA_snn_res.0.2", "orig.hospital", "RNA_snn_res.hospital")
ssigna_patterns = "cytotox|tfh|th17|treg.*chmied|interf|cycl|tfr"
out_dir = "./CD4T24_R1n2n3_sng_25p_"
pcs.comp = 16
norm_type = "RNA_snn_res"; verb = TRUE
rdims = list(umap = c('UMAP_1', 'UMAP_2'))
nres = "RNA_snn_res.0.2"

## Input - gsea
identnames = c(
  "0" = "Treg", "1" = "CTL", "5" = "CTL", "2" = "TCM", "3" = "TFH", "4" = "TH17"
)
prefix = "cd4t24_gsea"

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
mycells@meta.data[as.character(mycells@meta.data[, nres]) == "6", nres] <- "0"
mycells@meta.data[, nres] <- factormix(mycells@meta.data[, nres])
table(mycells@meta.data[, nres])
identnames <- identnames[names(identnames) %in% as.character(mycells@meta.data[, nres])]

# run fig2.R - X.X Signature calculation; and X.X GSEA metrics

### X.X Signature calculation on 24h ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T24_seurat.rdata"
cpatterns = c("RNA_snn_res.0.2", "orig.hospital", "RNA_snn_res.hospital")
ssigna_patterns = "tfr|tfh|patil"
out_dir = "./CD4T24_R1n2n3_sng_25p_"
pcs.comp = 16
norm_type = "RNA_snn_res"; verb = TRUE
rdims = list(umap = c('UMAP_1', 'UMAP_2'))

# Run fig2.R - X.X Signature calculation

#### X.X 24h Signatures ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T24_seurat.rdata"
nres = "RNA_snn_res.0.2"
signaturef = "./CD4T24_R1n2n3_sng_25p_signatures_16PCs/signatures.csv"
ssigna_patterns = "tfr|tfh|patil"
prefix = c("nf4_24_", "nf4_24_")
couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")

# Run fig2.R - 2.F Signatures

#### 4.X Cluster proportions per donor ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
# CD4T24
clustname = "../data/CD4T24_seurat.rdata"
nres = 'RNA_snn_res.0.2'
dcolumn = 'orig.donor'
subdiv = c("orig.hospital", "orig.sex")
sselect <- list(c("orig.hospital", "Yes", "No"), c("orig.virus2", "CV"))
corder = c('0', '1', '5', '3', '2', '4')
suffy = NULL
roder = NULL

add_props = list(
  newobject = NULL,
  filter = c("RNA_snn_res.0.2", "0"),
  vars = c("orig.donor", "tag_IL1R2"),
  totals = "orig.donor"
)
corder = c('0', '1', '5', '3', '2', '4', 'IL1R2+', 'IL1R2-')
suffy = "_il1r2"

corder = c('0', '1', '3', '2', 'IL1R2+', 'IL1R2-')
suffy = "_abcd_il1r2"

corder = c('0', '1', '3', '2')
suffy = "_abcd"
add_props = NULL

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
str(mycells@meta.data[, nres])
mycells@meta.data[as.character(mycells@meta.data[, nres]) == "6", nres] <- "0"
mycells@meta.data$Cluster <- factormix(mycells@meta.data[, nres])

# Run fig2.R - 2.F Cluster proportions per donor

### X.X 0n6 Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "data/CD4T0n6_seurat.rdata"
allconfigs <- list(
  vconfig = list(
    name = "nf4_0n6_tfh",
    genes = c('BTLA', 'CD200', 'POU2AF1'),
    sfilter = list(c('orig.stim_time', '0', '6')),
    ident = c('orig.stim_time'),
    ncols = 2, size = c(12, 8),
    order = c('0', '6')
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%#### Check TFH cells ####%%%%%%%%%%%%%%%%%%%%%%%%%%%#
setwdc(paste0(dirfig, "../tfh_exploration"))
### X.X Signature/GSEA vs Trajectory ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TFH, cytotoxicity, interferon, exhaustion signatures on cluster 5 vs 0

## Input - Signature
clustname = "../data/CD4T6_seurat.rdata"
cpatterns = c("RNA_snn_res.0.6", "orig.hospital", "RNA_snn_res.hospital")
ssigna_patterns = "tfh|cytotox|interf|xhaus"
out_dir = "CD4T6_R1n2n3_sng_25p_"
pcs.comp = 38
norm_type = "RNA_snn_res"; verb = TRUE
rdims = list(umap = c('UMAP_1', 'UMAP_2'))

# run fig2.R - X.X Signature calculation

## Input - GSEA
mynameis = c("cv_5vs0", "cv_5vs0n7")[2]
resnamef = paste0("../data/", mynameis, ".csv")
mymethod = c("liger", "fgsea")[2]
ssigna_patterns = "tfh|cytotox|interf|xhaus"

# run fig2.R - X.X GSEA

#### X.X Signature/gene vs or over Dimentions/pseudotime ####%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
mononame = "monocle_plot_selected_root.rdata"
signaturef = "../Figure_2/CD4T6_R1n2n3_sng_25p_signatures_38PCs/signatures.csv"
genes = c("IFI6", "ISG15", "ISG20", "IFITM1", "IFI35", "IFI44L")
nres = "RNA_snn_res.0.6"
sselect = list(c('orig.virus2', 'CV'), c("RNA_snn_res.0.6", '0', '5'))

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
monodims = theObjectSavedIn(mononame)
signaturedf <- readfile(signaturef, row.names = 1, stringsAsFactors = FALSE)

## Operations
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
str(monodims)
ddf <- remove.factors(joindf(monodims, mycells@meta.data[rownames(monodims), ]))
ddf <- joindf(ddf, signaturedf)
table(ddf[, nres])
table(ddf[scells, nres])
table(mycells@meta.data[scells, nres], ddf[scells, nres])
str(ddf)

# .libPaths('~/R/newer_packs_library/3.5') ## Using another R version !!!!!
# source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
# library(monocle)
# mononame = "../data/CD4T6_tfh_monocomponents.rdata"
# monobjname = "../data/CD4T6_tfh_monocle.rdata"
# cds <- theObjectSavedIn(monobjname)
# str(cds@reducedDimK) # nodes for pseudotime
# MPP_node_ids <- get_correct_root_state(cds, cell_phenotype = nres)
# cat('Ordering cells, root:', MPP_node_ids, '\n')
# comps <- data.frame(t(cds@reducedDimK)); colnames(comps) <- c('MC1', 'MC2')
# comps$Position <- rownames(comps)
#
# p <- ggplot(ddf, aes_string(x = 'MC1', y = 'MC2')) +
#   geom_point(alpha = 0.5) +
#   geom_text(data = comps, aes_string(x = 'MC1', y = 'MC2', label = "Position"), size = 1) +
#   labs(color = NULL) +
#   SetLegendPointsGG()
#
# pdf(paste0("monocle_plot_select_root.pdf"))
# print(p)
# dev.off()
# MPP_node_ids = "Y_137"
# cds <- orderCells(cds, root_pr_nodes = MPP_node_ids)
# pdf(paste0("monocle_plot_selected_root.pdf"))
# print(plot_cell_trajectory(cds))
# graphics.off()
# fname <- paste0('monocle_plot_selected_root.rdata')
# comps <- data.frame(t(cds@reducedDimS)); colnames(comps) <- c('MC1', 'MC2')
# comps$Pseudotime <- cds@phenoData@data$Pseudotime
# head(comps)
# save(comps, file = fname)

# Signatures on dimentions
# Scatter with signatures and genes; fitting line
# 'stats::loess()' is used for less than 1,000 observations;
# otherwise 'mgcv::gam()' is used with 'formula = y ~ s(x, bs =
# "cs")' with 'method = "REML"'. Somewhat anecdotally, 'loess'
# gives a better appearance, but is O(N^2) in memory, so does
# not work for larger datasets.
cluster_cols <- gg_color_hue(n = gtools::mixedsort(levels(mycells@meta.data[, nres])))

# Monocle trajectory map
mycols <- cluster_cols[unique(ddf[, nres])]
p <- ggplot(ddf, aes_string(x = 'MC1', y = 'MC2', color = nres)) +
  geom_point() +
  scale_color_manual(values = mycols) +
  labs(color = NULL) +
  SetLegendPointsGG()

pdf(paste0("monocle_plot.pdf"))
print(p)
dev.off()
pdf(paste0("monocle_plot_blank.pdf"))
print(shut_it(p))
dev.off()

void <- FetchData(
  object = mycells,
  vars = c("PRDM1", "PRF1", "ZBED2")
)
ddf <- joindf(ddf, void)
ddf$PRDM1 <- mycells@assays$RNA@data["PRDM1", rownames(ddf)]
summary(ddf$PRDM1)

# Signature on trajectory map
thisignas <- c("PRDM1", "PRF1", "ZBED2", "Pseudotime")
mycols <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
thisignas <- "TYPE1N2_INTERFERON_BROADREACT.Score"
mycols <- c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")
for(thisigna in thisignas){
  p <- ggplot(ddf, aes_string(x = "MC1", y = "MC2", color = thisigna)) +
    geom_point(size = 0.3) +
    scale_colour_gradientn(colours = mycols) +
    theme(legend.position = "right") + labs(color = NULL)

  fname <- paste0("monocle_plot_", gsub(".Score", "", thisigna))
  pdf(paste0(fname, ".pdf"))
  print(p)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"))
  print(shut_it(p))
  dev.off()
}

mycols <- cluster_cols[unique(ddf[scells, nres])]
for(colby in nres){
  for(xname in c("MC1", "Pseudotime", "MC2")){
    for(yname in c("CELL_CYTOTOXICITY_PATIL.Score", "TYPE1N2_INTERFERON_BROADREACT.Score", "PRDM1", "PRF1", "ZBED2")){
      p <- ggplot(ddf[scells, ], aes_string(x = xname, y = yname, color = colby)) +
        geom_point(alpha = 0.3, size = 0.2) +
        geom_smooth(mapping = aes_string(color = colby), se = FALSE) +
        scale_color_manual(values = mycols) +
        theme(legend.title = element_blank())
      pdf(paste0("trend_line_", colby, "_", xname, "_vs_", yname, ".pdf"))
      print(p)
      dev.off()
    }
  }
}

# Classifying Cytotoxic cells in all TFH
void <- add_gene_tag(
  lgenes = c("GZMB", "PRF1"),
  annot = ddf,
  mat = mycells@assays$RNA@data,
  v = TRUE
)
void$cytotoxic <- do.call(paste, c(void[, 1:2, drop = FALSE], sep = ""))
ddf <- joindf(ddf, void)
ddf$cytotoxic <- void$cytotoxic
table(ddf$cytotoxic, void[rownames(ddf), ]$cytotoxic, useNA = 'always')
ddf$cytotoxic <- ifelse(ddf$cytotoxic == "GZMB-PRF1-", "GZMB-PRF1-", "GZMB+orPRF1+")

for(colby in c("cytotoxic")){
  for(xname in c("MC1", "Pseudotime", "MC2")){
    for(yname in c("CELL_CYTOTOXICITY_PATIL.Score", "TYPE1N2_INTERFERON_BROADREACT.Score", "PRDM1", "PRF1", "ZBED2")){
      p <- ggplot(ddf[scells, ], aes_string(x = xname, y = yname, color = colby)) +
        geom_point(alpha = 0.3, size = 0.2) +
        geom_smooth(mapping = aes_string(color = colby), se = FALSE) +
        theme(legend.title = element_blank())
      pdf(paste0("trend_line_", colby, "_", xname, "_vs_", yname, ".pdf"))
      print(p)
      dev.off()
    }
  }
}
