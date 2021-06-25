#!/usr/bin/R

############
# Figure 2 #
############

# This script will create plots for figure 2 from cd4 data

source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_plotting.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_utilities.R')
library(Seurat)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/Figure_2'
setwdc(dirfig)
clustnamebk = "none"

# Global variables
colsname <- "../data/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

### Cluster equivalence ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_clustree.R')
## Input
out_dir = dirfig
fnames = c(
  "../data/CD4T6_metadata_deprecated.rdata",
  "../data/CD4T6_metadata.rdata"
)
clustname = "../data/CD4T6_seurat.rdata"

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
void <- try(generate_tree(
  dpath = fnames,
  patfile = 'metadata_',
  cprefix = c("RNA_snn_res", "0.6"),
  libname = 'origlib',
  sobject = mycells,
  outdir = paste0(dircheck(out_dir), 'clustree'),
  v = TRUE
))
clustnamebk <- clustname

#### 2.A Principal UMAP ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
nres = "RNA_snn_res.0.6"

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
mycells@meta.data$Cluster <- factormix(mycells@meta.data[, nres])
p <- DimPlot(
  object = mycells,
  reduction = 'umap',
  group.by = "Cluster",
  label = TRUE
)
pdf('f2a_umap.pdf')
print(p)
dev.off()
pdf('f2a_umap_blank.pdf')
print(shut_it(p, lays = "ext"))
dev.off()
clustnamebk = clustname

#### 2.B Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
nres = 'RNA_snn_res.0.6'

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
dfplot <- FetchData(
  object = mycells,
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster <- factormix(dfplot[, nres])
dfplot$virus <- make_title(dfplot$orig.virus)
dfplot$virus[dfplot$orig.virus2 == "HCV"] <- "HCV"
table(dfplot$virus)

# ---- Start generating the figures
scells <- sample_grp(annot = dfplot, cname = 'virus', v = TRUE)
p <- ggplot(dfplot[scells, ], aes(x = UMAP_1, y =  UMAP_2, color = Cluster)) +
  geom_point() +
  facet_wrap(facets = ~virus) +
  labs(color = NULL) +
  mytheme +
  SetLegendPointsGG()

pdf("f2b_umap_virus.pdf", width = 15, height = 15)
print(p)
dev.off()
pdf("f2b_umap_virus_blank.pdf", width = 15, height = 15)
print(shut_it(p))
dev.off()

pp <- plot_pct(x = dfplot, groups = c("virus", "Cluster"), orderby = "CV", normalise = TRUE, return_table = TRUE)
propdf <- pp$table
write.csv(propdf, file = "f2b_proportions_virus.csv")

stupsize <- c(5, 5, 5, 5, 5)
names(stupsize) <- c("BottleRocket1", "Rushmore1", "Darjeeling1", "FantasticFox1", "Zissou1")
for(i in 3){
  couls <- colorRampPalette(wesanderson::wes_palette(name = names(stupsize)[i], n = stupsize[[i]], type = "discrete"))(stupsize[[i]])
  if(names(stupsize)[i] %in% c("FantasticFox1", "Rushmore1", "Zissou1")) couls <- rev(couls)
  p2 <- pp$plot + coord_flip() + mytheme + scale_fill_manual(values = couls)
  fname <- paste0("f2b_proportions_virus_", names(stupsize)[i])
  pdf(paste0(fname, ".pdf"), width = 5, height = 15)
  print(p2)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 5, height = 15)
  print(shut_it(p2))
  dev.off()
}


p <- plot_pct(x = dfplot, groups = c("Cluster", "virus"), normalise = FALSE, type = "donut") +
  labs(fill = NULL) +
  blank_theme + theme(
    legend.position = "right",
    strip.text = element_text(face = 'bold', size = 10),
    axis.text.x = element_blank()
  )

pdf("sf2b_donut_virus.pdf", width = 10, height = 10)
print(p)
dev.off()
pdf("sf2b_donut_virus_blank.pdf", width = 10, height = 10)
print(shut_it(p))
dev.off()

scells <- getsubset(c("orig.virus2", "CV"), dfplot, v = TRUE)
p <- plot_pct(x = dfplot[scells, ], groups = c("Cluster", "orig.project_id"), normalise = FALSE) +
  labs(fill = NULL) + coord_flip() + mytheme

pdf("sf2a_bar_project_id.pdf", width = 5, height = 10)
print(p)
dev.off()
pdf("sf2a_bar_project_id_blank.pdf", width = 4, height = 10)
print(shut_it(p))
dev.off()
clustnamebk <- clustname

#### 2.C Cluster markers ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
marknames = "../data/CD4T6_markers.csv"
nres = 'RNA_snn_res.0.6'
dgenes = c(
  'IL21', 'POU2AF1', 'CD200', 'BTLA', 'IFNG', 'IL2', 'IL3', 'TNF', 'CSF2',
  'IL17A', 'IL17F', 'CCR6', 'CTSH', 'IL4I1', 'IFIT3', 'IFI44L', 'ISG15',
  'MX2', 'OAS1', 'PRF1', 'GZMB', 'GNLY', 'NKG7', 'XCL1', 'XCL2', 'CCR7',
  'IL7R', 'TCF7', 'MKI67', 'CDK1', 'CCNB1', 'UBE2C'
)
suffix = "markers_shared"
ntop = 200
equiv = c(
  '0' = '0', '6' = '5', '7' = '7', '1' = '1', '10' = '10', '9' = '2', '2' = '8',
  '5' = '3', '4' = '6', '8' = '9', '11' = '11', '3' = '4', '12' = '12'
)
identnames = c(
  "0" = "TFH", "5" = "TFH", "7" = "TFH",
  "1" = "TH1 PolyF", "10" = "TH1 PolyF",
  "2" = "TH17", "8" = "THIFNR",
  "3" = "THIFNR",
  "6" = "CTL", "9" = "CTL", "11" = "CTL-like",
  "4" = "Central Memory",
  "12" = "Cycling"
)

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
mygenes <- readfile(marknames, stringsAsFactors = FALSE, row.names = 1)

## Operations
mycells$Name = factor(identnames[as.character(mycells@meta.data[, nres])], levels = unique(identnames))
mycells$Cluster = factor(mycells@meta.data[, nres], levels = names(identnames))
table(mycells@meta.data[, nres], mycells$Name, useNA = 'always')
table(mycells@meta.data[, nres], mycells$Cluster, useNA = 'always')
mygenes$gene <- gsub("'", "", mygenes$gene_name)

head(mygenes[dgenes, 1:11])

tvar <- mygenes$Dpct > .1; table(tvar)
mygenes <- mygenes[tvar, ]
tvar <- mygenes$avg_logFC > .25; table(tvar)
mygenes <- mygenes[tvar, ]
dgenes[!dgenes %in% mygenes$gene]
# tvar <- !grepl("&", mygenes$sCluster); table(tvar)
# mygenes <- mygenes[tvar, ]
mygenes$cluster <- factor(mygenes$cluster, names(identnames))
mygenes <- mygenes[order(mygenes$cluster), ]
# mygenes$nclust <- stringr::str_count(mygenes$sCluster, "&")
topgenes <- get_top_n(x = mygenes, n = ntop)
suffy <- paste0(suffix, ifelse(nrow(topgenes) != nrow(mygenes), paste0("_top", ntop), ""))
# topgenes <- get_top_n(x = topgenes, n = ntop, orderby = 'nclust'); suffy <- paste0(suffy, "_nclustOrder")
genes <- gsub("'", "", topgenes$gene_name)
genes <- getfound(genes, rownames(mycells), v = TRUE)
genesl <- make_list(x = topgenes, colname = "cluster", col_objects = "gene_name")

### Heatmap
fname <- paste0("f2c_heatmap_mean_", suffy)
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R')
pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
custom_heatmap(
  object = mycells,
  rnames = genes,
  orderby = "Cluster",
  use_mean = "Cluster",
  sample_it = c(cname = "Cluster", maxln = "-1000"),
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
tvar <- lapply(genesl, function(x){ y <- x[gsub("'", "", x) %in% dgenes]; print(length(y)); print(head(y)) })

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
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)
fname <- "f2d_curtain"
pdf(paste0(fname, ".pdf"), width = 10, height = 13)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 8.5, height = 13)
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()

clustnamebk <- clustname

### 2.E Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "../data/CD4T6_seurat.rdata"
allconfigs <- list(
  vconfig = list(
    name = "f2e_tfh",
    genes = c('IL21', 'POU2AF1'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(10, 5),
    order = c('0', '5', '7', 'REST')
  ),
  vconfig = list(
    name = "sf2e_tfh",
    genes = c('CD200', 'BTLA', 'CXCL13'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('0', '5', '7', 'REST')
  ),
  vconfig = list(
    name = "f2e_th1",
    genes = c('IFNG', 'IL2'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(10, 5),
    order = c('1', '10', 'REST')
  ),
  vconfig = list(
    name = "sf2e_th1",
    genes = c("TNF", "CSF2", "IL3"),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('1', '10', 'REST')
  ),
  vconfig = list(
    name = "f2e_th17",
    genes = c('IL17A', 'CCR6'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(10, 5),
    order = c('8', '2', 'REST')
  ),
  vconfig = list(
    name = "sf2e_th17",
    genes = c('IL17F', 'CTSH', 'IL4I1'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('8', '2', 'REST')
  ),
  vconfig = list(
    name = "sf2e_ifnr",
    genes = c('IFIT3', 'IFI44L', 'ISG15'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('8', '3', 'REST')
  ),
  vconfig = list(
    name = "sf2e_cd4_ctl",
    genes = c('PRF1', 'GZMB', 'GNLY'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('6', '9', 'REST')
  )
)

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
for(vconfig in allconfigs){
  cat(vconfig[['name']], "\n")
  scells <- getsubset(vconfig[['sfilter']], mycells@meta.data, v = TRUE)
  mygenes <- getfound(vconfig[['genes']], rownames(mycells), v = TRUE)
  myvars <- c(mygenes, vconfig[['ident']])
  ddf <- FetchData(mycells, vars = myvars, cells = scells)
  ddf$tmp <- do.call(paste, c(ddf[, vconfig[['ident']], drop = FALSE], sep = "_"))
  if(any(vconfig[['order']] == "REST")){
    tvar <- as.character(ddf[, vconfig[['ident']]])
    ddf$Identity <- ifelse(!tvar %in% vconfig[['order']], "REST", tvar)
  }else{ ddf$Identity <- ddf$tmp }
  if(!is.null(vconfig[['order']])) ddf$Identity <- factor(ddf$Identity, vconfig[['order']])

  # Now individually... nightmare
  dname <- paste0(vconfig[['name']])
  dir.create(dname)
  for(i in 1:length(mygenes)){
    cat(mygenes[i], "\n")
    p <- violin(
      dat = ddf,
      xax = 'Identity',
      yax = mygenes[i],
      dots = FALSE,
      colour_by = 'pct'
    )
    fname <- paste0(vconfig[['name']], "/", mygenes[i])
    pdf(paste0(fname, ".pdf"), width = 5, height = 5)
    print(p)
    dev.off()
    pdf(paste0(fname, "_blank.pdf"), width = 5, height = 5)
    print(shut_it(p))
    dev.off()
  }
}
clustnamebk <- clustname

### X.X Signature calculation ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/gsea_signature.R')
load('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/data/signatures_cd8covid.rdata')
load('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/data/signatures.rdata')
load('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/data/signatures_tem.rdata')
## Input
clustname = "../data/CD4T6_seurat.rdata"
cpatterns = c("RNA_snn_res.0.6", "orig.hospital", "RNA_snn_res.hospital")
ssigna_patterns = "tfh|cytotox|th17|interf|chmied|cycl"
out_dir = "./CD4T6_R1n2n3_sng_25p_"
pcs.comp = 38
norm_type = "RNA_snn_res"; verb = TRUE
rdims = list(umap = c('UMAP_1', 'UMAP_2'))

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
tvar <- cpatterns[cpatterns %in% colnames(mycells[[]])]
mycells$RNA_snn_res.hospital = do.call(paste, c(mycells@meta.data[, tvar, drop = FALSE], sep = "_"))
mycells$RNA_snn_res.hospital[is.na(mycells$RNA_snn_res.hospital)] <- NA # make sure asthma is correct
mycells$RNA_snn_res.hospital[mycells$RNA_snn_res.hospital == "void"] <- NA
table(mycells$RNA_snn_res.hospital)

globalsign <- vlab_signatures
tvar <- read.csv("../data/tcell_activation_signatures.csv", stringsAsFactor = FALSE)
globalsign <- c(globalsign, as.list(tvar))
globalsign <- c(globalsign, signatures_cd8covid[!signatures_cd8covid %in% names(globalsign)])
ssignatures_list <- ssignatures_listt <- globalsign[grepl(ssigna_patterns, names(globalsign))]
ssignatures_list$tcell_cytotoxicity_guo <- c("CTSW", "GNLY", "GZMA", "GZMB", "GZMH", "IFNG", "KLRB1", "KLRD1", "KLRK1", "NKG7", "PRF1")
str(ssignatures_list)
ssignatures_list <- clean_feature_list(mat = mycells@assays$RNA@data, features = ssignatures_list, filterby = "p2", v = TRUE)
final_ssignatures_list <- lapply(names(ssignatures_list), function(x){
  y <- ssignatures_list[[x]]
  z <- ssignatures_listt[[x]]
  c(y, "Not found", z[!z %in% y])
})
str(final_ssignatures_list)
ddf <- vlist2df(final_ssignatures_list)
names(ddf) <- names(ssignatures_list)
write.csv(ddf, file = paste0("supp_table_gene_signatures_", basename(dirnamen(clustname, 3)), ".csv"))

void <- try(
  signature_scoring(
    object = mycells,
    prefix = paste0(out_dir, "signatures_", pcs.comp, 'PCs/'),
    lsignatures = rev(ssignatures_list),
    confounders = cpatterns,
    reductions = rdims,
    v = verb
  )
)

clustnamebk <- clustname

### X.X GSEA metrics ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/gsea_tests.R')
library(dplyr)
library(data.table)
## Input
mymetric = "Signal2Noise"
mymethod = c("liger", "fgsea")[2]
# CD4t6
prefix = "cd4t6_gsea"
clustname = "../data/CD4T6_seurat.rdata"
nres = 'RNA_snn_res.0.6'
identnames = c(
  "0" = "TFH", "5" = "TFH", "7" = "TFH",
  "1" = "TH1 PolyF", "10" = "TH1 PolyF",
  "2" = "TH17", "3" = "THIFNR",
  "8" = "THIFNR",
  "6" = "CTL", "9" = "CTL", "11" = "CTL-like",
  "4" = "Central Memory",
  "12" = "Cycling"
)

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
prefix <- dircheck(prefix)
dir.create(prefix)
tvar <- Matrix::rowMeans(mycells@assays$RNA@data[unlist(ssignatures_list), ])
length(tvar[tvar < 0.1])
expdata <- expm1(mycells@assays$RNA@data) * 100
genes <- rownames(expdata)#[Matrix::rowMeans(expdata) > 1]
length(genes)
scells <- sample_grp(mycells@meta.data, cname = nres, maxln = -3000, v = TRUE)
expdata <- as.matrix(expdata[genes, scells])
resnamef <- paste0(prefix, "metrics_per_", nres, ".csv")
if(file.exists(resnamef)){
  res <- read.csv(resnamef, row.names = 1, stringsAsFactor = FALSE, check.names = FALSE)
  res <- setrows(res)
}else{
  res <- data.frame(
    gene_name = paste0("'", genes),
    row.names = genes
  )
}
mygroups <- if(is.factor(mycells@meta.data[, nres])) levels(mycells@meta.data[, nres]) else unique(mycells@meta.data[, nres])
void <- list()
for(myclust in mygroups){
  mynameis <- paste0(myclust, "vsREST")
  cat(mynameis, "\n")
  if(!mynameis %in% colnames(res)){
    annot <- mycells@meta.data[scells, ]
    annot$tmp <- ifelse(as.character(annot[, nres]) == myclust, yes = myclust, no = "REST")
    annot <- annot[order(annot$tmp), ]
    if(unique(annot$tmp)[1] != myclust) annot <- annot[order(annot$tmp, decreasing = TRUE), ]
    tvar <- make_list(annot, "tmp", grouping = T)
    res$s2n <- gsea_metric(
      groups = tvar,
      mat = expdata[, rownames(annot)],
      metric = 'Signal2Noise',

      v = TRUE
    )
    colnames(res) <- gsub("^s2n$", mynameis, colnames(res))
    write.csv(res, file = paste0(prefix, "metrics_per_", nres, ".csv"))
  }
  void[[myclust]] <- gsea_tests(
    res = res,
    gene_name = "gene_name",
    lfc.type = mynameis,
    gsea_file = ssignatures_list,
    method = mymethod,
    myseed = 27,
    path = paste0(prefix, mymethod, '_', mynameis, "/"),
    plot_all = TRUE,
    v = TRUE
  )
}
write.csv(res, file = paste0(prefix, "metrics_per_", nres, ".csv"))
mygseas <- rbindlist(lapply(names(void), function(x){
  y <- void[[x]]
  cbind(comparison = strsplit(basename(x), "_a1")[[1]][1], y)
}))

# Now summary table
fwrite(mygseas, file = paste0(prefix, "a1_gsea_summary_", gsub(" .*", "", Sys.time()), ".txt"), sep = "\t")

# Summary plots
# mygseas <- mygseas[!mygseas$pathway %in% c("tcell_cytotoxicity_guo", "type1n2_interferon_broadreact", "treg_signature_simone"), ]
thistab <- mygseas[mygseas$padj <= 0.05, ]
round2 <- function(x) round(x, digits = 2)
nlog10 <- function(x) -log10(x)
plot_var <- list(
  padj = list(title = "Adjuste P-value", transf = nlog10, thresh = 0.05),
  NES = list(title = "Normalized Enrichment Score", transf = round2, thresh = -Inf)
)[2] # Select type

mysum <- dcast.data.table(thistab, comparison ~ pathway, value.var = names(plot_var))
mysum <- plot_var[[1]]$transf(data.frame(mysum[, -1], stringsAsFactors = FALSE, row.names = mysum$comparison))
tmysum <- mysum <- t(as.matrix(mysum))
# mysum[is.na(mysum)] <- 0 # use it if the plot fails
annoc <- data.frame(
  Cluster = names(identnames), Identity = unname(identnames), row.names = names(identnames)
)
annoc <- annoc[, sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
anncolist <- lapply(annoc, function(x) v2cols(x, gr.cols, v = TRUE) )
cr <- cc <- TRUE
library(pheatmap)
x <- pheatmap(
  mat               = mysum,
  cluster_rows      = cr,
  cluster_cols      = cc,
  scale             = 'none',
  border_color      = NA,
  show_colnames     = T,
  show_rownames     = T,
  main              = plot_var[[1]]$title,
  annotation_col    = annoc,
  annotation_colors = anncolist,
  na_col            = "#BEBEBE",
  annotation_legend = T,
  annotation_names_col = T,
  annotation_names_row = T,
  drop_levels       = TRUE,
  filename          = paste0(prefix, "a1_gsea_summary_", names(plot_var), "_", gsub(" .*", "", Sys.time()), ".pdf"),
  width = 10, height = 7
);

# in case it failed when not setting NAs to 0
mysum <- tmysum[x$tree_row$labels[x$tree_row$order], x$tree_col$labels[x$tree_col$order]]
cr <- cc <- FALSE


#### 2.F Signatures ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
nres = 'RNA_snn_res.0.6'
signaturef = "../Figure_2/CD4T6_R1n2n3_sng_25p_signatures_38PCs/signatures.csv"
ssigna_patterns = "tfh|cytotox|th17|interf|cyclying|chmied|cycl"
prefix = c("f2f_", "sf2e_")
couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
signaturedf <- readfile(signaturef, row.names = 1, stringsAsFactors = FALSE)

## Operations
tvar <- grepl(ssigna_patterns, colnames(signaturedf), ignore.case = TRUE)
ssignatures <- colnames(signaturedf)[tvar & grepl("Score", colnames(signaturedf))]
ddf <- signaturedf[, ssignatures]
colnames(ddf) <- ssignatures <- gsub(".Score", "", ssignatures)
ddf <- joindf(ddf, FetchData(mycells, vars = c("UMAP_1", "UMAP_2", nres)))
ddf$Cluster <- factormix(ddf[, nres])
ddf <- ddf[order(ddf$Cluster), ]
str(ddf)
sapply(ddf, function(x) sum(is.na(x)) )

for(ssig in ssignatures){
  cat(ssig, "\n")
  summary(ddf[, ssig]); quantile(ddf[, ssig], seq(.85, .99, 0.02))
  topit <- head(sort(ddf[, ssig], T), 3)
  topit <- topit[length(topit)]
  ddf$Signature <- ddf[, ssig];
  ddf$Signature[ddf$Signature > topit] <- topit
  p <- ggplot(ddf, aes(x = UMAP_1, y = UMAP_2, color = Signature)) +
    geom_point(size = 0.3) +
    scale_colour_gradientn(
      colours = couls
    ) +
    theme(legend.position = "right") + labs(color = NULL)
  fname <- paste0(prefix[1], ssig)
  pdf(paste0(fname, "_umap.pdf"), width = 8, height = 8)
  print(p)
  dev.off()
  pdf(paste0(fname, "_umap_blank.pdf"), width = 8, height = 8)
  print(shut_it(p))
  dev.off()

  means <- tapply(X = ddf$Signature, INDEX = ddf$Cluster, mean)
  ddf$mSignature <- means[as.character(ddf$Cluster)]
  p <- violin(
    dat = ddf,
    xax = "Cluster",
    yax = "Signature",
    dots = FALSE,
    colour_by = "mSignature"
  ) + theme(legend.title = element_blank())
  pdf(paste0(prefix[2], ssig, "_violin_mean.pdf"), width = 12, height = 8)
  print(p)
  dev.off()
}
clustnamebk <- clustname

#### S2.A QC Violin plots ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = c(
  "data/CD4T0n6_seurat.rdata",
  "../data/CD4T24_seurat.rdata"
)
low.filt.cells = c(800, 1500, -Inf)
high.filt.cells = c(4400, 20000, 10)
subs.names = c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
names(low.filt.cells) <- subs.names
names(high.filt.cells) <- subs.names

## Reading
lmycells <- lapply(clustname, theObjectSavedIn)
mycellscom <- if(length(lmycells) > 1){
  tvar <- lmycells[[1]]
  for(i in 2:length(lmycells)){
    tvar <- merge(tvar, lmycells[[i]])
  }
  tvar
}else{ lmycells }
rm(lmycells)
tvar <- reshape2::melt(table(mycellscom$orig.stim_time, mycellscom$origlib))
tvar <- tvar[tvar[, 3] > 0, ]
tvar <- tvar[order(tvar[, 1]), ]
mycellscom$Library <- factor(mycellscom$origlib, levels = unique(as.character(tvar[, 2])))
table(mycellscom$Library)

## Operations
mycols <- rep("#2e82ff", length(table(mycellscom$Library)))
names(mycols) <- names(table(mycellscom$Library))
for(qcvar in subs.names[1]){
  p <- VlnPlot(
    object = mycellscom,
    features = qcvar,
    group.by = "Library",
    cols = mycols
  ) + NoLegend()
  fname <- paste0("sf2a_qc_", qcvar)
  pdf(paste0(fname, ".pdf"), width = 12)
  print(rm_layer(p, "point"))
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 12)
  print(shut_it(p, "point|txt"))
  dev.off()
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%#### Check TFR cells ####%%%%%%%%%%%%%%%%%%%%%%%%%%%#
setwdc(paste0(dirfig, "../tfr_exploration"))
## Input
clustname = "../data/CD4T24_seurat.rdata"
genes = c("FOXP3", "IL1R2", "TNFRSF9", "TGFB1", "CCR8", "TNFRSF18", "TOP2A", "PDCD1", "TIGIT", "MKI67", "CTLA4", "IL10", "IL2RA", "BATF", "BCL6", "CXCR5")
nres = "RNA_snn_res.0.2"
sselect = list(c("RNA_snn_res.0.2", "0"))

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
mygenes <- getfound(genes, rownames(mycells), v = TRUE)
for(g in mygenes){
  cat(g, "\n")
  p <- FeaturePlot(
    object = mycells,
    features = g,
    cols = c("#fffeee", "#ffe080", "#ffc100", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000")
  )
  pdf(paste0(g, ".pdf"))
  print(p)
  dev.off()
}
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
ddf <- FetchData(
  object = mycells,
  cells = scells,
  vars = c(colnames(mycells[[]]), genes)
)
tvar <- table(ddf$RNA_snn_res.0.3)/table(mycells@meta.data$RNA_snn_res.0.3) > .1
tvar <- names(table(ddf$RNA_snn_res.0.3)[tvar])
tvar
ddf <- remove.factors(ddf[as.character(ddf$RNA_snn_res.0.3) %in% tvar, ])
void <- add_gene_tag(
  lgenes = "IL1R2",
  annot = ddf,
  mat = mycells@assays$RNA@data,
  v = TRUE
)
ddf <- joindf(ddf, void)
for(g in mygenes){
  cat(g, "\n")
  p <- violin(
    dat = ddf,
    xax = c("tag_IL1R2", "RNA_snn_res.0.3")[2],
    yax = g,
    colour_by = 'pct',
    chilli = FALSE
  )
  pdf(paste0("violin0.3_", g, ".pdf"))
  print(p)
  dev.off()
}
clustnamebk = clustname

### Tables per donor ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table)
## Input
clustname = "../data/CD4T24_seurat.rdata"
genes = c('IL1R2', 'CCR8')
sselect = list(c('RNA_snn_res.0.2', '0'), c('orig.donor', '-void'))

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
mygenes <- getfound(genes, rownames(mycells), v = TRUE)
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
void <- add_gene_tag(
  lgenes = mygenes,
  annot = mycells@meta.data[scells, ],
  mat = mycells@assays$RNA@data,
  v = TRUE
)
void$combn <- apply(void[, , drop = FALSE], 1, paste, collapse = "")
annot <- joindf(void, mycells@meta.data[scells, ])
table(annot[, 'orig.donor'])
ddf <- rbindlist(lapply(colnames(void), function(x){
  y <- table(annot[, c(x, "orig.donor") ])
  cbind(rownames(y), as.data.frame.matrix(y))
}))
ddfsumm <- data.frame(t(ddf[, -1]))
colnames(ddfsumm) <- as.character(ddf[[1]])
tvar <- paste0('Total_in_cluster', paste0(sselect[[1]][-1], collapse = "n"))
ddfsumm[, tvar] <- rowSums(ddfsumm[, 1:2])
ddfsumm$TotalDonor <- table(mycells@meta.data[, "orig.donor"])[rownames(ddfsumm)]
head(ddfsumm)

annot$dcluster <- do.call(paste, c(annot[, c("orig.donor"), drop = FALSE], sep = "-"))
mygroups <- make_list(
  x = annot,
  colname = "dcluster", grouping = TRUE
)
mystats <- get_stat_report(
  mat = mycells@assays$RNA@data,
  groups = mygroups,
  moments = c("mn", "p", "pmn", "md", "pmd"),
  rnames = mygenes,
  v = TRUE
)
tvar <- data.frame(t(mystats))
head(tvar)
if(any(grep("\\-", rownames(tvar)))){
  tvar$Donor <- gsub("\\-.*", "", rownames(tvar))
  tvar$Cluster <- gsub(".*\\-", "", rownames(tvar))
}else{
  tvar$Donor <- gsub("_[A-z].*", "", rownames(tvar))
  tvar$Cluster <- gsub(".*_([A-z].*)", "\\1", rownames(tvar))
}
head(tvar)
library(data.table)
mysum <- data.frame(dcast.data.table(data.table(tvar), Donor ~ Cluster, value.var = mygenes))
rownames(mysum) <- mysum[, 1]
mysum <- mysum[, -1]
headmat(mysum)
mysum <- joindf(ddfsumm, mysum)
fname <- paste0("summary_per_", gsub("orig.", "", "orig.donor"), "_", paste0(mygenes, collapse = "n"), ".csv")
write.csv(mysum, file = fname)

clustnamebk <- clustname
