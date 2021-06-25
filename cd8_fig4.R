#!/usr/bin/R

############
# Figure 4 #
############

# This script will create plots for figure 4 from cd8 data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd8/Figure_4'
setwdc(dirfig)
clustnamebk = ""

# Global variables
colsname <- "/home/ciro/covid19/info/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

### 3.B Heatmap cluster 0 and 2 hospital comparisons ###%%%%%%%%%%%%%%%%%%%%%%%
## Input
fnames = c(
  cluster_0 = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster0severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv",
  cluster_2 = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster2severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv"
)
selectss = list(c("RNA_snn_res.0.2", "0", "2"), c('orig.hospital', 'No', 'Yes'))
threshes <- c(pv = 0.05, fc = 0.25)
idents = c("RNA_snn_res.0.2", "orig.hospital")

## Reading
resl <- lapply(fnames, readfile, row.names = 1, stringsAsFactors = FALSE)

## Operations
genesl <- lapply(resl, getDEGenes, pv = threshes['pv'], fc = threshes['fc'], upreg = TRUE, v = TRUE)
genesl <- overlap_calc(
  groups_list = genesl, v = TRUE
)
str(genesl)
genes <- genesl[['cluster0~cluster2']]
mycells$Group <- do.call(paste, c(mycells@meta.data[, idents, drop = FALSE], sep = "_"))
scells <- getsubset(selectss, mycells@meta.data, v = TRUE)
table(mycells@meta.data[scells, c("orig.donor", "Group")])
fname <- paste0("f4a_heatmap")
source('/home/ciro/scripts/functions/pheatmapCorrection.R')
pdf(paste0(fname, ".pdf"), width = 10, height = 12, onefile = FALSE)
custom_heatmap(
  object = mycells[, scells],
  rnames = genes,
  orderby = c("Group", "pca"),
  use_mean = "orig.donor",
  # sample_it = c(cname = "Group", maxln = "-3000"),
  scale_row = TRUE,
  categorical_col = c("orig.severity", "orig.hospital", "orig.donor"),
  feature_order = "hclust",
  couls = NULL,
  hcouls = rev(c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c")),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = 1,
  verbose = TRUE,
  type = "pheat",
  show_rownames = FALSE,
  show_colnames = FALSE
)
graphics.off()

### 1.X crater plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
fnames = c(
  cluster0 = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster0severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv",
  cluster2 = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster2severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv"
)
selectss = list(c("RNA_snn_res.0.2", "0", "2"))
mylfcthresh = "0.25"
thesegenes = c(
  "NFKB1", "NFKB2", "REL", "RELB", "IL2RA", "STAT5A", "BHLHE40", "CRTAM", "PRDM1",
  "JUN", "BIRC3", "VIM", "BCL2L1", "VIM", "BAX", "BIM", "FAS"
)
degfilt = list(mean = list("<=0", NA), min_padj = list(">0.05", 1), cluster2.padj = list(">0.05", 1), cluster0.padj = list(">0.05", 1))[1:2]
outdir = "crater_plots/"

## Reading and Operations
onlythesegenes <- rownames(mycells)
onlythesegenes <- onlythesegenes[!grepl("^rps|^rpl|^mt-", onlythesegenes, ignore.case = TRUE)]
dir.create(outdir)
source('/mnt/BioHome/ciro/scripts/functions/crater_plot.R')
void <- crater_plot(
  fname1 = fnames,
  lfc = "log2FoldChange",
  pv = "padj",
  keep_nas = FALSE,
  topgenes = thesegenes,
  fedata = mycells@assays$RNA@data,
  fannot = mycells@meta.data,
  sample_filter = selectss,
  feature_filter = onlythesegenes,
  gene_filter = degfilt,
  lfcthresh = mylfcthresh,
  outputname = outdir,
  return_out = FALSE,
  v = TRUE
)
setwd(dirfig)

### 4.X Curtain plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
vconfig = list(
  name = "f4d_curtain_C0n2",
  genes = c('NFKB1', 'RELA', 'IL2RA', 'BHLHE40', 'PRDM1', 'CRTAM', 'SREBF2', 'BIRC3', 'BCL2L1', 'BCL2L11', 'FAS'),
  sfilter = list(c('RNA_snn_res.0.2', '0', '2'), c('orig.virus2', 'CV')),
  ident = c('RNA_snn_res.0.2', 'orig.hospital'),
  order = c('0_No', '0_Yes', '2_No', '2_Yes')
)
vconfig = list(
  name = "f4d_curtain_C0n2_noSRE",
  genes = c('NFKB1', 'RELA', 'IL2RA', 'BHLHE40', 'PRDM1', 'CRTAM', 'BIRC3', 'BCL2L1', 'BCL2L11', 'FAS'),
  sfilter = list(c('RNA_snn_res.0.2', '0', '2'), c('orig.virus2', 'CV')),
  ident = c('RNA_snn_res.0.2', 'orig.hospital'),
  order = c('0_No', '0_Yes', '2_No', '2_Yes')
)
vconfigs <- list(
  vconfig = list(
    name = "f4_curtain_nfkb_C",
    genes = c('NFKB1', 'NFKB2', 'REL', 'RELB'),
    sfilter = list(c('RNA_snn_res.0.2', NULL), c('orig.virus2', 'CV'), c('orig.hospital', 'Yes', 'No')),
    ident = c('RNA_snn_res.0.2', 'orig.hospital'),
    order = c('0_No', '0_Yes', '2_No', '2_Yes')
  ),
  vconfig = list(
    name = "f4_curtain_prmd1_bhl_jun_C",
    genes = c('PRDM1', 'BHLHE40', 'JUN'),
    sfilter = list(c('RNA_snn_res.0.2', NULL), c('orig.virus2', 'CV'), c('orig.hospital', 'Yes', 'No')),
    ident = c('RNA_snn_res.0.2', 'orig.hospital'),
    order = c('0_No', '0_Yes', '2_No', '2_Yes')
  ),
  vconfig = list(
    name = "f4_curtain_bcls_C",
    genes = c('BIRC3', 'BCL2L1', 'BCL2A1', 'BCL2L11', 'VIM', 'MCL1', 'FAS'),
    sfilter = list(c('RNA_snn_res.0.2', NULL), c('orig.virus2', 'CV'), c('orig.hospital', 'Yes', 'No')),
    ident = c('RNA_snn_res.0.2', 'orig.hospital'),
    order = c('0_No', '0_Yes', '2_No', '2_Yes')
  )
  vconfig = list(
    name = "f4_curtain_C",
    genes = c('CRTAM'),
    sfilter = list(c('RNA_snn_res.0.2', NULL), c('orig.virus2', 'CV'), c('orig.hospital', 'Yes', 'No')),
    ident = c('RNA_snn_res.0.2', 'orig.hospital'),
    order = c('0_No', '0_Yes', '2_No', '2_Yes')
  )
  vconfig = list(
    name = "f4_curtain_il2_C",
    genes = c('IL2RA', 'STAT5A'),
    sfilter = list(c('RNA_snn_res.0.2', NULL), c('orig.virus2', 'CV'), c('orig.hospital', 'Yes', 'No')),
    ident = c('RNA_snn_res.0.2', 'orig.hospital'),
    order = c('0_No', '0_Yes', '2_No', '2_Yes')
  )
)

## Reading
mycells <- theObjectSavedIn(clustname)

## Operations
for(vconfig in vconfigs){
  for(cl in c("0", "2")){
    vconfig[['sfilter']][[1]][2] <- cl
    scells <- getsubset(vconfig[['sfilter']], mycells[[]], v = TRUE)
    mycells$tmp <- do.call(paste, c(mycells[[]][, vconfig[['ident']], drop = FALSE], sep = "_"))
    if(any(vconfig[['order']] == "REST")){
      tvar <- as.character(mycells$tmp)
      mycells$Identity <- ifelse(!tvar %in% vconfig[['order']], "REST", tvar)
    }else{ mycells$Identity <- mycells$tmp }
    mycells$Identity <- factor(mycells$Identity, vconfig[['order']])
    table(mycells$tmp, mycells$Identity)
    scells <- scells[!is.na(mycells@meta.data[scells, ]$Identity)]
    unique(mycells@meta.data[scells, ]$Identity)

    # Dot-plots
    if(length(vconfig[['genes']]) > 1){
      p <- DotPlot(
        object = mycells[, scells],
        features = getfound(vconfig[['genes']], rownames(mycells), v = TRUE),
        group.by = 'Identity',
        cols = c('#fff4ba', '#ff0000'),
        col.min = -1.5, col.max = 1.5
      ) + coord_flip() +
        theme(
          axis.text.y = element_text(size = 13, face = "bold.italic"),
          axis.ticks.x = element_blank()
        ) + labs(y = NULL, x = NULL)
      fname <- paste0(vconfig[['name']], cl)
      pdf(paste0(fname, ".pdf"), width = 5, height = 5)
      print(p)
      graphics.off()
      pdf(paste0(fname, "_blank.pdf"), width = 5, height = 5)
      print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
      graphics.off()
    }

    # Violins
    myvars <- c(vconfig[['genes']], vconfig[['ident']], "Identity")
    ddf <- FetchData(mycells, vars = myvars, cells = scells)
    dname <- "f4_violins"
    dir.create(dname, showWarnings = FALSE)
    for(i in 1:length(vconfig[['genes']])){
      cat(vconfig[['genes']][i], "\n")
      p <- violin(
        dat = ddf,
        xax = 'Identity',
        yax = vconfig[['genes']][i],
        dots = FALSE,
        colour_by = 'pct'
      )
      fname <- paste0(dname, "/", gsub("_curtain", "", vconfig[['name']]), cl, vconfig[['genes']][i])
      pdf(paste0(fname, ".pdf"), width = 5, height = 5)
      print(p)
      dev.off()
      pdf(paste0(fname, "_blank.pdf"), width = 5, height = 5)
      print(shut_it(p))
      dev.off()
    }
  }
}
clustnamebk = clustname

### 4.X Scatter plots cluster 0 and 2 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
# in cluster 2
prefix = "f4_C2"
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
mycname = "orig.hospital"
sselect = list(c("RNA_snn_res.0.2", "2"), c("orig.hospital", "No", "Yes"))
genesdf = data.frame(
  g1 = c('IL2RA', 'FAS'),
  g2 = c('STAT5A', 'BCL2L1'),
  stringsAsFactors = FALSE
)
# in cluster 0
prefix = "f4_C0"
sselect = list(c("RNA_snn_res.0.2", "0"), c("orig.hospital", "No", "Yes"))
# run scatter section in fig3

#### S4.A,S4.B Volcano ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input A
prefix = "f3a_"
padjthr = 0.05
fcthr = .5
showgenes = NULL
resname = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster0severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv"
trimmer = "c300"

## Input B
prefix = "f3a_"
padjthr = 0.05
fcthr = .25
showgenes = NULL
resname = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster2severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv"
trimmer = "c300"

## Run Volcano part in fig3

### 4.C GSEA plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(pheatmap)
library(data.table)
source('/mnt/BioHome/ciro/scripts/gsea/gsea_liger.R')
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
signaturef = "../Figure_2/supptable_signatures.csv"
ssigna_patterns = "fuch|interf|ty_guo|_best"
ssigna_patterns = "sensus|interf|ty_guo|s_cullen|fuch|_best"
sselect = list(c("RNA_snn_res.0.2", "0", "2"), c("orig.hospital", "No", "Yes"))
idents = c("RNA_snn_res.0.2", "orig.hospital")
mygroups = list(c("0_Yes", "0_No"), c("2_Yes", "2_No"))
mymethod = "fgsea"

## Reading
if(clustname != clustnamebk || !exists("mycells")) mycells <- theObjectSavedIn(clustname)
signaturedf <- readfile(signaturef, row.names = 1, stringsAsFactors = FALSE)

## Operations
tvar <- grepl(ssigna_patterns, colnames(signaturedf), ignore.case = TRUE)
ssignatures <- colnames(signaturedf)[tvar]
ssignatures
ssignature_list <- as.list(signaturedf[, ssignatures])
ssignature_list <- lapply(as.list(ssignature_list), function(x){ # removing NAs and empty elements
  y <- x[!is.na(x)]; gsub("'| ", "", y[y != ""])
});

mycells$groups <- do.call(paste, c(mycells@meta.data[, idents, drop = FALSE], sep = "_"))
table(mycells$groups)
scells <- getsubset(sselect, mycells[[]], v = TRUE)
scells <- sample_grp(mycells@meta.data[scells, ], cname = 'groups', maxln = -3000, v = TRUE)
expdata <- expm1(mycells@assays$RNA@data[, scells]) * 100
genes <- rownames(expdata)[Matrix::rowMeans(expdata) > 1]
length(genes)
expdata <- as.matrix(expdata[genes, ])
resnamef <- paste0("metrics_per_", idents[1], ".csv")
if(file.exists(resnamef)){
  res <- read.csv(resnamef, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  res <- setrows(res)
}else{
  res <- data.frame(
    gene_name = paste0("'", genes),
    row.names = genes
  )
}
void <- list()
for(myclust in mygroups){
  myclust[2] <- ifelse(is.na(myclust[2]), "REST", myclust[2])
  mynameis <- paste0(myclust[1], "VS", myclust[2])
  cat(mynameis, "\n")
  if(!mynameis %in% colnames(res)){
    scells <- if(myclust[2] != "REST") getsubset(c("groups", myclust), mycells[[]], v = TRUE) else colnames(expdata)
    scells <- scells[scells %in% colnames(expdata)]
    annot <- mycells@meta.data[scells, ]
    annot$tmp <- ifelse(as.character(annot$groups) == myclust[1], yes = myclust[1], no = myclust[2])
    annot <- annot[order(annot$tmp), ]
    if(unique(annot$tmp)[1] != myclust[1]) annot <- annot[order(annot$tmp, decreasing = TRUE), ]
    tvar <- make_list(annot, "tmp", grouping = T)
    res$s2n <- gsea_metric(
      groups = tvar,
      mat = expdata[, rownames(annot)],
      metric = 'Signal2Noise',

      v = TRUE
    )
    colnames(res) <- gsub("^s2n$", mynameis, colnames(res))
    write.csv(res, file = resnamef)
  }
  void[[mynameis]] <- gsea_tests(
    res = res,
    gene_name = "gene_name",
    lfc.type = mynameis,
    gsea_file = ssignature_list,
    method = mymethod,
    myseed = 27,
    path = paste0(mymethod, '_', mynameis, "/"),
    plot_all = TRUE,
    v = TRUE
  )
}
mygseas <- rbindlist(lapply(names(void), function(x){
  y <- void[[x]]
  cbind(comparison = strsplit(basename(x), "_a1")[[1]][1], y)
}))
fwrite(mygseas, file = paste0("a1_gsea_summary_", gsub(" .*", "", Sys.time()), ".txt"), sep = "\t")

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
mysum[is.na(mysum)] <- 0
identnames = unlist(mygroups)
names(identnames) = unlist(mygroups)
annoc <- data.frame(
  Group = names(identnames), row.names = names(identnames)
)
annoc <- annoc[, sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
anncolist <- lapply(annoc, function(x) v2cols(x, gr.cols, v = TRUE) )
cr <- cc <- TRUE
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
  filename          = paste0("a1_gsea_summary_", names(plot_var), "_", gsub(" .*", "", Sys.time()), ".pdf"),
  width = 10, height = 7
);

# in case it failed
mysum <- tmysum[x$tree_row$labels[x$tree_row$order], x$tree_col$labels[x$tree_col$order]]
cr <- cc <- FALSE
