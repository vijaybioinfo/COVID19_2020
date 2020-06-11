#!/usr/bin/R

############
# Figure 2 #
############

# This script will create plots for figure 2 from cd4 data

source('/home/ciro/scripts/functions/handy_functions.R')
source('/home/ciro/scripts/seurat/plotting.R')
source('/home/ciro/scripts/seurat/utilities.R')
library(Seurat)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/F2'
setwd(dirfig)

# Global variables
colsname <- "/home/ciro/covid19/info/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

#### 2.A Principal UMAP ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"

## Reading
if(clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
mycells@meta.data$Cluster <- factormix(mycells@meta.data[, nres])
p <- DimPlot(
  object = mycells,
  reduction = 'umap',
  group.by = "Cluster",
  label = TRUE
) + scale_y_reverse()
pdf('f2a_umap.pdf')
print(p)
dev.off()
pdf('f2a_umap_blank.pdf')
print(shut_it(p, lays = "ext"))
dev.off()

#### 2.B Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.6'

## Reading
if(clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
dfplot <- FetchData(
  object = mycells,
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster <- factormix(dfplot[, nres])
dfplot$virus <- make_title(dfplot$orig.virus2)
# dfplot$virus[dfplot$orig.virus2 == "HCV"] <- "HCV"
table(dfplot$virus)

# ---- Start generating the figures
scells <- sample_grp(annot = dfplot, cname = 'virus', v = TRUE)
p <- ggplot(dfplot[scells, ], aes(x = UMAP_1, y =  UMAP_2, color = Cluster)) +
  geom_point() +
  facet_wrap(facets = ~virus) +
  labs(color = NULL) +
  mytheme + scale_y_reverse() +
  SetLegendPointsGG()

pdf("f2b_umap_virus.pdf", width = 15, height = 15)
print(p)
dev.off()
pdf("f2b_umap_virus_blank.pdf", width = 15, height = 15)
print(shut_it(p))
dev.off()

p2 <- plot_pct(x = dfplot, groups = c("virus", "Cluster"), orderby = -1) + coord_flip() + mytheme
dfplot$tmp <- factor(dfplot$Cluster, rev(c('6', '0', '4', '7', '12', '11', '10', '1', '8', '9', '3', '5', '2')))
propdf <- table_pct(df = dfplot, cnames = rev(c("virus", "Cluster")))
write.csv(propdf, file = "f2b_proportions_virus.csv")

pdf("f2b_proportions_virus.pdf", width = 5, height = 15)
print(p2)
dev.off()
pdf("f2b_proportions_virus_blank.pdf", width = 5, height = 15)
print(shut_it(p2))
dev.off()

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
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
marknames = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/markers/38PCs_RNA_snn_res.0.6_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
nres = 'RNA_snn_res.0.6'
dgenes = c(
  'IL21', 'POU2AF1', 'CD200', 'BTLA', 'IFNG', 'IL2', 'IL3', 'TNF', 'CSF2',
  'IL17A', 'IL17F', 'CCR6', 'CTSH', 'IL4I1', 'IFIT3', 'IFI44L', 'ISG15',
  'MX2', 'OAS1', 'PRF1', 'GZMB', 'GNLY', 'NKG7', 'XCL1', 'XCL2', 'CCR7',
  'IL7R', 'TCF7', 'MKI67', 'CDK1', 'CCNB1', 'UBE2C'
)
suffix = "markers_shared"
ntop = 200
identnames = c(
  "0" = "TFH", "6" = "TFH", "7" = "TFH",
  "1" = "TH1 PolyF", "10" = "TH1 PolyF",
  "9" = "TH17", "2" = "THIFNR",
  "5" = "THIFNR",
  "4" = "CTL", "8" = "CTL", "11" = "CTL-like",
  "3" = "Central Memory",
  "12" = "Cycling"
)

## Reading
if(clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
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
source('/home/ciro/scripts/functions/pheatmapCorrection.R')
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
lapply(genesl, function(x) x[gsub("'", "", x) %in% dgenes] )

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
fname <- "f2d_curtain"
pdf(paste0(fname, ".pdf"), width = 10, height = 13)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 10, height = 13)
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()

clustnamebk <- clustname

#### S2.C Volcano ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
padjthr = 0.05
fcthr = .25
showgenes = NULL
resname = "/home/ciro/large/covid19/results/dgea/CD4T6_R1n2_sng2_25p/comprs/incv_clusters/0vs6/results_0vs6_mastlog2cpm.csv"
trimmer = "c170"

## Reading
res <- readfile(resname, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
str(res)
colnames(res) <- gsub("FLU~HPIV~MPV", "NCV", colnames(res))

## Operations
# res$log2FoldChange <- -1*res$log2FoldChange # for 0 vs 6
dtype <- sub("Bmean", "", colnames(res)[grepl("Bmean", colnames(res))])
grps <- unlist(strsplit(gsub("results_|_mast.*", "", basename(resname)), "vs"))
grp2 <- grps[1]
tvar <- paste0(grps, "_mean", dtype)
modres <- data.frame(res[!is.na(res$padj), ], stringsAsFactors = F, check.names = FALSE)
modres$gene <- gsub("'", "", modres$gene)
modres$Mean <- round(log2(ifelse(modres$group == grp2, modres[, tvar[1]], modres[, tvar[2]]) + 1), 1)
modres$pcts <- modres$pct_diff
genes2plot <- mysignames <- getDEGenes(modres, pv = padjthr, fc = fcthr, gene_name = "gene", further = NULL, v = TRUE)
tvar <- rownames(modres)[modres[, grep("minExp", colnames(modres), value = TRUE)]]
genes2plot <- genes2plot[genes2plot %in% tvar]
modres$degs <- "Not_significant"
modres[modres$gene %in% genes2plot, ]$degs <- "DEG"
modres[!modres$gene %in% genes2plot, ]$Mean <- NA
modres[!modres$gene %in% genes2plot, ]$pcts <- 0
tvar <- cosmy(genes2plot, patties = "^rps|^rpl|^mt-|rp[0-9]{1,}-|^linc")
if(is.null(showgenes)) showgenes <- bordering(modres[tvar, ], cnames = "log2FoldChange", n = 25)
source('/home/ciro/scripts/functions/volcano_variant_color.R')
tvar <- -log10(modres$padj); tvar[is.infinite(tvar)] <- max(tvar[is.finite(tvar)]); summary(tvar)
for(trimit in trimmer){
  trimit <- ifelse(
    max(tvar[is.finite(tvar)])<as.numeric(gsub("c", "", trimit)),
    paste0("c", round(max(tvar[is.finite(tvar)]))),
    trimit
  )
  void <- try(volplot(
    modres,
    pvalth = padjthr,
    lfcth = fcthr,
    pvaltype = 'padj',
    lfctype = 'log2FoldChange',
    col_feature = "Mean",
    size_feature = "pcts",
    gene_name = 'gene',
    group = "degs",
    check_genes = list(text = parse_ens_name(showgenes)),
    titl = paste0("'", paste0(gsub("_", " ", grps), collapse = "' vs '"), "'"),
    return_plot = TRUE,
    clipp = trimit,
    v = TRUE
  ) + labs(size = "%", title = NULL))
  fname <- paste0(prefix, "volcano_", basename(dirname(resname)), "_trim", trimit)
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(void)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 10, height = 10)
  print(shut_it(void))
  dev.off()
}
clustnamebk <- clustname

### 2.E Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
allconfigs <- list(
  vconfig = list(
    name = "f2e_th1",
    genes = c('IFNG', 'IL2'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(10, 5),
    order = c('1', '10', 'REST')
  ),
  vconfig = list(
    name = "sf2d_th1",
    genes = c("TNF", "CSF2", "IL3"),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('1', '10', 'REST')
  ),
  vconfig = list(
    name = "f2e_tfh",
    genes = c('IL21', 'POU2AF1'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(10, 5),
    order = c('0', '6', '7', 'REST')
  ),
  vconfig = list(
    name = "f2e_tfh",
    genes = c('CD200', 'BTLA', 'CXCL13'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('0', '6', '7', 'REST')
  ),
  vconfig = list(
    name = "f2e_th17",
    genes = c('IL17A', 'CCR6'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(10, 5),
    order = c('2', '9', 'REST')
  ),
  vconfig = list(
    name = "sf2d_th17",
    genes = c('IL17F', 'CTSH', 'IL4I1'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('2', '9', 'REST')
  ),
  vconfig = list(
    name = "f2d_violin",
    genes = c('IFIT3', 'IFI44L', 'ISG15'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('2', '5', 'REST')
  ),
  vconfig = list(
    name = "f2d_violin_ctl",
    genes = c('PRF1', 'GZMB', 'GNLY'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12)),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(15, 5),
    order = c('4', '8', 'REST')
  )
)

## Reading
if(clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
for(vconfig in allconfigs){
  cat(vconfig[['name']], "\n")
  scells <- getsubset(vconfig[['sfilter']], mycells@meta.data, v = TRUE)
  myvars <- c(vconfig[['genes']], vconfig[['ident']])
  ddf <- FetchData(mycells, vars = myvars, cells = scells)
  ddf$tmp <- do.call(paste, c(ddf[, vconfig[['ident']], drop = FALSE], sep = "_"))
  if(any(vconfig[['order']] == "REST")){
    tvar <- as.character(ddf[, vconfig[['ident']]])
    ddf$Identity <- ifelse(!tvar %in% vconfig[['order']], "REST", tvar)
  }else{ ddf$Identity <- ddf$tmp }
  ddf$Identity <- factor(ddf$Identity, vconfig[['order']])

  # Now individually... nightmare
  dname <- paste0(vconfig[['name']])
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
    fname <- paste0(vconfig[['name']], "/", vconfig[['genes']][i])
    pdf(paste0(fname, ".pdf"), width = 5, height = 5)
    print(p)
    dev.off()
    pdf(paste0(fname, "_blank.pdf"), width = 5, height = 5)
    print(shut_it(p))
    dev.off()
  }
}
clustnamebk <- clustname

#### 2.F Signatures ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.6'
signaturef = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/signatures/38PCs/signatures.csv"
ssigna_patterns = "tfh|cytotox|th17|interf|cyclying"
couls = c("#fffeee", "#ffe080", "#ffc100", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000")

## Reading
if(clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
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
  ddf$Signature <- ddf[, ssig]
  p <- ggplot(ddf, aes(x = UMAP_1, y = UMAP_2, color = Signature)) +
    geom_point(size = 0.3) +
    scale_colour_gradientn(
      colours = couls
    ) +
    theme(legend.position = "right") + labs(color = NULL) + scale_y_reverse()
  fname <- paste0("f2f_", ssig)
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
  pdf(paste0("sf2e_", ssig, "_violin_mean.pdf"), width = 12, height = 8)
  print(p)
  dev.off()
}
clustnamebk <- clustname

ssignatures_list <- vlab_signatures[grepl(ssigna_patterns, names(vlab_signatures))]
str(ssignatures_list)
ssignatures_list <- lapply(ssignatures_list, function(x) x[x %in% rownames(mycells)] )
str(ssignatures_list)
final_ssignatures_list <- lapply(names(ssignatures_list), function(x){
  y <- ssignatures_list[[x]]
  z <- vlab_signatures[[x]]
  c(y, "Not found", z[!z %in% y])
})
str(final_ssignatures_list)
ddf <- vlist2df(ssignatures_list)
write.csv(ddf, file = "sup_table_12_gene_signatures.csv")

#### S2.A QC Violin plots ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = c(
  "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData",
  "/home/ciro/large/covid19/results/clustering/CD4T24_R1n2_sng_25p/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData"
)
low.filt.cells = c(800, 1500, -Inf)
high.filt.cells = c(4400, 20000, 10)
subs.names = c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
names(low.filt.cells) <- subs.names
names(high.filt.cells) <- subs.names

## Reading
mycells6 <- theObjectSavedIn(clustname[1])
mycells24 <- theObjectSavedIn(clustname[2])
mycellscom <- merge(mycells6, mycells24)
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
