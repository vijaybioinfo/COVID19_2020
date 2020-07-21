#!/usr/bin/R

############
# Figure 2 #
############

# This script will create plots for figure 2 from cd8 data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')
source('/home/ciro/scripts/seurat/utilities.R')
library(Seurat)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd8/Figure_2'
setwdc(dirfig)

# Global variables
colsname <- "/home/ciro/covid19/info/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

#### 2.A Principal UMAP ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
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
pdf('f2a_umap.pdf')
print(p)
dev.off()
pdf('f2a_umap_blank.pdf')
print(shut_it(p, lays = "ext"))
dev.off()

clustnamebk <- clustname

#### 2.B Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.2'

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

## Operations
dfplot <- FetchData(
  object = mycells,
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster <- factormix(dfplot[, nres])
dfplot$virus <- make_title(dfplot$orig.virus2)
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

source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')
pp <- plot_pct(x = dfplot, groups = c("virus", "Cluster"), orderby = "CV", normalise = TRUE, return_table = TRUE)
propdf <- pp$table
write.csv(propdf, file = "f2b_proportions_virus.csv")

p2 <- pp$plot + coord_flip() + mytheme
pdf("f2b_proportions_virus.pdf", width = 5, height = 15)
print(p2)
dev.off()
pdf("f2b_proportions_virus_blank.pdf", width = 5, height = 15)
print(shut_it(p2))
dev.off()

p <- plot_pct(x = dfplot, groups = c("Cluster", "virus"), type = "donut") +
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

table(dfplot[, c("origlib", "orig.project_id")])
scells <- getsubset(c("orig.virus2", "CV"), dfplot, v = TRUE)
table(dfplot[scells, c("origlib", "orig.project_id")])
table(dfplot[, c("origlib", "orig.virus2")])
table(dfplot[scells, c("origlib", "orig.project_id")])
p <- plot_pct(x = dfplot[, ], groups = c("Cluster", "orig.project_id"), normalise = FALSE) +
  labs(fill = NULL) + coord_flip() + mytheme

pdf("sf2a_bar_project_id.pdf", width = 5, height = 10)
print(p)
dev.off()
pdf("sf2a_bar_project_id_blank.pdf", width = 4, height = 10)
print(shut_it(p))
dev.off()

clustnamebk <- clustname

#### S2.A QC Violin plots ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
low.filt.cells = c(800, 1500, -Inf)
high.filt.cells = c(4400, 20000, 10)
subs.names = c('nFeature_RNA', 'nCount_RNA', 'percent.mt')
names(low.filt.cells) <- subs.names
names(high.filt.cells) <- subs.names

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)
tvar <- reshape2::melt(table(mycells$orig.stim_time, mycells$origlib))
tvar <- tvar[tvar[, 3] > 0, ]
tvar <- tvar[order(tvar[, 1]), ]
mycells$Library <- factor(mycells$origlib, levels = unique(as.character(tvar[, 2])))
table(mycells$Library)
table(mycells@meta.data[mycells@meta.data$orig.donor == "P08", c("origlib", "orig.donor")])

## Operations
mycols <- rep("#2e82ff", length(table(mycells$Library)))
names(mycols) <- names(table(mycells$Library))
for(qcvar in subs.names[1]){
  p <- VlnPlot(
    object = mycells,
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

#### 2.C Cluster markers ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
marknames = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/markers/25PCs_RNA_snn_res.0.2_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05_CPM.csv"
nres = 'RNA_snn_res.0.2'
dgenes = c(
  'MX1', 'IFIT3', 'IFI6', 'IRF7',
  'CRTAM',
  'IL2RA', 'PCNA',
  'IFNG', 'CCL4', 'XCL2', 'XCL1', 'LDHA', 'CCL3', 'TNF',
  # 'APOBEC3C', 'APOBEC3G',
  'CD27', 'CCR7', 'SELL', 'IL7R',
  'ZNF683'
)
suffix = "markers"
ntop = 200
identnames = c("1", "0", "2", "4", "3", "5", "6", "7")[c(-8)]
names(identnames) = identnames
sselect = c("RNA_snn_res.0.2", "-7")

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)
mygenes <- readfile(marknames, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

## Operations
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
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
  object = mycells[, scells],
  rnames = genes,
  orderby = "Cluster",
  use_mean = "Cluster",
  sample_it = c(cname = "Cluster", maxln = "-1000"),
  scale_row = TRUE,
  categorical_col = c("Cluster", "Name"),
  feature_order = TRUE,
  couls = NULL,
  hcouls = rev(c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c")),#c('yellow', 'black', 'blue'),
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
mycells$tmp <- factor(mycells$Cluster, unname(identnames))
p <- DotPlot(
  object = mycells[, scells],
  features = dgenes,
  group.by = "tmp",
  cols = c('#fff4ba', '#ff0000'),
  col.min = -1.5, col.max = 1.5
) + coord_flip() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold.italic"),
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)
fname <- "f2c_curtain"
pdf(paste0(fname, ".pdf"), width = 8, height = 13)
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = 8, height = 13)
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()

clustnamebk <- clustname

### X.X Signature calculation ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
cpatterns = c("RNA_snn_res.0.2", "orig.hospital", "RNA_snn_res.hospital")
out_dir = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/"
pcs.comp = 25
norm_type = "RNA_snn_res"; verb = TRUE
rdims = list(umap = c('UMAP_1', 'UMAP_2'))

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

## Operations
tvar <- cpatterns[cpatterns %in% colnames(mycells[[]])]
mycells$RNA_snn_res.hospital = do.call(paste, c(mycells@meta.data[, tvar, drop = FALSE], sep = "_"))
mycells$RNA_snn_res.hospital[is.na(mycells$RNA_snn_res.hospital)] <- NA # make sure asthma is correct
mycells$RNA_snn_res.hospital[mycells$RNA_snn_res.hospital == "void"] <- NA
table(mycells$RNA_snn_res.hospital)

load("/home/ciro/vdv/vijay_lab/data/signatures_cd8covid.rdata")
load("/mnt/BioAdHoc/Groups/vd-vijay/vijay_lab/data/signatures.rdata")
helped_signatures <- read.csv("/home/ciro/covid19/info/tcell_cd8story_signatures_helped_2020_06_24.csv", stringsAsFactor = FALSE)
helped_signatures <- lapply(as.list(helped_signatures), function(x){ # removing NAs and empty elements
  y <- x[!is.na(x)]; gsub("'| ", "", y[y != ""])
});
source("/home/ciro/scripts/functions/gene_name_convertion.R")
helped_signatures <- lapply(helped_signatures, mouse2human)
ssignatures_list <- c(signatures_cd8covid, vlab_signatures, helped_signatures)
ssignatures_list <- ssignatures_listt <- ssignatures_list[grepl(ssigna_patterns, names(ssignatures_list))]
source("/home/ciro/scripts/gsea/signatures_seurat.R")
ssignatures_list <- clean_feature_list(mat = mycells@assays$RNA@data, features = ssignatures_list, filterby = "p2", v = TRUE)
str(ssignatures_list)
ssignatures_list <- lapply(ssignatures_list, function(x) x[x %in% rownames(mycells)] )
str(ssignatures_list)
dgenes <- ssignatures_list[['dixhaust_consensus2_vj']]
final_ssignatures_list <- lapply(names(ssignatures_list), function(x){
  y <- ssignatures_list[[x]]
  z <- ssignatures_listt[[x]]
  c(y, "Not found", z[!z %in% y])
})
str(final_ssignatures_list)
ddf <- vlist2df(final_ssignatures_list)
names(ddf) <- names(ssignatures_list)
write.csv(ddf, file = "supptable_signatures.csv")

ssignatures_list <- list(
  apop_inh = c('MCL1', 'BCL2L1', 'BIRC3', 'BIRC2', 'BCL2A1', 'NFKB1', 'RELA'),
  apop_act = c('BAX', 'BCL2L11', 'FAS', 'CYCS', 'LMNA', 'NFKBIA')
)
ssignatures_list <- clean_feature_list(mat = mycells@assays$RNA@data, features = ssignatures_list, filterby = "p2", v = TRUE)

void <- try(
  signature_scoring(
    object = mycells,
    prefix = paste0(out_dir, "signatures/", pcs.comp, 'PCs/'),
    lsignatures = rev(ssignatures_list),
    confounders = cpatterns,
    reductions = rdims,
    v = verb
  )
)

#### 2.F Signatures - violins and umaps ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.2'
signaturef = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/signatures/25PCs/signatures.csv"
ssigna_patterns = "sus2|interf|virusres|ty_guo|cycl|hallmark|less_cul"
couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)
signaturedf <- readfile(signaturef, row.names = 1, stringsAsFactors = FALSE)

## Operations
tvar <- grepl(ssigna_patterns, colnames(signaturedf), ignore.case = TRUE)
ssignatures <- colnames(signaturedf)[tvar & grepl("Score", colnames(signaturedf))]
ssignatures
ddf <- signaturedf[, ssignatures]
colnames(ddf) <- ssignatures <- gsub(".Score", "", ssignatures)
ddf <- joindf(ddf, FetchData(mycells, vars = c("UMAP_1", "UMAP_2", nres)))
ddf$Cluster <- factormix(ddf[, nres])
ddf <- ddf[order(ddf$Cluster), ]
str(ddf)
sapply(ddf, function(x) sum(is.na(x)) )

for(ssig in ssignatures[4]){
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
  fname <- paste0("f2f_", casefold(ssig))
  pdf(paste0(fname, "_umap_x.pdf"), width = 8, height = 8)
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
  pdf(paste0("sf2e_", casefold(ssig), "_violin_mean.pdf"), width = 12, height = 8)
  print(p)
  dev.off()
}
clustnamebk <- clustname

#### 2.X Exhaustion ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
vconfig = list(
  name = "f2_exhaustion",
  genes = "supptable_signatures.csv",
  sfilter = list(c('RNA_snn_res.0.2', 0:7)),
  ident = c('RNA_snn_res.0.2'),
  ncols = 2, size = c(10, 13),
  order = c('1', 'REST')
)
vconfig = list(
  name = "f2_exhaustion",
  genes = c('TOX', 'HAVCR2', 'LAG3', 'PDCD1'),
  sfilter = list(c('RNA_snn_res.0.2', 0:7)),
  ident = c('RNA_snn_res.0.2'),
  ncols = 1, size = c(5, 5),
  order = c('1', 'REST')
)

## Reading
signaturedf <- readfile(signaturef, row.names = 1, stringsAsFactors = FALSE)
genesdf <- readfile(vconfig$genes, stringsAsFactor = FALSE, row.names = 1)

## Operations
vconfig$genes = genesdf[, 'dixhaust_consensus2_vj']
vconfig$genes <- vconfig$genes[!is.na(vconfig$genes)]
cat(vconfig[['name']], "\n")
scells <- getsubset(vconfig[['sfilter']], mycells@meta.data, v = TRUE)
myvars <- c(vconfig[['genes']], vconfig[['ident']])
ddf <- FetchData(mycells, vars = myvars, cells = scells)
ddf$tmp <- do.call(paste, c(ddf[, vconfig[['ident']], drop = FALSE], sep = "_"))
if(any(vconfig[['order']] == "REST")){
  tvar <- as.character(ddf[, vconfig[['ident']]])
  ddf$Identity <- ifelse(!tvar %in% vconfig[['order']], "REST", tvar)
}else{ ddf$Identity <- ddf$tmp }
mycells$Cluster <- ddf$Identity <- factor(ddf$Identity, vconfig[['order']])
table(ddf$Identity)

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
  fname <- paste0(dname, "/", vconfig[['genes']][i])
  pdf(paste0(fname, ".pdf"), width = 5, height = 5)
  print(p)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 5, height = 5)
  print(shut_it(p))
  dev.off()
}

p <- DotPlot(
  object = mycells[, scells],
  features = vconfig[['genes']],
  group.by = "Cluster",
  cols = c('#fff4ba', '#ff0000'),
  col.min = -1.5, col.max = 1.5
) + coord_flip() +
  theme(
    axis.text.y = element_text(size = 13, face = "bold.italic"),
    axis.ticks.x = element_blank()
  ) + labs(y = NULL, x = NULL)
fname <- paste0(dname, "/f2_curtain")
pdf(paste0(fname, ".pdf"), width = vconfig[['size']][1], height = vconfig[['size']][2])
print(p)
graphics.off()
pdf(paste0(fname, "_blank.pdf"), width = vconfig[['size']][1], height = vconfig[['size']][2])
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()



#### 2.E Violins ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
allconfigs <- list(
  vconfig = list(
    name = "f2e_cluster4",
    genes = c('IFNG', 'XCL1', 'XCL2', 'CCL3', 'CCL4', 'TNF'),
    sfilter = list(c('RNA_snn_res.0.2', 0:7)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 2, size = c(6, 10),
    order = c('4', 'REST')
  ),
  vconfig = list(
    name = "f2e_cluster1",
    genes = c('GZMB', 'GZMA', 'GZMH', 'GNLY'),
    sfilter = list(c('RNA_snn_res.0.2', 0:7)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 1, size = c(5, 20),
    order = c('1', 'REST')
  ),
  vconfig = list(
    name = "f2e_cluster6",
    genes = c('ZNF683'),
    sfilter = list(c('RNA_snn_res.0.2', 0:7)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 1, size = c(5, 5),
    order = c('6', 'REST')
  ),
  vconfig = list(
    name = "f2_exhaustion",
    genes = c('TOX', 'HAVCR2', 'LAG3', 'PDCD1'),
    sfilter = list(c('RNA_snn_res.0.2', 0:7)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 1, size = c(5, 5),
    order = c('1', 'REST')
  ),
  vconfig = list(
    name = "f2_cluster1_mx1_oas1",
    genes = c('MX1', 'OAS1'),
    sfilter = list(c('RNA_snn_res.0.2', 0:7)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 1, size = c(5, 5),
    order = c('1', 'REST')
  ),
  vconfig = list(
    name = "f2_cluster2_pcna_mxm7",
    genes = c('PCNA', 'MCM7'),
    sfilter = list(c('RNA_snn_res.0.2', 0:7)),
    ident = c('RNA_snn_res.0.2'),
    ncols = 1, size = c(5, 5),
    order = c('2', 'REST')
  )
)

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

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


### 2.E Scatter proportion changes ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
dscatter = "scatter"
nres = c("RNA_snn_res.0.2", "4")
genesdf = data.frame(
  g1 = c("XCL1", "IFNG", "IFNG"),
  g2 = c("XCL2", "TNF", "CCL4"),
  stringsAsFactors = FALSE
)

## Reading
if(clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
annot <- mycells@meta.data
expdata <- expm1(mycells@assays$RNA@data) * 100; gc()
expdata <- cts2cpm(expdata)
annot = annot[colnames(expdata), ]

## Operations
ddf <- FetchData(
  object = mycells,
  vars = c(nres[1], unlist(genesdf, use.names = FALSE))
)
ddf$Cluster <- ifelse(ddf[, nres[1]] %in% nres[-1], paste0(nres[-1], collapse = "n"), "REST")
table(ddf$Cluster)
dir.create(dscatter)
scells_global <- sample_grp(ddf, cname = "Cluster", v = TRUE)
for(cl in unique(ddf$Cluster)){
  firstfilt <- list(c("Cluster", cl))
  scells <- getsubset(firstfilt, ddf[scells_global, ], v = TRUE)
  fname <- paste0(dscatter, "/f2e_in_cluster", cl, "scatter")
  for(i in 1:nrow(genesdf)){
    cat("Using", unlist(genesdf[i, ]), "\n")
    p <- get_densities(
      mat = expdata[, scells],
      genes = unlist(genesdf[i, ]),
      log2t = TRUE, cuof = 1, usedp = TRUE,
      return_plot = TRUE, v = !TRUE
    )$scatter
    pdf(paste0(fname, "_", paste0(unlist(genesdf[i, ]), collapse = "_"), ".pdf"))
    print(p + viridis::scale_color_viridis(option = "magma"))
    dev.off()
    pdf(paste0(fname, "_", paste0(unlist(genesdf[i, ]), collapse = "_"), "_blank.pdf"))
    print(shut_it(p, lays = "text|dens") + viridis::scale_color_viridis(option = "magma"))
    dev.off()
  }
}

#### S2.C GSEA ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table)
library(pheatmap)
## Input
gseaname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/gsea/25PCs/a1_gsea_summary_2020-06-26.txt"
corder = c("1", "0", "2", "4", "3", "5", "6", "7")[-8]
corder = as.character(0:6)
rorder = c(
  "dixhaust_consensus2_vj", "type1n2_interferon_broadreact", "tcell_helpless_cullen",
  "tcell_cytotoxicity_guo", "hallmarkglycolisis_signature_gsea", "cell_cycling_best"
)
cr <- cc <- FALSE;
identnames <- paste0("C", as.character(0:7))
names(identnames) <- gsub("C", "", identnames)

## Reading
mygseas <- readfile(gseaname, header = 1, stringsAsFactors = FALSE)

## Operations
head(mygseas)
thistab <- mygseas[mygseas$padj <= 0.05, ]
round2 <- function(x) round(x, digits = 2)
nlog10 <- function(x) -log10(x)
plot_var <- list(
  padj = list(title = "Adjuste P-value", transf = nlog10, thresh = 0.05),
  NES = list(title = "Normalized Enrichment Score", transf = round2, thresh = -Inf)
)[2] # Select type

mysum <- dcast.data.table(data.table(thistab), comparison ~ pathway, value.var = names(plot_var))
mysum <- plot_var[[1]]$transf(data.frame(mysum[, -1], stringsAsFactors = FALSE, row.names = mysum$comparison))
tmysum <- mysum <- t(as.matrix(mysum))
if(is.null(corder)) mysum[is.na(mysum)] <- 0 else mysum <- mysum[, corder]
if(is.null(corder)) mysum[is.na(mysum)] <- 0 else mysum <- mysum[rorder, ]
# mysum[mysum <= plot_var[[1]]$transf(plot_var[[1]]$thresh)] <- 0
annoc <- data.frame(
  Cluster = names(identnames), Identity = unname(identnames), row.names = names(identnames)
)
annoc <- annoc[corder, sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
anncolist <- lapply(annoc, function(x) v2cols(x, gr.cols, v = TRUE) )
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
  filename          = paste0("f2c_gsea_summary.pdf"),
  width = 10, height = 7
);

# in case it failed
mysum <- tmysum[x$tree_row$labels[x$tree_row$order], x$tree_col$labels[x$tree_col$order]]
cr <- cc <- FALSE
