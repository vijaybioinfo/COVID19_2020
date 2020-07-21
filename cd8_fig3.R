#!/usr/bin/R

############
# Figure 3 #
############

# This script will create plots for figure 3 from cd8 data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd8/Figure_3'
setwdc(dirfig)

# Global variables
colsname <- "/home/ciro/covid19/info/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

#### 3.B UMAP per donor ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
sselect = list(c("orig.virus2", "CV"))
xaxes = c("UMAP_1", "UMAP_2", "RNA_snn_res.0.2", "orig.donor")
xaxes = c("UMAP_1", "UMAP_2", "orig.donor", "RNA_snn_res.0.2")

sselect = list(c("orig.virus2", "CV"), c("orig.donor", "P08"))
xaxes = c("UMAP_1", "UMAP_2", "RNA_snn_res.0.2", "origlib")

sselect = list(c("orig.virus2", "CV"), c("origlib", "CoVi17_8_M_24_9P_Gex"))

## Reading
mycells <- theObjectSavedIn(clustname)

## Operations
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
dfplot <- FetchData(
  object = mycells,
  var = xaxes,
  cells = scells
)
head(dfplot)

p <- ggplot(dfplot[scells, ], aes_string(x = xaxes[1], y =  xaxes[2], color = xaxes[3])) +
  geom_point(size = .3) +
  facet_wrap(facets = paste0("~", xaxes[4])) +
  labs(color = NULL) +
  mytheme +
  SetLegendPointsGG()

pdf(paste0("f3b_umap_", casefold(gsub("orig.", "", xaxes[4])), ".pdf"), width = 17, height = 15)
print(p)
dev.off()
pdf(paste0("f3b_umap_", casefold(gsub("orig.", "", xaxes[4])), "_blank.pdf"), width = 17, height = 15)
print(shut_it(p))
dev.off()

#### 3.C Cluster proportions per donor ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gr.cols["No", 1] <- "#8bd3fb"
gr.cols["Yes", 1] <- "#f091d4"
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.2'
dcolumn = 'orig.donor'
subdiv = "orig.hospital"
fname <- "f3c_traffic_donor"
sselect <- list(c(subdiv, "Yes", "No"), c("orig.virus2", "CV"))
corder = c('1', '0', '2', '3', '5')
roder = NULL
roder = c('P45', 'P24', 'P66', 'P49', 'P19', 'P26', 'P18', 'P40', 'P44', 'P61',
  'P12', 'P46', 'P20', 'P10', 'P29', 'P37', 'P07', 'P42', 'P22', 'P16', 'P01',
  'P25', 'P57', 'P17', 'P43', 'P06', 'P04', 'P05', 'P32', 'P64', 'P47', 'P15',
  'P03', 'P09', 'P08'
)

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

## Operations
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
# scells <- readfile("/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster1severity/YesvsNo/a1_rdata/g1.csv", stringsAsFactor = F)[,2]
annot <- remove.factors(mycells@meta.data[scells, ])
table(annot[, c(dcolumn, subdiv)])
annot[annot[, dcolumn] %in% "P0901", dcolumn] <- "P09" # combine P09 and P0901
table(annot[, c(dcolumn, subdiv)])
tvar <- make_list(annot, colname = dcolumn, col_objects = subdiv)
tvar <- reshape2::melt(lapply(tvar, function(x) paste0(unique(x), collapse = "-") ))
dgroup <- as.character(tvar[, 1])
names(dgroup) <- as.character(tvar[, 2])

ddfprops <- table(annot[, c(dcolumn, nres)])
# Proportions heatmap
matpct <- as.matrix(as.data.frame.matrix(prop.table(ddfprops, margin = 1)))
matpct <- round(matpct, 2)
matpct <- if(!is.null(roder)){
  suffix = "_supervised"
  matpct[roder, corder]
}else{
  suffix = "_unsupervised"
  matpct[, corder]
}
topz <- c(0, .7)
matpct[matpct > topz[2]] <- topz[2]; matpct[matpct < topz[1]] <- topz[1];
palettebreaks <- seq(from = 5, to = 70, by = 5)
annoc <- data.frame(
  Group = dgroup, row.names = names(dgroup)
)
annoc <- annoc[rownames(matpct), sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
anncolist <- lapply(annoc, function(x) v2cols(x, gr.cols, v = TRUE) )
mypalette <- colorRampPalette(colors = c("white", "red"), space = 'Lab')
mypalette <- colorRampPalette(colors = viridis::inferno(10), space = 'Lab')
mypalette <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 3, name = "YlGn"), space = 'Lab')
mypalette <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 3, name = "Greens"), space = 'Lab')
mypalette <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 3, name = "Blues"), space = 'Lab')
pdf(paste0("donor_heatmap", suffix, "_blues.pdf"), width = 10, height = 12, onefile = FALSE)
x <- NMF::aheatmap(
  x = matpct, annRow = annoc, annColors = anncolist,
  scale = 'none', Rowv = ifelse(is.character(roder), NA, TRUE), Colv = NA,
  col = mypalette(length(palettebreaks) - 1)
)
graphics.off()

#### 3.A Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.6'
cgroups = c('severity', "Mild", "NITU", "ITU")
orderit = c("NITU", "ITU")
cgroups = c('hospital', "No", "Yes")
orderit = -2
# 3, 2, 10, 6, 7, 5, 4, 13, 0, 1, 11, 9, 12, 8
sselect <- list(c("orig.severity", "ITU", "NITU", "Mild"))

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

## Operations
dfplot <- FetchData(
  object = mycells,
  cells = getsubset(sselect, mycells@meta.data, v = TRUE),
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster <- factormix(dfplot[, nres])
dfplot[, cgroups[1]] <- as.character(make_title(dfplot[, paste0("orig.", cgroups[1])]))
dfplot[, cgroups[1]] <- factor(dfplot[, cgroups[1]], make_title(cgroups[-1]))
table(dfplot[, cgroups[1]])
table(dfplot$orig.donor, dfplot[, cgroups[1]])

# ---- Start generating the figures
scells <- sample_grp(annot = dfplot, cname = cgroups[1], v = TRUE)
p <- ggplot(dfplot[scells, ], aes(x = UMAP_1, y =  UMAP_2, color = Cluster)) +
  geom_point() +
  facet_wrap(facets = ~severity) +
  labs(color = NULL) +
  mytheme +
  SetLegendPointsGG()

pdf(paste0("f3a_umap_", cgroups[1], ".pdf"), width = 15, height = 15)
print(p)
dev.off()
pdf(paste0("f3a_umap_", cgroups[1], "_blank.pdf"), width = 15, height = 15)
print(shut_it(p))
dev.off()

pp <- plot_pct(x = dfplot, groups = c(cgroups[1], "Cluster"), orderby = orderit, normalise = TRUE, return_table = TRUE)
p2 <- pp$plot + coord_flip() + mytheme
pdf(paste0("f3a_proportions_", cgroups[1], ".pdf"), width = 5, height = 15)
print(p2)
dev.off()
pdf(paste0("f3a_proportions_", cgroups[1], "_blank.pdf"), width = 5, height = 15)
print(shut_it(p2))
dev.off()

propdf <- pp$table
write.csv(propdf, file = paste0("f3a_proportions_", cgroups[1], ".csv"))

p <- plot_pct(x = dfplot, groups = c("Cluster", "severity"), type = "donut") +
  labs(fill = NULL) +
  blank_theme + theme(
    legend.position = "right",
    strip.text = element_text(face = 'bold', size = 10),
    axis.text.x = element_blank()
  )

pdf(paste0("sf3a_donut_", cgroups[1], ".pdf"), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0("sf3a_donut_", cgroups[1], "_blank.pdf"), width = 10, height = 10)
print(shut_it(p))
dev.off()

### 3.G Curtain plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
vconfigs <- list(
  vconfig = list(
    name = "f3g_curtain_cytotox",
    genes = c('GZMB', 'GZMH', 'GNLY', 'FASLG'),
    sfilter = list(c('RNA_snn_res.0.2',"1"), c('orig.hospital', 'No', 'Yes'), c('orig.virus2', 'CV')),
    ident = c('orig.hospital'),
    order = c('No', 'Yes')
  ),
  vconfig = list(
    name = "f3g_curtain_inflammatory",
    genes = c('CCL3', 'CCL4', 'TNF', 'CSF2', 'LTA', 'LTB'),
    sfilter = list(c('RNA_snn_res.0.2',"1"), c('orig.hospital', 'No', 'Yes'), c('orig.virus2', 'CV')),
    ident = c('orig.hospital'),
    order = c('No', 'Yes')
  ),
  vconfig = list(
    name = "f3g_curtain_tfs",
    genes = c('TBX21', 'BHLHE40', 'NFKB2', 'REL', 'FOS', 'JUNB'),
    sfilter = list(c('RNA_snn_res.0.2',"1"), c('orig.hospital', 'No', 'Yes'), c('orig.virus2', 'CV')),
    ident = c('orig.hospital'),
    order = c('No', 'Yes')
  )
)

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

## Operations
for(vconfig in vconfigs){
  scells <- getsubset(vconfig[['sfilter']], mycells[[]], v = TRUE)
  mycells$tmp <- do.call(paste, c(mycells[[]][, vconfig[['ident']], drop = FALSE], sep = "_"))
  if(any(vconfig[['order']] == "REST")){
    tvar <- as.character(mycells$tmp)
    mycells$Identity <- ifelse(!tvar %in% vconfig[['order']], "REST", tvar)
  }else{ mycells$Identity <- mycells$tmp }
  mycells$Identity <- factor(mycells$Identity, vconfig[['order']])
  table(mycells$tmp, mycells$Identity)

  # Dot-plots
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
  pdf(paste0(vconfig[['name']], ".pdf"), width = 4, height = 5)
  print(p)
  graphics.off()
  pdf(paste0(vconfig[['name']], "_blank.pdf"), width = 4, height = 5)
  print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
  graphics.off()
}

### 3.X Scatter plots cluster 1 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
prefix = "f2"
dscatter = "scatter"
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
mycname = "orig.hospital"
sselect = list(c("RNA_snn_res.0.2", "1"), c("orig.hospital", "No", "Yes"))
addaxisf = NULL
logit = TRUE
usedpcolor = TRUE # use double possitive as colour

# some genes
genesdf = data.frame(
  g1 = c('CCL3', 'TNF', 'CSF2', 'CSF2', 'TNF', 'TNF', 'LTA', 'CSF2', 'CSF2'),
  g2 = c('CCL4', 'CSF2', 'CCL3', 'CCL4', 'CCL3', 'CCL4', 'LTB', 'LTB', 'LTA'),
  stringsAsFactors = FALSE
)
# more genes
genesdf = data.frame(
  g1 = c('NFKB2', 'BHLHE40', 'BHLHE40'),
  g2 = c('BHLHE40', 'TBX21', 'RELA'),
  stringsAsFactors = FALSE
)

# more genes!
dscatter = "scatter_bhlhe40"
genesdf = data.frame(
  g1 = c('BHLHE40', 'BHLHE40', 'BHLHE40', 'BHLHE40', 'BHLHE40'),
  g2 = c('CSF2', 'TNF', 'CCL3', 'CCL4', 'GZMB'),
  stringsAsFactors = FALSE
)

# now a signature :0
logit = FALSE
usedpcolor = FALSE
dscatter = "scatter_apop"
sselect = list(c("RNA_snn_res.0.2", "0"), c("orig.hospital", "No", "Yes"))
addaxisf = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/signatures/25PCs/signatures.csv"
genesdf = data.frame(
  g1 = c('APOP_INH.Score'),
  g2 = c('APOP_ACT.Score'),
  stringsAsFactors = FALSE
)

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)
expdata <- expm1(mycells@assays$RNA@data) * 100; gc()
expdata <- cts2cpm(expdata)
if(!is.null(addaxisf)){
  addaxis <- readfile(addaxisf, row.names = 1, stringsAsFactors = FALSE)
  tvar <- unname(c(unlist(sselect), unlist(genesdf)))
  tvar <- tvar[tvar %in% colnames(addaxis)]
  # mycells@meta.data <- joindf(mycells@meta.data, addaxis[, tvar])
  if(any(!rownames(expdata) %in% colnames(addaxis)) && length(tvar) > 0){
    if(any(tvar %in% colnames(expdata))) tvar <- tvar[tvar %in% colnames(expdata)]
    expdata <- rbind(expdata, t(addaxis[colnames(expdata), tvar]))
    tailmat(expdata)
    expdata[tvar, ] <- abs(min(expdata[tvar, ])) + expdata[tvar, ]
  }
}

## Operations
ddf <- FetchData(
  object = mycells,
  vars = c(mycname, sapply(sselect, head, 1), unlist(genesdf, use.names = FALSE))
)
dir.create(dscatter)
scells <- getsubset(sselect, ddf, v = TRUE)
scells_global <- sample_grp(ddf[scells, ], cname = mycname, v = TRUE)
for(cl in as.character(unique(ddf[scells_global, mycname]))){
  firstfilt <- list(c(mycname, cl))
  scells <- getsubset(firstfilt, ddf[scells_global, ], v = TRUE)
  fname <- paste0(dscatter, "/", prefix, "_in_", casefold(gsub("orig.", "", mycname)), cl)
  for(i in 1:nrow(genesdf)){
    cat("Using", unlist(genesdf[i, ]), "\n")
    p <- get_densities(
      mat = expdata[, scells],
      genes = unlist(genesdf[i, ]),
      log2t = logit, cuof = 1, usedp = usedpcolor,
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

#### 3.F Volcano ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefix = "f3f_"
## Input
padjthr = 0.05
fcthr = .25
showgenes = c("CSF2", "TNF")
resname = "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2n3_sng_20p_25PCs_0.2R/comprs/cluster1severity/YesvsNo/results_YesvsNo_mastlog2cpm.csv"
trimmer = "c300"

## Reading
res <- readfile(resname, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
str(res)

## Operations
# res$log2FoldChange <- -1*res$log2FoldChange # in case you want to reverse it
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
showgenes <- unique(c(showgenes, bordering(modres[tvar, ], cnames = "log2FoldChange", n = 25)))
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
  fname <- paste0(prefix, "volcano_", basename(dirnamen(resname, 2)), "_trim", trimit)
  pdf(paste0(fname, ".pdf"), width = 10, height = 10)
  print(void)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 10, height = 10)
  print(shut_it(void))
  dev.off()
}

tvar <- modres[genes2plot, ]
genes <- rownames(tvar[order(-tvar$log2FoldChange), ])
mycells@meta.data$donor <- mycells@meta.data$orig.donor
mycells@meta.data$donor[mycells@meta.data$donor %in% "P0901"] <- "P09" # combine P09 and P0901
scells <- getsubset(list(c('orig.donor', '-void'), c('orig.hospital', 'No', 'Yes')), mycells@meta.data, v = TRUE)
fname <- paste0(prefix, "heatmap_", basename(dirnamen(resname, 2)))
source('/home/ciro/scripts/functions/pheatmapCorrection.R')
pdf(paste0(fname, "_alldegs.pdf"), width = 10, height = 12, onefile = FALSE)
custom_heatmap(
  object = mycells[, scells],
  rnames = genes,
  orderby = c("orig.hospital", "pca"),
  use_mean = "donor",
  sample_it = c(cname = "orig.hospital", maxln = "-1000"),
  scale_row = TRUE,
  categorical_col = c("orig.severity", "donor", "orig.hospital"),
  feature_order = TRUE,
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

### 3.X Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
vconfig = list(
  name = "f3g_violin_cluster1",
  genes = c('GZMB', 'GZMA', 'GNLY',
          c("CCL3", "TNF", "CSF2", "LTA", "CCL4", "LTB")),
  sfilter = list(c('RNA_snn_res.0.2', "1"), c('orig.hospital', 'No', 'Yes'), c('orig.virus2', 'CV')),
  ident = c('orig.hospital'),
  order = c('No', 'Yes')
)
vconfig = list(
  name = "f3_violin_cluster1",
  genes = c('GBP1', 'OAS1'),
  sfilter = list(c('RNA_snn_res.0.2', "1"), c('orig.hospital', 'No', 'Yes'), c('orig.virus2', 'CV')),
  ident = c('orig.hospital'),
  order = c('No', 'Yes')
)
vconfig = list(
  name = "f3_violin_cluster1_tfs",
  genes = c('TBX21', 'BHLHE40', 'NFKB2', 'REL', 'FOS', 'JUNB'),
  sfilter = list(c('RNA_snn_res.0.2',"1"), c('orig.hospital', 'No', 'Yes'), c('orig.virus2', 'CV')),
  ident = c('orig.hospital'),
  order = c('No', 'Yes')
)
vconfig = list(
  name = "f3_violin_cluster1_bhlhe40",
  genes = c('CSF2', 'TNF', 'CCL3', 'CCL4', 'GZMB'),
  sfilter = list(c('RNA_snn_res.0.2',"1"), c('orig.hospital', 'Yes')),
  ident = c('tag_BHLHE40'),
  order = c('BHLHE40+', 'BHLHE40-')
)

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)

## Operations
if(grepl("tag_", vconfig[['ident']])){
  genes <- gsub("tag_", "", vconfig[['ident']])
  void <- add_gene_tag(lgenes = genes, annot = mycells[[]], mat = mycells@assays$RNA@data, v = TRUE)
  head(void)
  mycells@meta.data <- joindf(mycells@meta.data, void)
}

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
table(ddf$tmp, ddf$Identity)

# Now individually... violins
dname <- paste0(vconfig[['name']])
dir.create(dname)
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

### S3.A Violin plots of signatures ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggpubr)
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
nres = 'RNA_snn_res.0.2'
signaturef = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/signatures/25PCs/signatures.csv"
ssigna_patterns = "sus2|interf|ty_guo|less_cul"
sselect = list(c("RNA_snn_res.0.2", "1"), c("orig.hospital", "No", "Yes"))
idents = "orig.hospital"
corder = c("No", "Yes")

## Reading
signaturedf <- readfile(signaturef, row.names = 1, stringsAsFactors = FALSE)

## Operations
tvar <- grepl(ssigna_patterns, colnames(signaturedf), ignore.case = TRUE)
ssignatures <- colnames(signaturedf)[tvar & grepl("Score", colnames(signaturedf))]
ssignatures
ddf <- signaturedf[, ssignatures]
colnames(ddf) <- ssignatures <- gsub(".Score", "", ssignatures)
head(ddf)
ddf <- joindf(ddf, FetchData(mycells, vars = c("UMAP_1", "UMAP_2", nres, idents)))
ddf[, idents] <- factor(ddf[, idents], corder)
unique(ddf[, idents])

scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
unique(ddf[scells, idents])
colby = c('pct', 'mean')[2]
for(ssig in ssignatures){
  cat(ssig, "\n")
  p <- violin(
    dat = ddf[scells, ],
    xax = idents,
    yax = ssig,
    dots = FALSE,
    colour_by = colby
  ) + stat_compare_means(method = "t.test")
  fname <- paste0("sf3a_", casefold(ssig), "_", colby)
  pdf(paste0(fname, ".pdf"), width = 5, height = 5)
  print(p)
  dev.off()
  pdf(paste0(fname, "_blank.pdf"), width = 5, height = 5)
  print(shut_it(p))
  dev.off()
}

### S3.B GSEA plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(pheatmap)
source('/mnt/BioHome/ciro/scripts/gsea/gsea_liger.R')
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD8T24_R1n2n3_sng_20p/clustering/zetInfo/clustCells25PCs_30Ks_0.06667JD.RData"
signaturef = "../Figure_2/supptable_signatures.csv"
ssigna_patterns = "sensus|interf|ty_guo|s_cullen"
idents = c("RNA_snn_res.0.2", "orig.hospital")
mymethod = "fgsea"
# Cluster 1
prefix = "cluster1"
sselect = list(c("RNA_snn_res.0.2", "1"), c("orig.hospital", "No", "Yes"))
mygroups = list(c("1_Yes", "1_No"))

## Reading
if(!exists("mycells")) mycells <- theObjectSavedIn(clustname)
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
resnamef <- paste0("metrics_per_", prefix, idents[1], ".csv")
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
      rnames = rownames(res),
      v = TRUE
    )
    colnames(res) <- gsub("^s2n$", mynameis, colnames(res))
    write.csv(res, file = resnamef)
  }
  void[[mynameis]] <- gsea_liger(
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
fwrite(mygseas, file = paste0("a1_gsea_summary_", prefix, gsub(" .*", "", Sys.time()), ".txt"), sep = "\t")

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

### 4.C Trajectory ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(monocle3)
## Input
clustname = "/home/ciro/large/covid19/results/trajectory/CD8T24_R1n2n3_sng_20p_25PCs_cv/object.rdata"

## Reading
cds <- theObjectSavedIn(clustname)

## Operations
p <- plot_cells(
  cds = cds,
  reduction_method = "UMAP",
  color_cells_by = "RNA_snn_res.0.2",
  group_cells_by = 'cluster',
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_branch_points = TRUE, # black circles
  label_roots = TRUE, # white circles
  label_leaves = TRUE, # gray circles
  graph_label_size = 4
) + theme_cowplot() + theme(legend.title = element_blank())
pdf(paste0("f3_trajectory_umap.pdf"), width = 8, height = 8)
print(p)
dev.off()
p2 <- shut_it(p, "text")
p2$layers <- p2$layers[1:2]
pdf(paste0("f3_trajectory_umap_blank.pdf"), width = 8, height = 8)
print(p2)
dev.off()
