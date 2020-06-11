#!/usr/bin/R

############
# Figure 3 #
############

# This script will create plots for figure 3 from cd4 data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/F3'
setwd(dirfig)

# Global variables
colsname <- "/home/ciro/covid19/info/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

#### S3.X Volcano ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
padjthr = 0.05
fcthr = .25
showgenes = c('PRF1', 'GZMB', 'ZBTB32', 'ZBED2', 'HAVCR2', 'LAG3', 'PDCD1', 'TIGIT', 'CD70', 'DOK5', 'DUSP4')
resname = "/home/ciro/large/covid19/results/dgea/CD4T6_R1n2_sng2_25p/comprs/incv_clusters/0vs6/results_0vs6_mastlog2cpm.csv"
prefix = "sf3a1_"
trimmer = c("c170", "c300", "c321", "c325")[3]
resname = "/home/ciro/large/covid19/results/dgea/CD4T6_R1n2_sng2_25p/comprs/incv_clusters/7vs6/results_7vs6_mastlog2cpm.csv"
prefix = "sf3a2_"
trimmer = c("c170", "c300", "c321", "c325")

### 3.X Curtain plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
vconfig = list(
  name = "f3b_curtain",
  genes = c('TIGIT', 'LAG3', 'HAVCR2', 'PDCD1', 'DUSP4', 'CD70', 'DOK5'),
  sfilter = list(c('RNA_snn_res.0.6', '0', '6', '7'), c('orig.virus2', 'CV')),
  ident = c('RNA_snn_res.0.6'),
  order = c('6', '0', '7')
)
vconfig = list(
  name = "f3c_curtain",
  genes = c('SLAMF7', 'CD72', 'GPR18'),
  sfilter = list(c('orig.virus2', 'CV')),
  ident = c('RNA_snn_res.0.6'),
  order = c('4', '8', 'REST')
)
vconfig = list(
  name = "f3d_curtain",
  genes = c('SLAMF7', 'CD72', 'GPR18'),
  sfilter = list(c('orig.virus2', 'CV')),
  ident = c('RNA_snn_res.0.6'),
  order = c('4', '8', 'REST')
)
vconfig = list(
  name = "f3d2_curtain",
  genes = c('HOPX', 'ZEB2', 'CD72', 'GPR18', 'SLAMF7'),
  sfilter = list(c('orig.virus2', 'CV')),
  ident = c('RNA_snn_res.0.6'),
  order = c('4', '8', 'REST')
)

## Reading
mycells <- theObjectSavedIn(clustname)

## Operations
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
pdf(paste0(vconfig[['name']], ".pdf"), width = 5, height = 5)
print(p)
graphics.off()
pdf(paste0(vconfig[['name']], "_blank.pdf"), width = 5, height = 5)
print(shut_it(p) + theme(legend.position = "right", legend.title = element_blank()))
graphics.off()

### 3.X Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
allconfigs <- list(
  vconfig = list(
    name = "f3b_violins",
    genes = c('ZBED2', 'ZBTB32'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '6', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(8, 4),
    order = c('6', '0', '7')
  ),
  vconfig = list(
    name = "f3d_violins",
    genes = c('ZEB2', 'HOPX'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(8, 4),
    order = c('4', '8', 'REST')
  ),
  vconfig = list(
    name = "sf3b_violins",
    genes = c('TIGIT', 'LAG3', 'HAVCR2', 'PDCD1', 'DUSP4', 'CD70', 'DOK5'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '6', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 4, size = c(3, 21),
    order = c('6', '0', '7')
  ),
  vconfig = list(
    name = "sf3c_violins",
    genes = c('PRF1', 'GZMB'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '6', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 1, size = c(8, 4),
    order = c('6', '0', '7')
  ),
  vconfig = list(
    name = "sf3e_violins",
    genes = c('SLAMF7', 'CD72', 'GPR18'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 3, size = c(12, 4),
    order = c('4', '8', 'REST')
  ),
  vconfig = list(
    name = "sf3f_violins",
    genes = c('CCL3', 'CCL4', 'CCL5', 'XCL1', 'XCL2'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 1, size = c(16, 4),
    order = c('4', '8', 'REST')
  )
)[c(4:5)]

## Reading
mycells <- theObjectSavedIn(clustname)

## Operations
sapply(allconfigs, function(x) all(x$genes %in% rownames(mycells)) )
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

### 3.X Co-expression/double positive proportions ###%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
fedata = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
fannot = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/metadata_38PCs_30Ks_0.06667JD.RData"
genes = c("PRF1", "GZMB")
selectss = list(c("RNA_snn_res.0.6", "0", "4", "6" , "7", "8")[c(1,2)])
genes = c("XCL1", "XCL2")
selectss = list(c("RNA_snn_res.0.6", "4", "8", "11")[c(1,2)])

## Reading
mycells <- readfile(fedata)
if(class(mycells) == "saver"){
  datatype <- "saver"
  expdata <- mycells$estimate; gc()
}
if(casefold(class(mycells)) == "seurat"){
  datatype <- "cpm"
  expdata <- expm1(mycells@assays$RNA@data) * 100; gc()
  annot <- mycells@meta.data
}else{
  annot <- readfile(fannot, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)
}

## Operations
expdata <- cts2cpm(expdata[, getsubset(c("RNA_snn_res.0.6", "0", "4", "6" , "7", "8"), annot, v = TRUE)])
dgenes <- remove.factors(create_pairs(data.frame(g = genes, stringsAsFactors = FALSE))[, 2:1])

scells <- getsubset(selectss, annot, v = TRUE)
for(i in 1:nrow(dgenes)){
  cat(commas(dgenes[i, ]), "\n")
  suffy <- paste0(datatype, "_", summary_subset(selectss))
  p <- get_densities(
    mat = expdata[unlist(dgenes[i, ]), scells],
    genes = unlist(dgenes[i, ]),
    log2t = TRUE, cuof = 1,
    return_plot = TRUE, v = !TRUE
  )$scatter
  fname <- paste0("sf3_coexpression_", paste0(unlist(dgenes[i, ]), collapse = "_"), suffy, ".pdf")
  pdf(fname)
  print(p + viridis::scale_color_viridis(option = "magma"))
  dev.off()
  pdf(gsub(".pdf", "blank.pdf", fname))
  print(shut_it(p, lays = "nsity|ext") + viridis::scale_color_viridis(option = "magma"))
  dev.off()
}

#### ST.3 QC per library ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metadata = read.csv("/home/ciro/covid19/info/libraries_annotation_2020_05_22.csv", stringsAsFactors = FALSE, row.names = 1)
qcfiles = c(
  gex = "/home/ciro/large/covid19/raw/NV030/COUNTS/a1_SingleCell5PE_libraries_summary.csv",
  ab_captute = "/home/ciro/large/covid19/raw/NV030/COUNTS/a1_SingleCell3v2or5_libraries_summary.csv",
  tcr = "/home/ciro/large/covid19/raw/NV030/COUNTS/a1_SingleCellVDJ_libraries_summary.csv"
)
summtabs <- lapply(names(qcfiles), function(x){
  y <- read.csv(qcfiles[x], stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
  cat(paste0(rownames(y)[!rownames(y) %in% rownames(metadata)], collapse = ","), "not found\n")
  z <- cbind(Library = rownames(y), metadata[rownames(y), ], y)
  write.csv(z, file = paste0("sup_table_qc_", x, "_libraries.csv"), row.names = FALSE)
  return(z)
})
headmat(summtabs[[2]])

### 3.X Revising polifunctionality ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "/home/ciro/large/covid19/results/clustering/CD4T6_R1n2_sng2_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData"
genes = c("IFNG", "IL2", "TNF")

## Reading
mycells <- theObjectSavedIn(clustname)
annot <- mycells@meta.data
expdata <- expm1(mycells@assays$RNA@data) * 100; gc()
expdata <- cts2cpm(expdata[, getsubset(c("RNA_snn_res.0.6", "1", "10"), annot, v = TRUE)])
annot = annot[colnames(expdata), ]

## Operations
void <- add_gene_tag(
  lgenes = genes,
  annot = annot,
  mat = expdata,
  thresh = 1,
  v = TRUE
)
annot <- joindf(annot, void)
dgenes <- rep(genes, 2)
dir.create('polyfunctionality')
mycells$Cluster_virus <- paste0(mycells$RNA_snn_res.0.6, "_", mycells$orig.virus2)
ddf <- FetchData(
  object = mycells,
  vars = c("Cluster_virus", genes),
  cells = getsubset(c("orig.virus2", "CV", "FLU"), annot, v = TRUE)
)
p <- violins(
  dat = ddf,
  xax = "Cluster_virus",
  yax = genes,
  dots = TRUE,
  combine = FALSE,
  extra_theme = RotatedAxis()
)
pdf(paste0("polyfunctionality/violins.pdf"), height = 14)
print(plot_grid(plotlist = p[-length(p)], ncol = 1))
dev.off()
for(cl in c("1", "10")[1]){
  for(virus in c("", "CV", "FLU")[1]){
    firstfilt <- list(c("RNA_snn_res.0.6", cl), c("orig.virus2", virus))[ifelse(virus == "", -2, 1:2)]
    scells <- getsubset(firstfilt, annot, v = TRUE)
    fname <- paste0("polyfunctionality/in_cluster", cl, virus)
    genes_props <- reshape2::melt(table(void[scells, ]))
    write.csv(genes_props, file = paste0(fname, ".csv"))

    for(i in 1:length(genes)[1]){
      newfilt <- c(firstfilt, gettag(dgenes[i]))
      scells1 <- getsubset(newfilt, annot, v = TRUE)
      cat("Using", length(scells1), "cells\n")
      p <- get_densities(
        mat = expdata[, scells1],
        genes = dgenes[(i+1):(i+2)],
        log2t = TRUE, cuof = 1,
        return_plot = TRUE, v = !TRUE
      )$scatter
      pdf(paste0(fname, "_", paste0(dgenes[i:(i+2)], collapse = "_"), ".pdf"))
      print(p + viridis::scale_color_viridis(option = "magma"))
      dev.off()
      pdf(paste0(fname, "_", paste0(dgenes[i:(i+2)], collapse = "_"), "_blank.pdf"))
      print(shut_it(p, lays = "text|dens") + viridis::scale_color_viridis(option = "magma"))
      dev.off()
    }
  }
}
