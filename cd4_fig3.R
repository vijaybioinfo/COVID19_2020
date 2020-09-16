#!/usr/bin/R

############
# Figure 3 #
############

# This script will create plots for figure 3 from cd4 data

source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
library(Seurat)
library(cowplot)
theme_set(theme_cowplot())

dirfig <- '/home/ciro/large/covid19/results/a1_final_figures_cd4/Figure_3'
setwdc(dirfig)
clustnamebk <- "none"

# Global variables
colsname <- "../data/global_colours.csv"
gr.cols <- readfile(colsname, stringsAsFactors = FALSE, check.name = FALSE, row.names = 1)

#### 3.X Proportions ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
nres = 'RNA_snn_res.0.6'
prefix = c('f3a_', 'f3a_')

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
scells <- getsubset(list(c("orig.virus2", "CV"), c('orig.hospital', 'Yes', 'No')), mycells@meta.data, v = TRUE)
dfplot <- FetchData(
  object = mycells,
  cells = scells,
  vars = c(colnames(mycells@meta.data), c('UMAP_1', 'UMAP_2'))
)
dfplot$Cluster <- factormix(dfplot[, nres])
dfplot$hospital <- ifelse(dfplot$orig.hospital == "Yes", "Hospital", "Non-hospital")
table(dfplot$hospital)

# ---- Start generating the figures
scells <- sample_grp(annot = dfplot, cname = 'hospital', v = TRUE)
p <- ggplot(dfplot[scells, ], aes(x = UMAP_1, y =  UMAP_2, color = Cluster)) +
  geom_point() +
  facet_wrap(facets = ~hospital) +
  labs(color = NULL) +
  mytheme +
  SetLegendPointsGG()

pdf(paste0(prefix[1], "umap_hospital.pdf"), width = 15, height = 9.5)
print(p)
dev.off()
pdf(paste0(prefix[1], "umap_hospital_blank.pdf"), width = 15, height = 9.5)
print(shut_it(p))
dev.off()

pp <- plot_pct(x = dfplot, groups = c("hospital", "Cluster"), orderby = "Hospital", normalise = TRUE, return_table = TRUE)
propdf <- pp$table
write.csv(propdf, file = paste0(prefix[1], "proportions_hospital.csv"))

p <- pp$plot + labs(fill = NULL) + coord_flip() + mytheme

pdf(paste0(prefix[2], "bar_hospital.pdf"), width = 5, height = 10)
print(p)
dev.off()
pdf(paste0(prefix[2], "bar_hospital_blank.pdf"), width = 4, height = 10)
print(shut_it(p))
dev.off()

tvar <- dfplot[getsubset(list(c("orig.virus2", "CV"), c('orig.hospital', 'Yes')), dfplot, v = TRUE), ]
table(tvar[, c("orig.sex", "Cluster")], useNA = 'always')
p <- plot_pct(
  x = tvar,
  groups = c("orig.sex", "Cluster"), orderby = "Male", normalise = TRUE
) + labs(fill = NULL) + coord_flip() + mytheme
pdf(paste0(prefix[2], "bar_sex.pdf"), width = 5, height = 10)
print(p)
dev.off()

clustnamebk <- clustname

#### S3.X Volcano ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
prefix = "sf3b"
padjthr = 0.05
fcthr = .25
showgenes = c('PRF1', 'GZMB', 'ZBTB32', 'ZBED2', 'HAVCR2', 'LAG3', 'PDCD1', 'TIGIT', 'CD70', 'DOK5', 'DUSP4', 'add')
resname = "../data/cv_5vs0.csv"
trimmer = "c350"

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
if(is.null(showgenes) || any(showgenes %in% "add")) showgenes <- c(showgenes, bordering(modres[tvar, ], cnames = "log2FoldChange", n = 10))
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/volcano.R')
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


### 3.X Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "../data/CD4T6_seurat.rdata"
allconfigs <- list(
  vconfig = list(
    name = "sf3c_violins_tfh",
    genes = c('TIGIT', 'LAG3', 'HAVCR2', 'PDCD1', 'DUSP4', 'CD70', 'DOK5'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '5', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 4, size = c(3, 21),
    order = c('5', '0', '7')
  ),
  vconfig = list(
    name = "sf3d_violins_tfh",
    genes = c('PRF1', 'GZMB'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '5', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 1, size = c(8, 4),
    order = c('5', '0', '7')
  ),
  vconfig = list(
    name = "sf3f_ctl_violins",
    genes = c('CCL3', 'CCL4', 'CCL5', 'XCL1', 'XCL2'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 1, size = c(16, 4),
    order = c('6', '9', 'REST')
  ),
  vconfig = list(
    name = "f3b_tfs_violins",
    genes = c('ZBTB32', 'ZBED2'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '5', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(8, 4),
    order = c('5', '0', '7')
  ),
  vconfig = list(
    name = "f3d_violins",
    genes = c('HOPX', 'ZEB2', 'CD72', 'GPR18', 'SLAMF7'),
    sfilter = list(c('RNA_snn_res.0.6', 0:12), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 2, size = c(8, 4),
    order = c('6', '9', 'REST')
  ),
  vconfig = list(
    name = "nf3_violins",
    genes = c('PRDM1'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '5', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    ncols = 1, size = c(8, 4),
    order = c('0', '5', '7')
  )
)

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

# Run fig2.R - 2.E Violin plots

clustnamebk <- clustname

### 3.X Curtain plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
allconfigs <- list(
  vconfig = list(
    name = "f3b_curtain",
    genes = c('TIGIT', 'LAG3', 'HAVCR2', 'PDCD1', 'DUSP4', 'CD70', 'DOK5'),
    sfilter = list(c('RNA_snn_res.0.6', '0', '5', '7'), c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    order = c('5', '0', '7')
  ),
  vconfig = list(
    name = "sf3e_curtain",
    genes = c('HOPX', 'ZEB2', 'CD72', 'GPR18', 'SLAMF7'),
    sfilter = list(c('orig.virus2', 'CV')),
    ident = c('RNA_snn_res.0.6'),
    order = c('6', '9', 'REST')
  )
)

## Reading
if((!exists("mycells")) || clustname != clustnamebk) xtheObjectSavedIn(clustname)

## Operations
for(vconfig in allconfigs){
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
}

clustnamebk <- clustname

### 3.X Co-expression/double positive proportions ###%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
fannot = "../data/CD4T6_metadata.rdata"

prefix = "f3c_scatter_ctox_"
genes = c("PRF1", "GZMB")
selectss = list(list(c("RNA_snn_res.0.6", "5")), list(c("RNA_snn_res.0.6", "0")), list(c("RNA_snn_res.0.6", "7")))

prefix = "sf3_scatter_ctox_"
genes = c("XCL1", "XCL2")
selectss = list(list(c("RNA_snn_res.0.6", "6")), list(c("RNA_snn_res.0.6", "9")), list(c("RNA_snn_res.0.6", "11")))

prefix = "../Figure_2/sf2_scatter_polyfunc_"
genes = c("TNF", "IL2")
selectss = list(c("RNA_snn_res.0.6", "1"), c("tag_IFNG", "IFNG+"))

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- readfile(clustname)
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

if(any(grepl("tag_", unlist(selectss)))){
  mygenes <- unlist(selectss)[grepl("tag_", unlist(selectss))]
  mygenes <- gsub("tag_|_[0-9]{1,}", "", mygenes)
  mygenes <- getfound(mygenes, rownames(expdata))
  if(length(mygenes) > 0){
    void <- add_gene_tag(
      lgenes = mygenes,
      annot = annot,
      mat = expdata,
      thresh = 1,
      v = TRUE
    )
    annot <- joindf(annot, void)
  }
}

## Operations
dgenes <- remove.factors(create_pairs(data.frame(g = genes, stringsAsFactors = FALSE))[, 2:1])
for(i in 1:nrow(dgenes)){
  cat(commas(dgenes[i, ]), "\n")
  if(!is.list(unlist(selectss, recursive = FALSE))) selectss <- list(selectss)
  for(sselectss in selectss){
    scells <- getsubset(sselectss, annot, v = TRUE)
    suffy <- make.names(paste0(datatype, "_", summary_subset(sselectss)))
    p <- get_densities(
      mat = expdata[unlist(dgenes[i, ]), scells],
      genes = unlist(dgenes[i, ]), usedp = TRUE,
      log2t = TRUE, cuof = 1,
      return_plot = TRUE, v = !TRUE
    )$scatter
    fname <- paste0(prefix, paste0(unlist(dgenes[i, ]), collapse = "_"), "_", suffy, ".pdf")
    pdf(fname)
    print(p + viridis::scale_color_viridis(option = "magma"))
    dev.off()
    pdf(gsub(".pdf", "blank.pdf", fname))
    print(shut_it(p, lays = "nsity|ext") + viridis::scale_color_viridis(option = "magma"))
    dev.off()
  }
}

clustnamebk <- clustname

#### ST.5 ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(openxlsx)
## Input
colnamestype <- c(gene_name = "TEXT", hange = "NUMBER", pvalue = "SCIENTIFIC", padj = "SCIENTIFIC", "_mean" = "NUMBER", "_exprFra" = "NUMBER")
fnames = c(
  "5 vs 0" = "../data/cv_5vs0.csv",
  "5 vs 7"= "../data/cv_5vs7.csv",
  "Cluster 5 Severe vs Mild" = "../data/cluster5_YesvsNo.csv",
  "6 vs 9" = "../data/cv_6vs9.csv"
)
title_name <- "Differential expression analysis"
prefix = "supp_table6"
select_features <- function(x) getDEGenes(x = x, pv = 0.05, fc = 0.25, v = !TRUE)

## Reading
mytables <- lapply(fnames, read.csv, row.names = 1, stringsAsFactors = FALSE)
if(is.null(mytables)) names(mytables) <- basename(dirname(fnames))

## Operations
workfile <- createWorkbook(
  creator = "Ciro",
  title = title_name
)

## create and add a style to the column headers
body_style <- createStyle(
  fontSize = 12, fontColour = "#000000",
  numFmt = "GENERAL",
  border = "TopBottomLeftRight", borderColour = "#000000", borderStyle = "thin",
  fgFill = "#FFFFFF",
  halign = "center", valign = "center",
)
feature_style <- createStyle(
  fontSize = 12, fontColour = "#000000",
  numFmt = "TEXT",
  border = "LeftRight", borderColour = "#000000", borderStyle = "thick",
  fgFill = "#D9D9D9",
  halign = "center", valign = "center",
  textDecoration = "italic"
)
create_header_style <- function(f = "GENERAL"){
  createStyle(
    fontSize = 12, fontColour = "black",
    numFmt = f,
    border = "TopBottomLeftRight", borderColour = "#000000", borderStyle = "thick",
    fgFill = "#BFBFBF",
    halign = "center", valign = "center",
    textDecoration = "bold"
  )
}

for(tname in names(mytables)){
  cat(tname, "\n")
  addWorksheet(wb = workfile, sheetName = tname)
  ddf <- mytables[[tname]]
  ddf <- ddf[select_features(x = ddf), ]
  ddf$gene_name <- gsub("'", "", ddf$gene_name)
  taken_cols <- found_partial(names(colnamestype), colnames(ddf))
  cat("Columns:", commas(taken_cols, Inf), "\n")
  ddf <- ddf[, taken_cols]
  for(colname in colnames(ddf)){
    myformat <- colnamestype[[which(sapply(names(colnamestype), grepl, colname))]]
    if(myformat == 'NUMBER') ddf[, colname] <- round(ddf[, colname], 2)
  }
  str(ddf)
  writeData(wb = workfile, sheet = tname, x = paste0(title_name, ". ", tname))
  writeDataTable(
    wb = workfile,
    sheet = tname,
    x = ddf,
    startRow = 4
  )

  addStyle(wb = workfile, sheet = tname, style = body_style, rows = 4:(nrow(ddf)+4), cols = 1:ncol(ddf), gridExpand = TRUE, stack = FALSE)
  addStyle(wb = workfile, sheet = tname, style = feature_style, rows = 4:(nrow(ddf)+4), cols = 1, gridExpand = TRUE, stack = FALSE)
  for(colname in colnames(ddf)){
    myformat <- colnamestype[[which(sapply(names(colnamestype), grepl, colname))]]
    addStyle(wb = workfile, sheet = tname,
      style = create_header_style(f = myformat),
      rows = 4, cols = which(colnames(ddf) %in% colname), gridExpand = TRUE, stack = FALSE
    )
  }
  setColWidths(wb = workfile, sheet = tname, cols = 1:ncol(ddf), widths = 16)
  setColWidths(wb = workfile, sheet = tname, cols = 1:ncol(ddf), widths = 16)
}
saveWorkbook(wb = workfile, file = paste0(prefix, ".xlsx"), overwrite = TRUE)

#### ST.3 QC per library ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metadata = read.csv("../data/libraries_annotation.csv", stringsAsFactors = FALSE, row.names = 1)
qcfiles = "/home/ciro/large/covid19/raw/NV031/COUNTS" # 'cellranger counts' output
qc_groups <- make_list(x = metadata, colname = "Library_type")
cellranger_order <- list(
  Gex = list(
    per_cell = c(
      'Estimated Number of Cells', 'Mean Reads per Cell', 'Median Genes per Cell',
      'Number of Reads', 'Valid Barcodes', 'Fraction Reads in Cells',
      'Total Genes Detected', 'Median UMI Counts per Cell'
    ),
    mapping = c(
      'Reads Mapped to Genome', 'Reads Mapped Confidently to Genome', 'Reads Mapped Confidently to Intergenic Regions',
      'Reads Mapped Confidently to Intronic Regions', 'Reads Mapped Confidently to Exonic Regions',
      'Reads Mapped Confidently to Transcriptome', 'Reads Mapped Antisense to Gene'
    ),
    seq = c(
      'Sequencing Saturation', 'Q30 Bases in Barcode', 'Q30 Bases in RNA Read',
      'Q30 Bases in RNA Read 2', 'Q30 Bases in Sample Index', 'Q30 Bases in UMI'
    )
  ),
  TCR = list(
    per_cell = c(
      'Estimated Number of Cells', 'Mean Read Pairs per Cell', 'Number of Cells With Productive V-J Spanning Pair',
      'Number of Read Pairs', 'Valid Barcodes', 'Fraction Reads in Cells', 'Median TRA UMIs per Cell', 'Median TRB UMIs per Cell'
    ),
    mapping = c(
      'Reads Mapped to Any V(D)J Gene', 'Reads Mapped to TRA', 'Reads Mapped to TRB', 'Mean Used Read Pairs per Cell'
    ),
    seq = c(
      'Q30 Bases in Barcode', 'Q30 Bases in RNA Read 1', 'Q30 Bases in RNA Read 2', 'Q30 Bases in Sample Index', 'Q30 Bases in UMI'
    ),
    tcr_info = c(
      'Cells With Productive V-J Spanning Pair', 'Cells With Productive V-J Spanning (TRA, TRB) Pair',
      'Paired Clonotype Diversity', 'Cells With TRA Contig', 'Cells With TRB Contig',
      'Cells With CDR3-annotated TRA Contig', 'Cells With CDR3-annotated TRB Contig',
      'Cells With V-J Spanning TRA Contig', 'Cells With V-J Spanning TRB Contig',
      'Cells With Productive TRA Contig', 'Cells With Productive TRB Contig'
    )
  ),
  CITE = list(
    per_cell = c(
      'Estimated Number of Cells', 'Antibody: Mean Reads per Cell', 'Antibody: Valid Barcodes'
    ),
    ab_info = c(
      'Antibody: Fraction Antibody Reads', 'Antibody: Fraction Antibody Reads Usable', 'Antibody: Antibody Reads Usable per Cell',
      'Antibody: Fraction Reads in Barcodes with High UMI Counts', 'Antibody: Fraction Unrecognized Antibody',
      'Antibody: Antibody Reads in Cells', 'Antibody: Median UMIs per Cell (summed over all recognized antibody barcodes)'
    ),
    seq = c(
      'Antibody: Sequencing Saturation', 'Antibody: Q30 Bases in Barcode', 'Antibody: Q30 Bases in Antibody Read',
      'Antibody: Q30 Bases in Sample Index', 'Antibody: Q30 Bases in UMI'
    )
  )
)
summtabs <- lapply(names(qc_groups), function(x){
  cat("Group:", x, "\n")
  fnames <- paste0(dircheck(qcfiles), qc_groups[[x]], "/outs/metrics_summary.csv")
  names(fnames) <- qc_groups[[x]]
  tvar <- !file.exists(fnames)
  if(any(tvar)) cat(paste0(names(fnames)[tvar], collapse = ","), "not found\n")
  fnames <- fnames[file.exists(fnames)]
  y <- t(sapply(fnames, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  y <- data.frame(y, stringsAsFactors = FALSE, check.names = FALSE)
  tvar <- !rownames(y) %in% rownames(metadata)
  if(any(tvar)) cat(paste0(rownames(y)[tvar], collapse = ","), "not found\n")
  cnames <- unlist(cellranger_order[[x]], use.names = FALSE); tvar <- cnames %in% colnames(y)
  if(any(!tvar)) cat(paste0(cnames[!tvar], collapse = "\n"), "\n");
  cnames <- cnames[tvar]
  z <- cbind(Library = rownames(y), metadata[rownames(y), ], y[, cnames])
  # write.csv(z, file = paste0("sup_table_qc_", x, "_libraries.csv"), row.names = FALSE)
  data.table::fwrite(z, file = paste0("sup_table_qc_", x, "_libraries.csv"), row.names = FALSE)
  return(z)
})
headmat(summtabs[[2]])

#### 3.E Gene expression ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "../data/CD4T6_seurat.rdata"
genes = c("CCL3", "CCL4", "CCL5", "XCL1", "XCL2")
nres = "RNA_snn_res.0.6"
prefix = "f3e_scatter_expression_"

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
ddf <- FetchData(mycells, vars = c(nres, genes))
for(g in genes){
  cat(g, "\n")
  p <- FeaturePlot(
    object = mycells,
    features = g,
    cols = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
    # cols = c('#fffeee', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
    # cols = c("#fffeee", "#ffe080", "#ffc100", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000")
  )
  pdf(paste0(prefix, g, ".pdf"))
  print(p)
  dev.off()
  pdf(paste0(prefix, g, "_blank.pdf"))
  print(p)
  dev.off()
}
clustnamebk = clustname

#### X.X 0n6 UMAP ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "data/CD4T0n6_seurat.rdata"
nres = "orig.stim_time"

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
mycells@meta.data$Cluster <- factor(mycells@meta.data[, nres], c("6", "0"))
p <- DimPlot(
  object = mycells,
  reduction = 'umap',
  group.by = "Cluster",
  label = TRUE
) + scale_color_brewer(palette = "Set1")
pdf('nf3_0n6_umap.pdf')
print(p)
dev.off()
pdf('nf3_0n6_umap_blank.pdf')
print(shut_it(p, lays = "ext"))
dev.off()
clustnamebk = clustname

### X.X Signature calculation on 0n6 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/clustering_utilities.R')
## Input
clustname = "data/CD4T0n6_seurat.rdata"
cpatterns = c("RNA_snn_res.0.6", "orig.hospital", "RNA_snn_res.hospital", "orig.stim_time")
ssigna_patterns = "tfh|tiva|patil"
out_dir = "./CD4T0n6_R1n2n3_sng_25p_"
pcs.comp = 30
norm_type = "RNA_snn_res"; verb = TRUE
rdims = list(umap = c('UMAP_1', 'UMAP_2'))

# Run fig2.R - X.X Signature calculation

#### X.X 0n6 Signatures ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
clustname = "data/CD4T0n6_seurat.rdata"
nres = 'orig.stim_time'
signaturef = "./CD4T0n6_R1n2n3_sng_25p_signatures_30PCs/signatures.csv"
ssigna_patterns = "tfh|tiva|patil"
prefix = c("nf3_0n6_", "nf3_0n6_")
couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")

# Run fig2.R - 2.F Signatures

#### X.X 0n6 UMAP expression ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
myset = "nf3_umapCD4T0n6"
clustname = "data/CD4T0n6_seurat.rdata"
genes <- c('PRF1', 'GZMB', 'NKG7', 'GNLY')
# CD4T6
myset = "f3_umapCD4T6"
clustname = "../data/CD4T6_seurat.rdata"
genes <- c('PRF1', 'GZMB', 'CCL3', 'CCL4', 'CCL5', 'XCL1', 'XCL2', 'NKG7', 'GNLY')

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)
genes <- getfound(genes, rownames(mycells), v = TRUE)

## Operations
dname <- paste0(myset, "/")
dir.create(dname)
for(gg in genes){
  cat(gg, "\n")
  fname <- paste0(dname, gg, ".pdf")
  p <- FeaturePlot(
    object = mycells,
    reduction = 'umap',
    features = gg
  ) + scale_colour_gradientn(colours = couls)
  pdf(fname)
  print(p)
  dev.off()
}

clustnamebk <- clustname

### X.X 0n6 Violin plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(grid)
## Input
clustname = "data/CD4T0n6_seurat.rdata"
allconfigs <- list(
  vconfig = list(
    name = "nf3_0n6_cytokines",
    genes = c('IFNG', 'IL2', 'TNF', 'CSF2', 'IL3', 'CXCL13', 'IL21'),
    sfilter = list(c('orig.stim_time', '0', '6')),
    ident = c('orig.stim_time'),
    ncols = 2, size = c(12, 8),
    order = c('0', '6')
  ),
  vconfig = list(
    name = "nf3_0n6_activation",
    genes = c('MIR155HG', 'TNFRSF4', 'TNFRSF18', 'TNFRSF9'),
    sfilter = list(c('orig.stim_time', '0', '6')),
    ident = c('orig.stim_time'),
    ncols = 2, size = c(12, 8),
    order = c('0', '6')
  )
)

# run fig2.R - 2.E Violin plots

#### 3.X Cluster proportions per donor ####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Input
# CD4T6
clustname = "../data/CD4T6_seurat.rdata"
nres = 'RNA_snn_res.0.6'
dcolumn = 'orig.donor'
subdiv = c("orig.hospital", "orig.sex")
sselect <- list(c(subdiv, "Yes", "No"), c("orig.virus2", "CV"))
corder = c('0', '7', '5', '6', '2', '3')
suffy = NULL
roder = NULL
add_props = NULL

## Reading
if((!exists("mycells")) || clustname != clustnamebk) mycells <- theObjectSavedIn(clustname)

## Operations
scells <- getsubset(sselect, mycells@meta.data, v = TRUE)
annot <- remove.factors(mycells@meta.data[scells, ])
table(annot[, c(dcolumn, subdiv)])
dgroup <- sapply(subdiv, function(subd){
  tvar <- make_list(annot, colname = dcolumn, col_objects = subd)
  tvar <- reshape2::melt(lapply(tvar, function(x) paste0(unique(x), collapse = "-") ))
  tmp <- as.character(tvar[, 1])
  names(tmp) <- as.character(tvar[, 2])
  tmp
})

ddfprops <- table(annot[, c(dcolumn, nres)])
# Proportions heatmap
matpct <- as.matrix(as.data.frame.matrix(prop.table(ddfprops, margin = 1)))
if(!is.null(add_props)){
  tscells <- scells
  if(!is.null(add_props$newobject)){
    add_cells <- theObjectSavedIn(add_props$newobject)
    add_annot <- add_cells@meta.data
  }else{
    add_annot <- annot
  }
  scells2 <- getsubset(add_props$filter, add_annot, v = TRUE)
  if(any(grep("tag_", add_props$vars))){
    tmp <- add_props$vars[grep("tag_", add_props$vars)]
    void <- add_gene_tag(
      lgenes = gsub("tag_", "", tmp),
      annot = add_annot,
      mat = if(!is.null(add_props$newobject)) add_cells@assays$RNA@data else mycells@assays$RNA@data,
      thresh = 0,
      v = TRUE
    )
    add_annot <- joindf(void, add_annot)
  }
  matpctnew <- table(add_annot[scells2, c(add_props$vars)])
  if(is.null(add_props$totals)){
    matpctnew <- as.matrix(as.data.frame.matrix(prop.table(matpctnew, margin = 1)))
  }else{
    matpctnew <- apply(matpctnew, 2, function(x) x / table(add_annot[tscells, c(add_props$totals)])[names(x)] )
  }
  matpct <- joindf(data.frame(matpct, check.names = FALSE), data.frame(matpctnew, check.names = FALSE))
}
matpct <- round(matpct, 2)
matpct <- if(!is.null(roder)){
  suffix = paste0(suffy, "_supervised")
  matpct[roder, corder]
}else if(!is.null(corder)){
  suffix = paste0(suffy, "_unsupervised")
  matpct[, corder]
}else{ matpct }
topz <- c(0, .7)
matpct[matpct > topz[2]] <- topz[2]; matpct[matpct < topz[1]] <- topz[1];
palettebreaks <- seq(from = 5, to = 70, by = 5)
annoc <- data.frame(dgroup); colnames(annoc) <- casefold(sub("orig.", "", colnames(annoc)), upper = FALSE)
annoc <- annoc[rownames(matpct), sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
anncolist <- lapply(annoc, function(x) v2cols(x, gr.cols, v = TRUE) )
mypalette <- colorRampPalette(colors = c("white", "red"), space = 'Lab')
pdf(paste0("donor_heatmap", suffix, ".pdf"), width = 10, height = 12, onefile = FALSE)
x <- NMF::aheatmap(
  x = matpct, annRow = annoc, annColors = anncolist,
  scale = 'none', Rowv = ifelse(is.character(roder), NA, TRUE), Colv = NA,
  col = mypalette(length(palettebreaks) - 1)
)
graphics.off()

clustnamebk <- clustname
