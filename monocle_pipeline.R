#!/share/apps/R/3.5/bin/Rscript
##############
# MONOCLE 3 ##
##############

## Monocle 3 has been re-engineered to analyze large, complex single-cell datasets.
# The algorithms at the core of Monocle 3 are highly scalable and can handle millions of cells.
# Monocle 3 will add some powerful new features that enable the analysis of organism- or embryo-scale experiments:
# A better structured workflow to learn developmental trajectories.
# Support for the UMAP algorithm to initialize trajectory inference.
# Support for trajectories with multiple roots.
# Ways to learn trajectories that have loops or points of convergence.
# Algorithms that automatically partition cells to learn disjoint or parallel trajectories using ideas from "approximate graph abstraction".
# A new statistical test for genes that have trajectory-dependent expression.
# A 3D interface to visualize trajectories and gene expression.

#### Check arguments ####
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
.libPaths(new = "/mnt/BioHome/ciro/R/newer_packs_library/3.5/")
defargs <- list(
  mycellsf = '',
  rawdata = '',
  outdir = '/mnt/BioScratch/cramirez/',
  catg = 'res.0.6',
  catgfilt = 'mycolumn~myclass',
  smpcells = FALSE,
  cpmconv = FALSE,
  gcolour = 'no_file',
  topgenes = 16,
  myseed = 27,
  # specific paramaters #
  redtype = 'tSNE', #c("DDRTree", "ICA", "tSNE", "UMAP")
  seudisps = TRUE,
  pcs_comp = 50,
  resMF = NULL,
  use_obj = TRUE,
  genesetf = "no_file",
  ctstype = 'raw',
  filtcells = FALSE
)
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')

flaggy <- vector()
for (varb in names(defargs)) {
  cat(paste0(varb, ': '))
  if(varb %in% names(args)){
    checkmate <- args[[varb]] == ifelse(is.null(defargs[[varb]]), 'NULL', defargs[[varb]])
    flaggy <- c(flaggy, checkmate); aclass <- class(defargs[[varb]])
    aclass <- ifelse(aclass == 'NULL', 'character', aclass); special <- grepl("\\(|,", args[[varb]])
    if(special){
      assign(varb, eval(expr = parse(text = args[[varb]])))
    }else{ assign(varb, methods::as(args[[varb]], aclass)) }
  }else{ assign(varb, defargs[[varb]]); flaggy <- c(flaggy, FALSE) }
  print(eval(as.name(varb)))
  args[[varb]] <- eval(as.name(varb))
}; names(flaggy) <- names(defargs)
set.seed(myseed)

# extra parameters
rnloc = 1
topgenes = unique(c("FOXP3", "TNFRSF18", "CCR8", "IL1R2", "CXCL13", "BCL6", "CXCR5", "TCF7", "TOP2A", "TNFRSF9", "ICOS", "CD28", "TIGIT", "CTLA4",
  'IL2RA', 'FOXP3', 'CD4', 'CTLA4', 'TNFRSF18', 'TNFRSF9', 'PDCD1', 'MKI67',
  "Hmgb2", "H2afz", "Birc5", "Cdc20", "Tpx2", "Ccr9", "Cd2", "Bcl2", "Gzmb", "Nkg7", "Klra9", "Ccl5", "Cxcr3"))

outdir <- dircheck(outdir)

# specific modifications
if(!dir.exists(paste0(outdir, '../..'))){
  cat('Parent doesn\'t exist: output at /mnt/BioScratch/cramirez/\n')
  outdir <- '/mnt/BioScratch/cramirez/'
}
setwdc(outdir)

### Log file ####
tvar <- 'mono'
tvar <- ifelse(isTRUE(seudisps[1]), 'seu', tvar)
tvar <- ifelse(file.exists(as.character(genesetf[1])), 'gs', tvar)
if((tvar %in% c('seu', 'gs')) && isTRUE(use_obj)) tvar <- paste0(tvar, "obj")
fsufix <- paste0(c(redtype, tvar, paste0(pcs_comp, 'PC'), ctstype), collapse = '_')
dir.create('logs'); hashy <- digest::digest(paste0(args, collapse = ""), "md5", serialize = FALSE)
dir.create('data/')
log.file <- paste0('logs/monocle3_', hashy, '.log')
if(!file.exists(log.file)) file.create(log.file)
out.file <- file(log.file, open = 'wt')
sink(out.file) ; sink(out.file, type = 'message')
cat('Date and time:\n'); st.time <- timestamp()
#####

#### Spcific functions #### ----------------------------------------------------
# a helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type = NULL){
  catgcells <- as.character(pData(cds)[, cell_phenotype])
  if(is.null(root_type)){
    root_type <- names(head(sort(table(catgcells)), 1))
  }
  cell_ids <- which(as.character(catgcells) == root_type)

  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
#### ####### ######### #### ----------------------------------------------------

#### Loading dependencies #### -------------------------------------------------
depend <- c('Seurat', 'monocle', 'VGAM', 'dplyr', 'reticulate', 'gridExtra', 'cowplot')
load_packs(depend)
#### Presenting parameteres #### -----------------------------------------------
present_params(parms = c(
  'mycellsf' = 'Meta data',
  'rawdata' = 'Count data',
  'outdir' = 'Output directory',
  'catg' = 'Category',
  'catgfilt' = 'Filter',
  'smpcells' = 'Sample cells',
  'cpmconv' = 'CPM conversion',
  'gcolour' = 'Colours file',
  'topgenes' = 'nGenes to plot',
  'myseed' = 'Seed',
  # specific paramaters #
  'redtype' = 'Dimentional reduction',
  'seudisps' = 'Seurat genes',
  'pcs_comp' = 'Components to use',
  'resMF' = 'Regression Fromula',
  'use_obj' = 'Use given object',
  'ctstype' = 'Count type',
  'filtcells' = 'Filter cells'
), flaggy = flaggy, fast_track = TRUE)

### Preparing data #### -------------------------------------------------------
cat('-------- Reading data --------\n')
cat('Getting annotation data - ')
if(file.exists(mycellsf) && !dir.exists(mycellsf)){
  mycells <- readfile(mycellsf, v = TRUE, row.names = rnloc, stringsAsFactors = FALSE, check.names = FALSE)
  if(casefold(class(mycells)) == 'seurat'){ # Dealing with the seurat object
    annot <- mycells@meta.data
  }else{
    annot <- mycells; rm(mycells)
  }
}else{ stop('Annotation file must be seurat object or .csv') }
cat('Meta data:\n')
print(str(annot))

# subsetting annotation
catgfilt <- translist(catgfilt)
if(any(sapply(catgfilt, head, 1) %in% colnames(annot))){
  annot <- annot[getsubset(catgfilt, annot, v = T), ]
  print(sapply(sapply(catgfilt, head, 1), function(x) table(annot[, x[1]])))
}

catg <- unlist(translist(catg))
catg <- getfound(catg, colnames(annot), v = T)
catg <- catg[sapply(catg, function(x) length(table(annot[, x])) > 1 )]
if(isTRUE(smpcells)){
  tvar <- min(table(annot[, catg[1]]))
  if(is.numeric(smpcells)){
    if(smpcells < tvar) tvar <- smpcells; smpcells <- TRUE
  }
  cat('Sampling cells <------------\n')
  set.seed(myseed)
  mynames <- lapply(mixedsort(unique(annot[, catg[1]])), function(x){
    sample(getsubset(c(catg[1], x), annot), tvar)
  })
  annot <- annot[mynames, ]
}

cat('\nReading expression data - ')
if(exists('mycells')){
  cat('from given object ')
}else if(file.exists(rawdata) && !dir.exists(rawdata)){ # Dealing with the seurat object
  cat('Loading from ')
  mycells <- readfile(rawdata, v = TRUE, row.names = rnloc)
}else if(dir.exists(rawdata)){
  cat('from folder ')
  rawdata <- setdir(rawdata, pats = c('outs', 'filtered_gene_bc_matrices_mex', cd, 'hg19'), v = F)
  cat(rawdata, '\n'); mycells <- Seurat::Read10X(data.dir = rawdata)
}else{
  stop('10X data not found. Must be CSV, RData (if seurat object) or a folder.')
}

if(casefold(class(mycells)) == 'seurat'){
  DefaultAssay(mycells) <- "RNA"
  ctstype <- if(ctstype %in% c('tpm', 'fpkm', 'raw', 'counts')) "counts" else ctstype
  cat("slot", ctstype, 'from Seurat object', DefaultAssay(mycells), '\n')
  cts <- GetAssayData(mycells, slot = ctstype)
}else{
  cts <- mycells
}

cts <- cts[, getfound(rownames(annot), colnames(cts), element = 'cells', v = TRUE), drop = FALSE]
print(t(headmat(cts)))
cat('Genes:', nrow(cts), '\nCells:', ncol(cts), '\n\n')

cat('Getting colours:', gcolour,'\n')
if(file.exists(gcolour)){
  allcolours <- readfile(gcolour, row.names = 1, stringsAsFactors = F)
}else{ allcolours <- NULL }
cell_type_color <- lapply(annot[, catg, drop = FALSE], function(x){
  v2cols(select = x, sour = allcolours, fw = 'gg')
})

#### MAIN PROGRAM #### ---------------------------------------------------------
options(warn = -1) # Turn off warning message globally
# Pass TRUE if you want to see progress output on some of Monocle 3's operations
void <- DelayedArray:::set_verbose_block_processing(TRUE)
# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size = 1000e6)

norm_pseu = list( # setting data type distributions
  raw_umi_counts = c(negbinomial.size(), "log", 0.1),
  tpm_fpkm = c(tobit(Lower = 0.1), "log", 0.1),
  scale_data_logtpm_logfpkm = c(gaussianff(), 'none', 0)
)
norm_pseu <- norm_pseu[[ names(norm_pseu)[grepl(ctstype, names(norm_pseu))][1] ]]
print(norm_pseu)

cat('Creating Monocle object\n')
pd <- new("AnnotatedDataFrame", data = annot)
cat("Duplicated features", sum(duplicated(row.names(cts))), "\n")
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(cts), row.names = row.names(cts)))
if(casefold(class(mycells)) == 'seurat' && use_obj && FALSE){ # importCDS cant deal with Seurat object
  mycells = SubsetData(mycells, cells.use = rownames(annot))
  mycells@meta.data <- cbind_repcol(annot, mycells@meta.data)
  cat('Importing from', class(mycells), '\n'); cds <- importCDS(mycells)
  cat('Updating\n'); cds <- updateCDS(cds)
}else{
  if(class(cts) %in% c("data.frame", "matrix")) {
    cat('- sparse matrix conversion\n'); cts <- as(as.matrix(cts), "sparseMatrix")
  }
  cds <- newCellDataSet(cellData = cts, phenoData = pd, featureData = fd,
    lowerDetectionLimit = norm_pseu[[length(norm_pseu)]], expressionFamily = norm_pseu[[1]])
}

# save.image(file = 'all.RData')

cat('Family:', cds@expressionFamily@vfamily, '\n')
if(grepl('tobit', cds@expressionFamily@vfamily) && TRUE){
  cat('Relative to absolute expression conversion\n')
  fname <- paste0('data/rpc.RData')
  if(!file.exists(fname)){
    rpc_matrix <- relative2abs(cds, method = "num_genes")
    save(rpc_matrix, file = fname)
  }else{ load(fname) }
  cds <- newCellDataSet(rpc_matrix, phenoData = pd, featureData = fd,
    lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
}
cat('------------ Step 1: Normalize and pre-process the data ------------\n')
cat(as.character(Sys.time()), '\n')
if(grepl('negbinomial', cds@expressionFamily@vfamily)){
  # BiocManager::valid()
  cat('Estimate size factors and dispersions\n')
  cds <- estimateSizeFactors(cds); cat(as.character(Sys.time()), '\n')
  cds <- estimateDispersions(cds); cat(as.character(Sys.time()), '\n')
}
if(isTRUE(resMF != 'NULL')){
  resMF <- unlist(translist(resMF)); resMF <- resMF[resMF %in% colnames(pData(cds))]
  if(length(resMF)) resMF <- paste0("~", paste(resMF, collapse = "+")) else resMF <- NULL
  cat('Regression:', resMF, '\n')
}else{ resMF <- NULL }

# if((isTRUE(seudisps[1]) || file.exists(as.character(genesetf[1]))) && !isTRUE(use_obj)){
if(is.na(genesetf[2])) genesetf[2] <- Inf
if((isTRUE(seudisps[1]) || file.exists(as.character(genesetf[1]))) && casefold(class(mycells)) != 'seurat'){
  # you also need stats if you only have a list of genes
  cat('Getting HVG from Seurat\n')
  mycells <- CreateSeuratObject(cts, meta.data = annot)
  # if(ctstype != "data" && min(cts[, 1:3]) > 0) mycells <- NormalizeData(object = mycells, verbose = TRUE)
  if(ctstype != "data") mycells <- NormalizeData(object = mycells, verbose = TRUE)
  if(DefaultAssay(mycells) != "RNA") mycells <- SCTransform(object = mycells, verbose = TRUE)
  if(!length(VariableFeatures(object = mycells))) mycells <- FindVariableFeatures(object = mycells, verbose = TRUE)
  if(!is.infinite(as.numeric(genesetf[2]))) VariableFeatures(mycells) <- head(VariableFeatures(mycells), as.numeric(genesetf[2]))
  fname <- paste0('data/', fsufix, '_seurat.RData')
  save(mycells, file = fname)
}
cat('Choosing genes: ')
if(isTRUE(seudisps[1]) && casefold(class(mycells)) == 'seurat'){
  cat("S\n")
  # ordering_genes <- mycells@var.genes
  # disp_table <- mycells@hvg.info
  ordering_genes <- VariableFeatures(mycells)
  disp_table <- HVFInfo(mycells)
  colnames(disp_table) <- c('mean_expression', 'variance', 'gene.dispersion.scaled')
  disp_table$dispersion_empirical <- disp_table$variance
  disp_table$gene_id <- rownames(disp_table)
  disp_table$dispersion_fit <- min(disp_table[ordering_genes, ]$dispersion_empirical)
  used_genes <- 'Seurat'
}else if(!is.null(cds@dispFitInfo[["blind"]])){
  cat("M\n")
  disp_table <- dispersionTable(cds)
  tvar <- ifelse(is.numeric(genesetf[1]), as.numeric(genesetf[1]), 0.1)
  #ordering_genes <- as.character(subset(disp_table, mean_expression >= tvar)$gene_id)
  ordering_genes <- head(as.character(subset(disp_table, mean_expression >= tvar)$gene_id), genesetf[2])
  used_genes <- 'Monocle'
}else{ cat('Warning: no unsupervised gene selection\n') }
if(file.exists(as.character(genesetf[1]))){
  cat("From file:", genesetf, "\n")
  mytab <- readfile(genesetf, stringsAsFactors = FALSE, row.names = 1, v = TRUE); mytab$X123 <- rownames(mytab)
  tvar <- which(sapply(mytab, function(x) any(x %in% rownames(cds)) ))[1]
  print(head(mytab, 3)); tvar <- mytab[, tvar]
  if("gene_name" %in% colnames(mytab)) tvar <- gsub("'", "", mytab$gene_name)
  rownames(mytab) <- tvar
  print(head(mytab, 20))
  ordering_genes <- getfound(rownames(mytab), rownames(cds))
  used_genes <- 'Selected'
}
if(!exists('disp_table')){
  cat("Building dispersion table\n")
  disp_table <- data.frame(matrix(rep(1:nrow(cds), 4), ncol = 4, byrow = F))
  rownames(disp_table) <- disp_table$gene_id <- rownames(cds)
  colnames(disp_table) <- c('gene_id', 'mean_expression', 'dispersion_fit', 'dispersion_empirical')
}
# if(length(topgenes)){ ordering_genes <- topgenes; used_genes <- 'Selected genes' }
cat(length(ordering_genes), 'from', casefold(used_genes), '\n')

pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))
distr <- log10(pData(cds)$Total_mRNAs)
upper_bound <- 10^(mean(distr) + 2*sd(distr)); lower_bound <- 10^(mean(distr) - 2*sd(distr))
tvar <- pData(cds)$Total_mRNAs > lower_bound & pData(cds)$Total_mRNAs < upper_bound | (!filtcells)
cat('Using', sum(tvar), 'of', length(tvar), 'cells\n')
print(cbind(
  Total = table(annot[rownames(pData(cds)), catg[1]]),
  Filtered = table(annot[rownames(pData(cds))[tvar], catg[1]])
))

fname <- paste0('filt', filtcells, '_', cds@expressionFamily@vfamily, ctstype, '_',
  used_genes, length(ordering_genes), '_', sum(tvar), '.pdf')
if(!file.exists(fname) && FALSE){
  cat("Fitlering plots\n")
  void <- lapply(catg, function(catgy){
    ggplot(pData(cds), aes_string("Total_mRNAs", color = catgy)) + geom_density() +
      geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)
  })
  g <- list(ggplot(pData(cds), aes_string("Total_mRNAs")) + geom_density() +
    geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound))
  cds <- cds[, tvar]
  g[[length(g)+1]] <- ggplot(pData(cds), aes_string("Total_mRNAs")) + geom_density() +
    geom_vline(xintercept = lower_bound) + geom_vline(xintercept = upper_bound)
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(log(exprs(cds))))))
  g[[length(g)+1]] <- qplot(value, geom = "density", data = melted_dens_df) +
    stat_function(fun = dnorm, size = 0.5, color = 'red') +
    xlab(paste0("Standardized ", "log(", ctstype, ")")) + ylab("Density")
  # g[[length(g)+1]] <- plot_pc_variance_explained(cds, return_all = FALSE,
  #     use_existing_pc_variance = TRUE, norm_method = norm_pseu[[2]])
  g[[length(g)+1]] <- qplot(mean_expression, dispersion_empirical, data = disp_table,
    log = ifelse(grepl("Seu", used_genes), "", "xy"), color = I("darkgrey")) + labs(title = used_genes) +
    geom_line(aes(y = dispersion_fit),  color = "red") +
    geom_point(data = subset(disp_table, gene_id %in% ordering_genes), color = "black")
  pdf(fname, height = 12, width = 16);
  print(plot_grid(plotlist = g)); print(plot_grid(plotlist = void)); graphics.off()
}

cat("Set ordering filter", as.character(Sys.time()), '\n')
cds <- setOrderingFilter(cds, ordering_genes)
cat('Genes', commas(rownames(cds)[fData(cds)$use_for_ordering]), '\n')
cat('Preprocess CDS\n')
cds <- preprocessCDS(cds, num_dim = pcs_comp, norm_method = norm_pseu[[2]], residualModelFormulaStr = resMF, verbose = TRUE)
cat(as.character(Sys.time()), '\n')

fname <- paste0('data/', fsufix, '_data.RData')
if(!file.exists(fname)){
  cat('------------ Step 2: Reduce the dimensionality of the data ------------\n')
  cds <- reduceDimension(cds, reduction_method = redtype, max_components = 2, verbose = TRUE)#, perplexity = 10)

  if(redtype %in% c('UMAP', 'tSNE')){
    cat('------------ Step 3: Partition the cells into supergroups ------------\n')
    cds <- partitionCells(cds)
  }

  cat('------------ Step 4: Learn the principal graph ------------\n')
  cds <- learnGraph(cds, RGE_method = ifelse(redtype == "DDRTree", redtype, 'SimplePPT'))

  MPP_node_ids <- get_correct_root_state(cds, cell_phenotype = catg[1])
  cat('Ordering cells, root:', MPP_node_ids, '\n')
  cds <- orderCells(cds, root_pr_nodes = MPP_node_ids)
  save(cds, file = fname)
}else{
  cat('------------ Step 2 to 5: from object ------------\n')
  cds <- theObjectSavedIn(fname)
}
fname <- paste0('data/', fsufix, '_components.RData')
comps <- data.frame(t(cds@reducedDimS)); colnames(comps) <- c('MC1', 'MC2')
comps$Pseudotime <- cds@phenoData@data$Pseudotime
head(comps)
if(!file.exists(fname)) save(comps, file = fname)

cat('------------ Step 5: Visualize the trajectory ------------\n')
pdfsize <- round(log10(nrow(annot)) * 3.5, 1)
pdf(paste0(fsufix, 'trajectory', pdfsize, '.pdf'), height = pdfsize, width = pdfsize)
void <- lapply(catg, function(catgy){
  print(plot_cell_trajectory(cds, color_by = catgy, cell_size = 2) +
    scale_color_manual(values = cell_type_color[[catgy]]))
})
print(plot_cell_trajectory(cds))
graphics.off()
pdf(paste0(fsufix, 'trajectory_size0.8.pdf'), height = pdfsize, width = pdfsize)
void <- lapply(catg, function(catgy){
  print(plot_cell_trajectory(cds, color_by = catgy, cell_size = 0.8) +
    scale_color_manual(values = cell_type_color[[catgy]]))
})
graphics.off()

if(is.numeric(topgenes)){
  topgenes <- head(ordering_genes, topgenes)
}
topgenes = c(
  'PRF1', 'GZMB', 'GNLY', 'ZBTB32', 'ZBED2', 'BTLA', 'CXCL13', 'DUSP4',
  'LAG3', 'HAVCR2', 'PDCD1', 'DOK5', 'CD70', 'IL21', 'POU2AF1', 'CD200',
  'PRDM1', 'IL2RB2', 'MAF', 'IRF7'
)
genes <- getfound(topgenes, levels(cds@featureData$gene_short_name), element = 'genes', v = T)
if(length(genes) > 0){
  fname <- paste0('genes_dot_plot_', fsufix,'.pdf')
  if(!file.exists(fname)){
    void <- lapply(catg, function(catgy){
      plot_markers_by_group(cds, genes, group_by = catgy, ordering_type = 'maximal_on_diag')
    }); ssize <- fitgrid(void)
    pdf(fname, height = ssize[2]*ssize[3]+5, width = ssize[1]*ssize[3]+5)
    print(plot_grid(plotlist = void))
    graphics.off()
  }

  cds_subset <- cds[head(genes, 7), ]
  for(thiscat in c("Pseudotime", catg)){
    pdf(paste0(fsufix, '_genes_in_pseudotime_', thiscat, '.pdf'), height = 12, width = 10)
    void <- try(print(plot_genes_in_pseudotime(cds_subset, color_by = thiscat))); cat('.')
    graphics.off()
  };

  tvar <- exprs(cds[genes, ]); tmp <- t(reducedDimS(cds))
  df <- data.frame(cbind(tmp, as.matrix(t(tvar[, rownames(tmp)]))), stringsAsFactors = F)
  colnames(df) <- c('Dim1', 'Dim2', genes)
  dff <- reshape2::melt(df, id.vars = colnames(df)[1:2])
  couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000')
  dff$log2CTS_1 <- log2(dff$value + 1)
  dir.create(paste0(fsufix, '_genes_on_dimentions'))
  for(g in genes){
    ddf <- remove.factors(dff[dff$variable == g, ])
    g1 <- ggplot(ddf, aes(x = Dim1, y = Dim2, color = log2CTS_1)) +
      geom_point(size = 0.8) + scale_color_gradientn(colours = couls) +
      facet_wrap(~ variable, ncol = fitgrid(g)[2])
    pdf(paste0(fsufix, '_genes_on_dimentions/', g, '.pdf'), height = 12, width = 12);
    print(g1)
    graphics.off()
  }
  tvar <- cbind_repcol(df[, 1:2], pData(cds)[, "Total_mRNAs", drop = F])
  pdf(paste0(fsufix, '_genes_on_dimentions/total_mRNAs.pdf'), height = 12, width = 12)
  print(ggplot(tvar, aes(Dim1, Dim2, color = log2(Total_mRNAs))) +
  geom_point() + scale_color_gradientn(colours = couls))
  graphics.off()
}

#### MAIN PROGRAM #### ---------------------------------------------------------

cat('\nDONE\nStarting time:\n', st.time, '\n')
cat('Finishing time:\n'); timestamp()
cat('\n\n*******************************************************************\n')
cat('SESSION INFO:\n') ; print(sessionInfo())
cat('*******************************************************************\n\n')
cat('Pipeline successfully finished!\n')
# Closing output writing
if(file.exists('Rplots.pdf')) file.remove('Rplots.pdf')
sink(type = 'message')
sink()
# void <- try(savehistory(file = ".Rhistory"))
