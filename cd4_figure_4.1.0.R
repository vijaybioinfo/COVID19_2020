############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -------------   Figure 4    -------------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    -------------   Figure 4    -------------    ############\n')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')

# ---> About the script.
# Version: 0
# Subversion: 2

# ---> About the paper.
# Paper name: Single-cell transcriptomic analysis of SARS-CoV-2 reactive CD4 + T cells
# Authors: Meckiff, et al, 2020.

############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############
# Script to reproduce the figures shown as part of figure 4 for the paper **Single-cell transcriptomic analysis of SARS-CoV-2 reactive CD4 + T cells**.

############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n\n')
library(Seurat)
library(ggplot2)
library(reshape2)
library(english)
library(VennDiagram)
library(stringr)
library(gtools)
library(tidyr)
library(UpSetR)
library(data.table)
library(parallel)
# Path to where the handy functions file was stored. Change accordingly.
source('/home/vfajardo/scripts/functions/R_handy_functions.R')

############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ------------   Expansion    -------------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    ------------   Expansion    -------------    ############\n')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')

############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

# 1 -------------------------------------------------------------------->
# Name: Output clonotypes types values barplots across tag values.
# Description:
# Given a set of tag values defined in the seurat object within the environment, this function outputs to an output directory the frequency for expanded and non-expanded clonotypes from the cells across those tag values.
# Arguments ------------------------->
# output.path
# tag.values - Tag values to consider.
# cells.input - Logical indicating wether the output should be considered at the clonotypes (unqiue clonotype values per tag value) or cell (non-unique clonotypes) level.
# Function:

output.barplot <- function(output.path, tag.values, cells.input=TRUE, seurat.obj, tmp.tag){
  tag.clons.nos <- sapply(X=tag.values, FUN=function(tmp.tag.val){
    val.clons.info <- seurat.obj@meta.data[(seurat.obj@meta.data[, tmp.tag]==tmp.tag.val & !is.na(seurat.obj@meta.data[, tmp.tag])), c('clonotype.tag', 'TCR.tag')]
    val.clons <- val.clons.info$clonotype.tag
    # If this is required for cell level, we don't need to restrict the analysis to unique clonotypes.
    if(!cells.input) val.clons <- unique(val.clons)
    # No matter what, we need to filter NA values.
    val.clons <- val.clons[!is.na(val.clons)]
    val.clons.no <- length(val.clons)
    # Take clonotypes with expanded tags.
    val.exp.clons <- val.clons.info$clonotype.tag[val.clons.info$TCR.tag=='Expanded']
    val.exp.clons <- val.exp.clons[!is.na(val.exp.clons)]
    val.exp.clons.no <- ifelse(test=cells.input, yes=length(val.exp.clons), no=length(unique(val.exp.clons)))
    # And same for non-expanded.
    val.non.exp.clons <- val.clons.info$clonotype.tag[val.clons.info$TCR.tag=='Non-expanded']
    val.non.exp.clons <- val.non.exp.clons[!is.na(val.non.exp.clons)]
    val.non.exp.clons.no <- ifelse(test=cells.input, yes=length(val.non.exp.clons), no=length(unique(val.non.exp.clons)))
    # Return results
    to.output <- c(val.clons.no, val.exp.clons.no, val.non.exp.clons.no)
    return(to.output)
  })
  tag.clons.nos <- as.data.frame(tag.clons.nos, stringsAsFactors=FALSE)
  # Get a backup of the table to output.
  to.output <- as.data.frame(t(tag.clons.nos), stringsAsFactors=FALSE)
  colnames(to.output) <- c('total', 'abundant.size', 'non.abundant.size')
  to.output$abundant.prop <- to.output$abundant.size / to.output$total
  to.output$non.abundant.prop <- to.output$non.abundant.size / to.output$total
  tmp.file.name <- paste0(output.path, '/TagByClonotypeRawData.csv')
  write.csv(x=to.output, file=tmp.file.name)
  # Get tidy data.
  tag.clons.nos$type <- c('All', 'Expanded', 'Non-expanded')
  # Save a bakcup of this stricture to get a specific tag values order.
  tmp.df <- tag.clons.nos
  # Then, go on tidying data up.
  tag.clons.nos <- gather(data=tag.clons.nos, key=tag.value, value=clonotypes, tag.values)
  # About the tag.
  if(grepl(x=tmp.tag, pattern='\\.tag')) tag.lab <- str_replace(string=tmp.tag, pattern='\\.tag', replacement='') else tag.lab <- tmp.tag
  # @ Clonotypes percentages (across tag values) barplots.
  # Filter certain clonotype type if required.
  tag.clons.nos <- tag.clons.nos[tag.clons.nos$type!='All', ]
  # Set tag values order according to expansion level.
  rownames(tmp.df) <- tmp.df$type
  tmp.df$type <- NULL
  tag.vals.order <- colnames(tmp.df)[order(tmp.df['Expanded', ] / tmp.df['All', ], decreasing=FALSE)]
  tag.clons.nos$tag.value <- factor(x=tag.clons.nos$tag.value, levels=tag.vals.order)
  # ---> ggplot and output.
  # Normal.
  tmp.ggplot <- ggplot(data=tag.clons.nos, aes(x=tag.value, y=clonotypes, fill=type)) + geom_bar(stat='identity', position='fill', width=0.8) + scale_y_continuous(labels = scales::percent_format()) + scale_fill_manual(values=expansion.cols) + coord_flip() + theme(panel.background=element_blank()) + xlab(tag.lab) + ylab('Clonotypes\' proportion')
  tmp.file.name <- paste0(output.path, '/TagByClonotypeExpansionFilledBarplot.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
  # Blank style
  tmp.file.name <- paste0(output.path, '/TagByClonotypeExpansionFilledBarplot_Blank.pdf')
  pdf(file=tmp.file.name, height=10)
  print(tmp.ggplot + blank.complement)
  dev.off()
  # ---> Output.
  # Return main table.
  return(to.output)
}

# 2 -------------------------------------------------------------------->
# Name: Main TCR analysis.
# Description:
# Main program. To be applied to each data set.
# Workflow:
# 1. Annotation of seurat object with TCR data based on cellranger vdj (aggr) output.
# 2. Display of cell expansion on UMAP (dim reduction).
# 3. Calculation and output of cell and clone expansion across groups for a set of tags.
# Arguments ------------------------->
# @ reports.path - Absolute path to directory to save reports.
# @ seurat.obj - Seurat object loaded in the current session.
# @ clons.info - clontypes.csv file loaded in the current session.
# @ cells.clons.info - filtered_contig_annotations.csv file loaded in the current session.
# @ tag.values - Tag values to consider.
# cells.input - Logical indicating wether the output should be considered at the clonotypes (unqiue clonotype values per tag value) or cell (non-unique clonotypes) level.
# Function:

do.tcr.analysis <- function(reports.path, this.seurat.obj, clons.info, cells.clons.info, tags.of.int){
  ### ------------------------- Main program ------------------------- ###
  cat('### ------------------------- Main program ------------------------- ###\n\n')

  ### ------------------ Mark clonotypes on UMAP/tSNE ----------------- ###
  cat('### ------------------ Mark clonotypes on UMAP/tSNE ----------------- ###\n')
  dim.reduction.path <- paste0(reports.path, '/dim_reduction')
  create.dir(dir.path=dim.reduction.path, path.desc='Dimensionality reduction')

  # ---> Add tags to seurat object.
  cat('# ---> Add tags to seurat object.\n')
  # Add TCR info to seurat object meta data.
  # For each cell, we'll get:
  #   TRA.tag, TRB.tag, TCR.tag (combination of TRA and TRB) each with three possible values: Expanded, Non-expanded or NA
  #   CMC (chain-related multiplet clonotype) tag.
  #   Chains sequences tags, for both, aa and nt, levels.

  clons.info.for.cells <- as.data.frame(t(sapply(X=Cells(this.seurat.obj), FUN=function(cell){
    # Chains info.
    what.chain <- cells.clons.info$chain[cells.clons.info$barcode==cell]
    # Clonotype ID.
    clon.id <- unique(cells.clons.info$raw_clonotype_id[cells.clons.info$barcode==cell])
    if(length(clon.id)>0){
      # Expansion info.
      if(clons.info$frequency[clons.info$clonotype_id==clon.id]>=expansion.thold) expanded <- 'Expanded' else expanded <- 'Non-expanded'
      # Chains sequences.
      TRA.nt.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_nt', chain.ptn='TRA')
      TRB.nt.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_nt', chain.ptn='TRB')
      TRA.aa.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_aa', chain.ptn='TRA')
      TRB.aa.chains <- get.cdr3s.from.table(cdr3s.table=clons.info[clons.info$clonotype_id==clon.id, ], cdr3.col='cdr3s_aa', chain.ptn='TRB')
      # Clonotype size.
      clon.size <- unique(clons.info[clons.info$clonotype_id==clon.id, 'frequency'])
      clon.prop <- unique(clons.info[clons.info$clonotype_id==clon.id, 'proportion'])
      # CMC-related info.
      if(length(TRA.nt.chains) > 1 | length(TRB.nt.chains) > 1) cmc.tag <- TRUE else cmc.tag <- FALSE
      if(length(unique(what.chain))==2){
        to.return <- c(clon.id, rep(expanded, times=3), cmc.tag, paste0(TRA.nt.chains, collapse=';'), paste0(TRB.nt.chains, collapse=';'), paste0(TRA.aa.chains, collapse=';'), paste0(TRB.aa.chains, collapse=';'), clon.size, clon.prop)
      }else{
        # The TCR tag will indicate if there's any chain marked as expanded.
        if('TRA' %in% what.chain) to.return <- c(clon.id, expanded, NA, expanded, cmc.tag, paste0(TRA.nt.chains, collapse=';'), NA, paste0(TRA.aa.chains, collapse=';'), NA, clon.size, clon.prop) else to.return <- c(clon.id, NA, expanded, expanded, cmc.tag, NA, paste0(TRB.nt.chains, collapse=';'), NA, paste0(TRB.aa.chains, collapse=';'), clon.size, clon.prop)
      }
    }else{
      to.return <- rep(NA, times=11)
    }
    names(to.return) <- c('clonotype.tag', 'TRA.tag', 'TRB.tag', 'TCR.tag', 'CMC.tag', 'TRA.nt.chains.tag', 'TRB.nt.chains.tag', 'TRA.aa.chains.tag', 'TRB.aa.chains.tag', 'clon.size.tag', 'clon.proportion.tag')
    return(to.return)
  })))
  clons.info.for.cells[, 'CMC.tag'] <- as.logical(clons.info.for.cells[, 'CMC.tag'])
  clons.info.for.cells[, 'clonotype.tag'] <- as.character(clons.info.for.cells[, 'clonotype.tag'])
  clons.info.for.cells[, 'clon.size.tag'] <- as.numeric(as.character(clons.info.for.cells[, 'clon.size.tag']))
  clons.info.for.cells[, 'clon.proportion.tag'] <- as.numeric(as.character(clons.info.for.cells[, 'clon.proportion.tag']))
  this.seurat.obj@meta.data <- cbind(this.seurat.obj@meta.data, clons.info.for.cells)

  # ---> Depict new tags on Dimensionality reduction maps.
  cat('# ---> Depict new tags in a DimPlot.\n')
  # -@ Clone size.
  cells.to.depict <- Cells(this.seurat.obj)[!is.na(this.seurat.obj@meta.data[, 'clon.size.tag']) & this.seurat.obj@meta.data[, 'clon.size.tag'] > 2 & (this.seurat.obj@meta.data[, 'orig.virus2']=='CV')]
  # Normal style.
  tmp.file.name <- paste0(dim.reduction.path, '/CloneSizeOnFeaturePlot.pdf')
  pdf(file=tmp.file.name)
  print(FeaturePlot(this.seurat.obj, reduction='umap', feature='clon.size.tag', cells=cells.to.depict, min.cutoff=3, max.cutoff=up.thold) + labs(title='Clone size') + scale_alpha(0.7) + scale_colour_gradientn(trans='log2', colours=this.color.scale, limits=c(3, up.thold)))
  dev.off()
  # Blank style.
  tmp.file.name <- paste0(dim.reduction.path, '/CloneSizeOnFeaturePlot_Blank.pdf')
  pdf(file=tmp.file.name)
  print(FeaturePlot(this.seurat.obj, reduction='umap', feature='clon.size.tag', cells=cells.to.depict, min.cutoff=3, max.cutoff=up.thold) + scale_alpha(0.7) + scale_colour_gradientn(trans='log2', colours=this.color.scale, limits=c(3, up.thold)) + theme(line=element_blank(), title=element_blank(), axis.text=element_blank()))
  dev.off()

  if(!is.null(tags.of.int)){
  ### --------------------- Tag-specific analyses --------------------- ###
  cat('### --------------------- Tag-specific analyses --------------------- ###\n')
  tags.path <- paste0(reports.path, '/tag_specific_analysis')
  create.dir(dir.path=tags.path, path.desc='Tags-specific analyses')
  # Regarding TCR info.
  for(tmp.tag in tags.of.int){
    cat('---> Process for tag', tmp.tag, '\n')
    # ------ TCR expansion accross tag values ------ #
    # Get all tag values (different than NA).
    tag.values <- mixedsort(unique(as.character(this.seurat.obj@meta.data[, tmp.tag])), decreasing=FALSE)
    tag.values <- tag.values[!is.na(tag.values)]

    if(is.factor(this.seurat.obj@meta.data[, tmp.tag])) this.seurat.obj@meta.data[, tmp.tag] <- factor(x=this.seurat.obj@meta.data[, tmp.tag], levels=tag.values)

    # Skip in case there aren't at least two tag values to evaluate.
    if(length(tag.values)<2) next

    # This tag path
    tmp.tag.path <- paste0(tags.path, '/', tmp.tag)
    create.dir(dir.path=tmp.tag.path, path.desc=tmp.tag)

    # ---> Clone size across cells in the same group.
    clones.tag.path <- paste0(tmp.tag.path, '/clone_size')
    create.dir(dir.path=clones.tag.path, path.desc=paste0(tmp.tag, ' for clone size'))
    # Box plot.
    tmp.ggplot <- ggplot(data=this.seurat.obj@meta.data, aes_string(y='clon.size.tag', x=tmp.tag, fill=tmp.tag)) + geom_boxplot(outlier.shape=NA, alpha=0.7) + scale_y_continuous(trans='log10') + theme(legend.position='none') + ylab('Clone size') + labs(title=paste0('Clone size across groups'))
    tmp.file.name <- paste0(clones.tag.path, '/CloneSizeAcrossTagValuesBoxplots.pdf')
    pdf(tmp.file.name)
    print(tmp.ggplot)
    dev.off()
    # Violin plot.
    tmp.ggplot <- ggplot(data=this.seurat.obj@meta.data, aes_string(y='clon.size.tag', x=tmp.tag, fill=tmp.tag)) + geom_violin(alpha=0.7) + scale_y_continuous(trans='log10') + theme(legend.position='none') + ylab('Clone size') + labs(title=paste0('Clone size across groups'))
    tmp.file.name <- paste0(clones.tag.path, '/CloneSizeAcrossTagValuesViolinPlots.pdf')
    pdf(tmp.file.name)
    print(tmp.ggplot)
    dev.off()

    # ---> Expansion for cells across tag values.
    tmp.cells.tag.path <- paste0(tmp.tag.path, '/cells_expansion')
    create.dir(dir.path=tmp.cells.tag.path, path.desc=paste0(tmp.tag, ' for cells expansion'))
    # Tag by clonotypes info for cells.
    # i.e., we don't care about unique clonotypes.
    tmp <- output.barplot(output.path=tmp.cells.tag.path, tag.values=tag.values, cells.input=TRUE, seurat.obj=this.seurat.obj, tmp.tag=tmp.tag)

    # ---> Expansion for clonotypes across tag values.
    tmp.clons.tag.path <- paste0(tmp.tag.path, '/clons_expansion')
    create.dir(dir.path=tmp.clons.tag.path, path.desc=paste0(tmp.tag, ' for clonotypes expansion'))
    tmp <- output.barplot(output.path=tmp.clons.tag.path, tag.values=tag.values, cells.input=FALSE, seurat.obj=this.seurat.obj, tmp.tag=tmp.tag)
  }
  cat('Tag-specific analyses finished!\n')
  }

  cat('\n\n')
  ### ------------------------- Seurat object ------------------------- ###
  cat('### ------------------------- Seurat object ------------------------- ###\n')
  objs.path <- paste0(reports.path, '/seurat_objects')
  create.dir(dir.path=objs.path, path.desc='Seurat Objects')
  tmp.file.name <- paste0(objs.path, '/SeuratObj.RDS')
  saveRDS(object=this.seurat.obj, file=tmp.file.name)
  cat('Seurat object saved!\n')

  # Return seurat object.
  return(this.seurat.obj)
}

# 3 -------------------------------------------------------------------->
# Name: Depict clusters' expansion.
# Description:
# Outputs (to dim_reduction folder for a given data set) a file depicting clone size across clusters through a Violin plot.
# Arguments ------------------------->
# @ this.seurat.obj - Seurat object loaded in the current session.
# @ clusters.tag - Tag (defined in seurat object's meta data) for cells' clusters annotations.
# Function:

# get.clusts.expansion <- function(this.seurat.obj, clusters.tag){
#   # Directory to dimensionality reduction results.
#   tmp.path <- paste0(reports.path, '/dim_reduction')
#   if(!dir.exists(tmp.path)) stop('Dimensionality reduction path does not exist.')
#   # Clone size across cells per cluster.
#   # Set data
#   tmp.df <- this.seurat.obj@meta.data
#   # Clusters in order (numeric).
#   new.levels <- mixedsort(unique(as.character(tmp.df[, clusters.tag])))
#   levels(tmp.df[, clusters.tag]) <- new.levels
#   # ggplot and output
#   tmp.file.name <- paste0(tmp.path, '/CloneSizeAcrossClustersVlnPlot.pdf')
#   tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill=clusters.tag)) + geom_violin(alpha=0.7, adjust=7/8) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme(legend.position='none', panel.background=element_blank())
#   pdf(file=tmp.file.name, width=12)
#   print(tmp.ggplot)
#   dev.off()
#   # Mock return.
#   return(NA)
# }

# 4 -------------------------------------------------------------------->
# Name: Get mean clone size per cluster.
# Description:
# Ouptuts (within the general reports folder) a table depicting mean clone size per cluster for all cells and only for CV-reactive cells.
# Arguments ------------------------->
# @ this.seurat.obj - Seurat object loaded in the current session.
# @ clusters.tag - Tag (defined in seurat object's meta data) for cells' clusters annotations.
# Function:

get.cspc <- function(reports.path, seurat.obj){
  # Temporary path.
  tmp.path <- paste0(reports.path, '/general_reports')
  if(!dir.exists(tmp.path)) dir.create(tmp.path)
  # Get and output.
  meta.data <- as.data.table(seurat.obj@meta.data)
  # Only CV-reactive cells.
  to.output <- meta.data[orig.virus2=='CV', .('Mean clone size'=mean(clon.size.tag, na.rm=TRUE)), by=.(Cluster=RNA_snn_res.0.6)]
  tmp.file.name <- paste0(tmp.path, '/MeanCloneSizePerClusterCVOnly.csv')
  fwrite(file=tmp.file.name, x=to.output)
  # All cells.
  to.output <- meta.data[, .('Mean clone size'=mean(clon.size.tag, na.rm=TRUE)), by=.(Cluster=RNA_snn_res.0.6)]
  tmp.file.name <- paste0(tmp.path, '/MeanCloneSizePerClusterAll.csv')
  fwrite(file=tmp.file.name, x=to.output)
}


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############
cat('### ----------------------- General Arguments ----------------------- ###\n')
# ---> For tag-specific analysis.
# @ Tags to depict specifically.
gen.tags.of.int <- 'hospitalization.tag'
# ---> Others.
# @ Expansio threshold.
expansion.thold <- 2
# @ Color scales and blank style complement.
this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.
blank.complement <- theme(line=element_blank(), text=element_blank(), legend.position='none')
# @ Expanded vs non-expanded colors.
# expansion.cols <- c('#fb8072', '#d9d9d9')
expansion.cols <- binary.col.scale
names(expansion.cols) <- c('Expanded', 'Non-expanded')
# @ Upper threshold for clone size
up.thold <- 100


############    -----------------------------------------    ############
############    ------------   Expansion    -------------    ############
############    ---------   6 hrs. data set    ----------    ############
############    -----------------------------------------    ############
cta('############    ------------   Expansion    -------------    ############\n')
cat('############    ---------   6 hrs. data set    ----------    ############\n')

cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# ---> File and paths definitions.
# All should be absolute paths.
# @ Reports path (where results should be deposited).
reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2/CD4s_paper/final_figures/figure_4/figure_4_panels/expansion_6hrs'
if(!dir.exists(reports.path)) dir.create(reports.path)
# @ csv file depicting relationships between clones and cells. Aggregation of filtered_contig_annotations.csv files output by cellranger vdj from each independent sample.
cells.clons.info.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-16-2020/aggr_vdj/COVID_GenEx_Batches-1-2_CD4_STIM_6/COVID_GenEx_Batches-1-2_CD4_STIM_6_05-28-2020/filtered_contig_annotations_aggr.csv'
# @ csv file describing clone's specific info. Same, file result of the aggregation of a file output by cellranger vdj from each independent sample, in this case 'clonotypes.csv'
clons.info.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-16-2020/aggr_vdj/COVID_GenEx_Batches-1-2_CD4_STIM_6/COVID_GenEx_Batches-1-2_CD4_STIM_6_05-28-2020/clonotypes_aggr.csv'
# @ Seurat object file processed for all basic workflow steps (in brief, cell normalization (LogNormalize), cell normalization, picking of variable features, clustering and dimensionality reduction for UMAP embeddings).
seurat.obj.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2_CD4_STIM_6/COVID_GenEx_Batches-1-2_CD4_STIM_6_05-28-2020_qc-std_var-25_pc_38_hto-all_stim_stim_virus_all_stimtime_six_ciros/seurat_objects/SeuratObj.RDS'
# @ Tags of interest.
tags.of.int <- c('RNA_snn_res.0.6', gen.tags.of.int)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj <- readRDS(file=seurat.obj.file)
# Check that for the seurat object, we have dim. reductions methods applied.
if(!typeof(seurat.obj@reductions$umap)=="S4") stop('There is not either UMAP or tSNE applied to the input seurat object.')
# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

### ---------------------- Data preprocessing ----------------------- ###
# ---> UMAP coordinates.
# Flip UMAP embbedings by changing UMAP 2 dimension sign.
seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_2'] <- -seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_2']

# ---> Cells-clonotypes relationships info (filtering).
# Keep track of relationships deleted due to a chain different than alpha or beta.
chains.info <- as.data.frame(table(cells.clons.info$chain))
colnames(chains.info) <- c('chain', 'freq')
# Also, keep track of how many relationships we had prefiltering.
pre.filtering.cells <- sum(Cells(seurat.obj) %in% cells.clons.info$barcode)
pre.filtering.rels <- nrow(cells.clons.info)
# Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
# Conditions:
# * For barcodes called as cells.
# * For clonotype contigs marked with high confidence.
# * For chains ultimately defined as TRA or TRB.
# * For productive clonotypes (see 10X support website for a broader definition of 'productive').
# Also, we'll take out this info since it has been taken into consideration already.
cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]

### ------------------------- Main program ------------------------- ###
cat('### ------------------------- Main program ------------------------- ###\n\n')
# ---> General workflow.
seurat.obj <- do.tcr.analysis(reports.path=reports.path, this.seurat.obj=seurat.obj, clons.info=clons.info, cells.clons.info=cells.clons.info, tags.of.int=tags.of.int)

# ---> Clone size across clusters.
get.clusts.expansion(this.seurat.obj=seurat.obj, clusters.tag='RNA_snn_res.0.6')

# ---> Mean clone size per cluster (only CV-reactive cells).
get.cspc(reports.path=reports.path, seurat.obj=seurat.obj)

############    -----------------------------------------    ############
############    ------------   Expansion    -------------    ############
############    ---------   24 hrs. data set    ---------    ############
############    -----------------------------------------    ############
cta('############    ------------   Expansion    -------------    ############\n')
cat('############    ---------   24 hrs. data set    ---------    ############\n')

### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# ---> File and paths definitions.
# All should be absolute paths.
# @ Reports path (where results should be deposited).
reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2/CD4s_paper/final_figures/figure_4/figure_4_panels/expansion_24hrs'
if(!dir.exists(reports.path)) dir.create(reports.path)
# @ csv file depicting relationships between clones and cells. Aggregation of filtered_contig_annotations.csv files output by cellranger vdj from each independent sample.
cells.clons.info.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-16-2020/aggr_vdj/COVID_GenEx_Batches-1-2_CD4_STIM_24/COVID_GenEx_Batches-1-2_CD4_STIM_24_05-28-2020/filtered_contig_annotations_aggr.csv'
# @ csv file describing clone's specific info. Same, file result of the aggregation of a file output by cellranger vdj from each independent sample, in this case 'clonotypes.csv'
clons.info.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-16-2020/aggr_vdj/COVID_GenEx_Batches-1-2_CD4_STIM_24/COVID_GenEx_Batches-1-2_CD4_STIM_24_05-28-2020/clonotypes_aggr.csv'
# @ Seurat object file processed for all basic workflow steps (in brief, cell normalization (LogNormalize), cell normalization, picking of variable features, clustering and dimensionality reduction for UMAP embeddings).
seurat.obj.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2_CD4_STIM_24/COVID_GenEx_Batches-1-2_CD4_STIM_24_05-28-2020_qc-std_var-25_pc_16_hto-all_stim_stim_virus_all_stimtime_twentyfour_ciros/seurat_objects/SeuratObj.RDS'
# @ Tags of interest.
tags.of.int <- c('RNA_snn_res.0.2', gen.tags.of.int)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj <- readRDS(file=seurat.obj.file)
# Check that for the seurat object, we have dim. reductions methods applied.
if(!typeof(seurat.obj@reductions$umap)=="S4") stop('There is not either UMAP or tSNE applied to the input seurat object.')
# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

### ---------------------- Data preprocessing ----------------------- ###
# ---> Cells-clonotypes relationships info (filtering).
# Keep track of relationships deleted due to a chain different than alpha or beta.
chains.info <- as.data.frame(table(cells.clons.info$chain))
colnames(chains.info) <- c('chain', 'freq')
# Also, keep track of how many relationships we had prefiltering.
pre.filtering.cells <- sum(Cells(seurat.obj) %in% cells.clons.info$barcode)
pre.filtering.rels <- nrow(cells.clons.info)
# Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
# Conditions:
# * For barcodes called as cells.
# * For clonotype contigs marked with high confidence.
# * For chains ultimately defined as TRA or TRB.
# * For productive clonotypes (see 10X support website for a broader definition of 'productive').
# Also, we'll take out this info since it has been taken into consideration already.
cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]

cat('\n\n')
### ------------------------- Main program ------------------------- ###
cat('### ------------------------- Main program ------------------------- ###\n\n')
# General workflow.
seurat.obj <- do.tcr.analysis(reports.path=reports.path, this.seurat.obj=seurat.obj, clons.info=clons.info, cells.clons.info=cells.clons.info, tags.of.int=tags.of.int)

# Clone size across clusters.
get.clusts.expansion(this.seurat.obj=seurat.obj, clusters.tag='RNA_snn_res.0.2')

############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ----------   Clones' sharing    ---------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    ----------   Clones\' sharing    ---------    ############')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')

### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# ---> vdj data.
aggr.vdj.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-16-2020/aggr_vdj/COVID_GenEx_Batches-1-2_CD4_STIM/COVID_GenEx_Batches-1-2_CD4_STIM_05-29-2020'
cells.clons.info.file <- paste0(aggr.vdj.path, '/filtered_contig_annotations_aggr.csv')
clons.info.file <- paste0(aggr.vdj.path, '/clonotypes_aggr.csv')
file.exists(cells.clons.info.file)
file.exists(clons.info.file)
# ---> Seurat objects.
seurat.obj.six.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2_CD4_STIM_6/COVID_GenEx_Batches-1-2_CD4_STIM_6_05-28-2020_qc-std_var-25_pc_38_hto-all_stim_stim_virus_all_stimtime_six_ciros/seurat_objects/SeuratObj.RDS'
seurat.obj.tf.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2_CD4_STIM_24/COVID_GenEx_Batches-1-2_CD4_STIM_24_05-28-2020_qc-std_var-25_pc_16_hto-all_stim_stim_virus_all_stimtime_twentyfour_ciros/seurat_objects/SeuratObj.RDS'
file.exists(seurat.obj.six.file)
file.exists(seurat.obj.tf.file)
# ---> Program-specific parameters.
cells.clons.file.suffix <- 16:20
names(cells.clons.file.suffix) <- 1:5
exp.thold <- 2
clusters.tag.six <- 'RNA_snn_res.0.6'
clusters.tag.tf <- 'RNA_snn_res.0.2'
donor.id.tag <- 'donor.id.tag'
# Clusters to look at.
six.cluster.sets <- paste0(paste0('cluster.', 0:12), '.six')
tf.cluster.sets <- paste0(paste0('cluster.', 0:5), '.twentyfour')
# ---> Path definitons.
reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2/CD4s_paper/final_figures/figure_4/figure_4_panels/clones_sharing'
dir.exists(reports.path)
#dir.create(reports.path)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj.six <- readRDS(file=seurat.obj.six.file)
seurat.obj.tf <- readRDS(file=seurat.obj.tf.file)
# Check that the tags of interest are already defined in seurat objects.
all(c(donor.id.tag, clusters.tag.tf) %in% colnames(seurat.obj.tf[[]])) & all(c(donor.id.tag, clusters.tag.six) %in% colnames(seurat.obj.six[[]]))

# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
cat('### ---------------------- Data Preprocessing ---------------------- ###\n')
# ---> Seurat object.
tf.cell.names <- Cells(seurat.obj.tf)
for(suffix in names(cells.clons.file.suffix)) tf.cell.names <- str_replace(string=tf.cell.names, pattern=paste0('-', suffix), replacement=paste0('-', cells.clons.file.suffix[suffix]))
seurat.obj.tf <- RenameCells(object=seurat.obj.tf, new.names=tf.cell.names)

# ---> Clonotypes info.
# Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
# Conditions:
# * For barcodes called as cells.
# * For clonotype contigs marked with high confidence.
# * For chains ultimately defined as TRA or TRB.
# * For productive clonotypes (see 10X support website for a broader definition of 'productive').
# Also, we'll take out this info since it has been taken into consideration already.
cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]

cat('Data has been preprocessed!\n')

cat('\n\n')
### ------------------------- Main program ------------------------- ###
cat('### ------------------------- Main program ------------------------- ###\n\n')

cat('\n\n')
### -------------------------- Clones info ------------------------- ###
### ----------------------- Tidying data up! ----------------------- ###
cat('### -------------------------- Clones info ------------------------- ###\n')

# ---> Combine TCR data from both data sets.
# Take the overall metadata
meta.data <- bind_rows(seurat.obj.six@meta.data, seurat.obj.tf@meta.data, .id='data.set')
meta.data[, 'data.set'] <- ifelse(test=meta.data[, 'data.set']==1, yes='six', no='twentyfour')
rownames(meta.data) <- c(rownames(seurat.obj.six@meta.data), rownames(seurat.obj.tf@meta.data))
meta.data$barcode <- rownames(meta.data)
# Then, combine clones data with other meta data.
cells.clons.info <- as.data.table(merge(x=cells.clons.info, y=meta.data, by='barcode', all.x=TRUE, all.y=FALSE))

# ---> Clones info table.
# Filter out relationships with no clustering info available (no genex info).
cells.clons.info.tmp <- cells.clons.info[!(is.na(cells.clons.info[, RNA_snn_res.0.6])|is.na(cells.clons.info[, RNA_snn_res.0.2])), ]
# Create a clusters tag that reflects what data set the cluster comes from.
cells.clons.info.tmp[, clusters.tag:=RNA_snn_res.0.6]
tmp.idxs <- cells.clons.info.tmp[, data.set]!='six'
cells.clons.info.tmp[tmp.idxs, clusters.tag:=RNA_snn_res.0.2]
cells.clons.info.tmp[1:100, .(data.set, clusters.tag, RNA_snn_res.0.6, RNA_snn_res.0.2)]
cells.clons.info.tmp[, clusters.tag:=paste(clusters.tag, data.set, sep='.')]
# Get clone size across clusters as well as other tags info.
clons.table <- cells.clons.info.tmp[, .SD[, .(N=uniqueN(barcode)), by=.(cluster=clusters.tag, donor=donor.id.tag, peptide=cv.peptide.to.stim.tag, virus=virus.tag, vaccinne=vaccination.tag)], by=.(clone.id=raw_clonotype_id)]
# Spread according to cluster value.
clons.table <- spread(data=clons.table, key=cluster, value=N, fill=0, sep='.')
clons.table <- as.data.frame(clons.table)
# Iterate per unique clone ID to compute sharing between groups for each tag.
clone.ids <- unique(clons.table[, 'clone.id'])
cluster.cols <- colnames(clons.table)[grepl(x=colnames(clons.table), pattern='cluster.')]
tag.cols <- colnames(clons.table)[!grepl(x=colnames(clons.table), pattern='cluster.')]
clons.table <- lapply(X=clone.ids, FUN=function(clone.id){
  # --> Clone-specific info.
  idxs <- clons.table[, 'clone.id']==clone.id
  clone.info <- clons.table[idxs, ]
  # --> Get tags recovered for each clone.
  tags.info <- apply(X=clone.info[, tag.cols], MARGIN=2, FUN=function(tmp.tag){ tmp.tag <- unique(tmp.tag); tmp.tag <- tmp.tag[!is.na(tmp.tag)]; tmp.tag <- sort(tmp.tag); paste0(tmp.tag, collapse=';') })
  # --> Total size per individual.
  clusters.info <- apply(X=clone.info[, cluster.cols], MARGIN=2, FUN=sum)
  # --> Clone size.
  clone.size <- sum(clusters.info)
  names(clone.size) <- 'clone.size'
  # --> Format and output.
  clone.info <- c(tags.info, clusters.info, clone.size)
  return(clone.info)
})
clons.table <- Reduce(x=clons.table, f=rbind)
clons.table <- as.data.frame(clons.table)

for(tmp.thold in seq(from=100, to=1000, by=100)){
# ---> Save tidy data.
# And we reload in order to capture integer info.
tmp.file.name <- paste0(reports.path, '/ClonesTidyData.csv')
write.csv(file=tmp.file.name, x=clons.table, quote=FALSE, row.names=FALSE)
clons.table <- read.csv(file=tmp.file.name, stringsAsFactors=FALSE)

six.cluster.sets <- paste0(paste0('cluster.', 0:12), '.six')
tf.cluster.sets <- paste0(paste0('cluster.', 0:5), '.twentyfour')

cat('\n\n')
############    ----------   Clones' sharing    ---------    ############
############    ----------   on a cell basis    ---------    ############
cat('############    ----------   Clones\' sharing    ---------    ############\n')

 # ----> Data loading, arguments definition and preprocessing.
 # Data.
tmp.file.name <- paste0(reports.path, '/ClonesTidyData.csv')
clons.table <- read.csv(file=tmp.file.name, stringsAsFactors=FALSE)
# Agruments.
six.cluster.sets <- paste0(paste0('cluster.', 0:12), '.six')
tf.cluster.sets <- paste0(paste0('cluster.', 0:5), '.twentyfour')
# Clone size per data set.
clons.table$clone.size.six <- rowSums(clons.table[, six.cluster.sets])
clons.table$clone.size.twentyfour <- rowSums(clons.table[, tf.cluster.sets])
# Remove non-COVID reactive TCRs.
tmp.table <- str_split(string=clons.table$virus, pattern=';', simplify=TRUE)
row.to.keep <- apply(X=tmp.table, MARGIN=1, FUN=function(tmp.clone) return('CV' %in% tmp.clone))
clons.table <- clons.table[row.to.keep, ]

# ---> Set elements per cluster according to clone size.
# Elements and relabeling.
for(col.idx in c(six.cluster.sets, tf.cluster.sets)){
  col.vals <- clons.table[, col.idx]
  col.vals[col.vals <= 1] <- 0
  col.vals[col.vals > 1] <- 1
  col.vals <- as.integer(col.vals)
  # Add new label.
  if(str_detect(string=col.idx, pattern='six')) new.lab <- paste0('6.C', str_extract(string=col.idx, pattern='\\d+')) else new.lab <- paste0('24.C', str_extract(string=col.idx, pattern='\\d+'))
  clons.table[, new.lab] <- col.vals
  # Remove old label.
  clons.table[, col.idx] <- NULL
}
# New labels.
six.cluster.sets <- paste0('6.C', str_extract(string=six.cluster.sets, pattern='\\d+'))
tf.cluster.sets <- paste0('24.C', str_extract(string=tf.cluster.sets, pattern='\\d+'))

# ---> Create clone-basis structure.
# Remove clones not seen in the 6H data set with a clone size larger than 1 for at least one of the clusters.
clons.to.keep <- apply(X=clons.table, MARGIN=1, FUN=function(clone.info){
  to.eval <- sum(clone.info[six.cluster.sets]>0)
  return(to.eval>0)
})
clons.table <- clons.table[clons.to.keep, ]
# Then, get cell-basis structure of the data. In brief, a clone row is repeated as many times as its clone size.
clones.info.cb <- mclapply(X=1:nrow(clons.table), FUN=function(clone.row.id){
  clones.clusts.info <- clons.table[clone.row.id, six.cluster.sets]
  times.to.rep <- clons.table[clone.row.id, 'clone.size.six']
  clones.clusts.info <- mclapply(X=1:times.to.rep, FUN=function(x) return(clones.clusts.info))
  clones.clusts.info <- rbindlist(clones.clusts.info)
  return(clones.clusts.info)
})
clones.info.cb <- rbindlist(clones.info.cb)
# Proof checking we got it right.
# Total clone size should equal the total amount of rows we got in the processed data.
sum(clons.table[, 'clone.size.six'])
nrow(clones.info.cb)

# ---> Clone size definition across intersections.
# This took me ages, damn, it's not gonna be in the paper!
# # Define clone-specific sets.
# clone.sets <- colnames(clones.info.cb)
#
# # Get clone sets frequency from the cell basis info.
# clone.sets.cb.freq <- apply(X=clones.info.cb, MARGIN=1, FUN=function(clone.counts){
#   to.ouput <- clone.sets[clone.counts>0]
#   to.ouput <- paste0(to.ouput, collapse=';')
#   return(to.ouput)
# })
# clone.sets.cb.freq <- as.data.frame(table(clone.sets.cb.freq), stringsAsFactors=FALSE)
#
# # Get clone sets frequency from the clone basis info.
# clone.sets.tb.freq <- apply(X=clons.table, MARGIN=1, FUN=function(clone.info){
#   clone.counts <- clone.info[clone.sets]
#   to.ouput <- clone.sets[clone.counts>0]
#   to.ouput <- paste0(to.ouput, collapse=';')
#   return(to.ouput)
# })
# # Add it to the original table.
# clons.table$clones.set <- clone.sets.tb.freq
# # Subset info to keep clone size and set temp ID (and only for the clonotypes of interest -freq. > 1-).
# short.clons.table <- clons.table[, c('clone.size.six', 'clones.set')]
#
# # Also add degree.
# clone.sets.info <- clone.sets.cb.freq
# colnames(clone.sets.info) <- c('set', 'freq')
# clone.sets.info$degree <- (str_count(string=clone.sets.info$set, pattern=';') + 1)
#
# # Set order and set ID.
# clone.sets.info <- as.data.table(clone.sets.info)
# setorderv(x=clone.sets.info, c('degree', 'freq'), order=c(1, -1))
# clone.set.ids <- as.character(1:(nrow(clone.sets.info)))
# names(clone.set.ids) <- clone.sets.info[, set]
#
# # Add clone ID to clonotype-specific info.
# short.clons.table$clone.set.id <- factor(x=clone.set.ids[short.clons.table$clones.set], levels=clone.set.ids)
#
# tmp.ggplot <- ggplot(data=short.clons.table, aes(x=clone.set.id, y=clone.size.six)) + geom_boxplot(fill='#0074e6') + theme(axis.text.x=element_text(angle=80)) + scale_y_log10() + xlab('Intersection ID') + ylab('Clone size (log10)')
# tmp.file.name <- paste0(reports.path, '/ClonesSharingAcrossClusters6HrsDataSetOnCellBasisClonalExpansion.pdf')
# pdf(file=tmp.file.name, width=12)
# print(tmp.ggplot)
# dev.off()

# ---> Clonal sharing across clusters.
# ---> 6 hrs. clusters.
# Normal version.
for(tmp.scale in c('log10', 'log2', 'identity')){
  tmp.file.name <- paste0(reports.path, '/ClonesSharingAcrossClusters6HrsDataSetOnCellBasisScale', tmp.scale, '.pdf')
  pdf(file=tmp.file.name, width=12)
  print(upset(data=clones.info.cb,
    sets=six.cluster.sets,
    nintersects=NA, # Manual counting (when using NA as the parameter value) allowed us to determine that there are only 32 intersection sets with an intersection size larger than 2
    query.legend="top",
    order.by=c("freq", "degree"),
    keep.order=TRUE,
    # mainbar.y.label='Clones Intersection Size',
    # sets.x.label="Total clones (freq. > 2)",
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    # text.scale=c(0, 0, 0, 0, 0, 0), # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
    scale.intersections=tmp.scale
  ))
  dev.off()
}

# Blank version.
for(tmp.scale in c('log10', 'log2', 'identity')){
  tmp.file.name <- paste0(reports.path, '/ClonesSharingAcrossClusters6HrsDataSetOnCellBasisScale', tmp.scale, '_Blank.pdf')
  pdf(file=tmp.file.name, width=12)
  print(upset(data=clones.info.cb,
    sets=six.cluster.sets,
    nintersects=NA, # Manual counting (when using NA as the parameter value) allowed us to determine that there are only 32 intersection sets with an intersection size larger than 2
    query.legend="top",
    order.by=c("freq", "degree"),
    keep.order=TRUE,
    # mainbar.y.label='Clones Intersection Size',
    # sets.x.label="Total clones (freq. > 2)",
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    text.scale=c(0, 0, 0, 0, 0, 0), # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
    scale.intersections=tmp.scale
  ))
  dev.off()
}

############    -----------------------------------------    ############
### ----------------------------- Cites ----------------------------- ###
############    -----------------------------------------    ############
# ---> UpSetR
# Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R Package for the Visualization of Intersecting Sets and their Properties doi: https://doi.org/10.1093/bioinformatics/btx364
