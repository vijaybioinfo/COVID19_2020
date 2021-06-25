############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -------------   Figure 3    -------------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    -------------   Figure 3    -------------    ############\n')
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

cat('\n\n')
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
# Name: Set TCR annotations for seurat object.
#   First part of the main TCR analysis.
# Description:
# Workflow:
# 1. Annotation of seurat object with TCR data based on cellranger vdj (aggr) output. (HERE)
# 2. Display of cell expansion on UMAP (dim reduction). (NOT HERE)
# 3. Calculation and output of cell and clone expansion across groups for a set of tags. (NO HERE)
# Arguments ------------------------->
# @ reports.path - Absolute path to directory to save reports.
# @ this.seurat.obj - Seurat object loaded in the current session.
# @ clons.info - clontypes.csv file loaded in the current session.
# @ cells.clons.info - filtered_contig_annotations.csv file loaded in the current session.
# Function:

do.tcr.annotations <- function(reports.path, this.seurat.obj, clons.info, cells.clons.info){
  ### ------------------------- Main program ------------------------- ###
  cat('### ------------------------- Main program ------------------------- ###\n\n')

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
# Name: TCR analysis: Clonal expansion and tag-specific analysis.
#   Second part of the main TCR analysis.
# Description:
# Main program. To be applied to each data set.
# Workflow:
# 1. Annotation of seurat object with TCR data based on cellranger vdj (aggr) output. (NOT HERE)
# 2. Display of cell expansion on UMAP (dim reduction). (HERE)
# 3. Calculation and output of cell and clone expansion across groups for a set of tags. (HERE)
# Arguments ------------------------->
# @ reports.path - Absolute path to directory to save reports.
# @ this.seurat.obj - Seurat object loaded in the current session.
# @ tag.values - Tag values to consider.
# Function:

do.tcr.analysis <- function(reports.path, this.seurat.obj, tags.of.int){
  ### ------------------------- Main program ------------------------- ###
  cat('### ------------------------- Main program ------------------------- ###\n\n')

  ### ------------------ Mark clonotypes on UMAP/tSNE ----------------- ###
  cat('### ------------------ Mark clonotypes on UMAP/tSNE ----------------- ###\n')
  dim.reduction.path <- paste0(reports.path, '/dim_reduction')
  create.dir(dir.path=dim.reduction.path, path.desc='Dimensionality reduction')

  # Depict clone size distribution.
  tmp.ggplot <- ggplot(data=this.seurat.obj@meta.data, aes(x=clon.size.tag)) + geom_density(alpha=0.8, fill='#8b0000') + scale_x_log10() + theme_minimal() + xlab('Clone size') + ylab('Density')
  tmp.file.name <- paste0(dim.reduction.path, '/CloneSizeDistribution.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()

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

  # Mock return.
  return(NA)
}

# 4 -------------------------------------------------------------------->
# Name: Depict clusters' expansion.
# Description:
# Outputs (to dim_reduction folder for a given data set) a file depicting clone size across clusters through a Violin plot.
# Arguments ------------------------->
# @ this.seurat.obj - Seurat object loaded in the current session.
# @ clusters.tag - Tag (defined in seurat object's meta data) for cells' clusters annotations.
# Function:

get.clusts.expansion <- function(this.seurat.obj, clusters.tag, reports.path, clusters.to.keep=NULL, covid.only=FALSE, adjust.val=1.6, size.thold=10, width.val=12){
  # ---> Set data.
  # Clone size across cells per cluster.
  # Set data
  tmp.df <- this.seurat.obj@meta.data
  # Clusters in order (numeric).
  new.levels <- mixedsort(unique(as.character(tmp.df[, clusters.tag])))
  tmp.df[, clusters.tag] <- factor(x=tmp.df[, clusters.tag], levels=new.levels)
  # Turn dataframe into datatable.
  tmp.df <- as.data.table(tmp.df)
  # Create a new tag.
  tmp.df[, cluster.hospitalization.tag:=paste0(get(clusters.tag), '.', hospitalization.tag)]

  # ---> Expansion proportions per cluster.
  # First, get rid of cells from viruses different than COVID if necessary.
  if(covid.only) tmp.df <- tmp.df[virus.tag=='CV19' & !is.na(hospitalization.tag)]
  # @ Add proportion of clonally expanded cells per cluster.
  # For all of the cells from each cluster.
  cell.props <- tmp.df[, .(all.exp.cells.prop=sum(clon.size.tag>1, na.rm=TRUE)/.N), by=.(cluster=get(clusters.tag))]
  colnames(cell.props)[1] <- clusters.tag
  tmp.df <- merge(x=tmp.df, y=cell.props, by=clusters.tag, all.x=TRUE)
  # For cells from each cluster-hospital status pair.
  cell.props <- tmp.df[, .(hosp.exp.cells.prop=sum(clon.size.tag>1, na.rm=TRUE)/.N), by=.(cluster=cluster.hospitalization.tag)]
  colnames(cell.props)[1] <- 'cluster.hospitalization.tag'
  tmp.df <- merge(x=tmp.df, y=cell.props, by='cluster.hospitalization.tag', all.x=TRUE)
  # @ Add proportion of clonally expanded clones per cluster.
  # For all of the cells from each cluster.
  clones.props <- tmp.df[, .(no.cells=.N), by=.(cluster=get(clusters.tag), clone=clonotype.tag)]
  clones.props <- clones.props[, .(all.exp.clones.prop=sum(no.cells>1)/.N), by=cluster]
  colnames(clones.props)[1] <- clusters.tag
  tmp.df <- merge(x=tmp.df, y=clones.props, by=clusters.tag, all.x=TRUE)
  # For cells from each cluster-hospital status pair.
  clones.props <- tmp.df[, .(no.cells=.N), by=.(cluster=cluster.hospitalization.tag, clone=clonotype.tag)]
  clones.props <- clones.props[, .(hosp.exp.clones.prop=sum(no.cells>1)/.N), by=cluster]
  colnames(clones.props)[1] <- 'cluster.hospitalization.tag'
  tmp.df <- merge(x=tmp.df, y=clones.props, by='cluster.hospitalization.tag', all.x=TRUE)
  # @ Add clone size mean across clusters.
  # For all of the cells from each cluster.
  cs.means <- tmp.df[, .(all.clone.size.mean=mean(clon.size.tag, na.rm=TRUE)), by=.(cluster=get(clusters.tag))]
  colnames(cs.means)[1] <- clusters.tag
  tmp.df <- merge(x=tmp.df, y=cs.means, by=clusters.tag, all.x=TRUE)
  # For cells from each cluster-hospital status pair.
  cs.means <- tmp.df[, .(hosp.clone.size.mean=mean(clon.size.tag, na.rm=TRUE)), by=.(cluster=cluster.hospitalization.tag)]
  colnames(cs.means)[1] <- 'cluster.hospitalization.tag'
  tmp.df <- merge(x=tmp.df, y=cs.means, by='cluster.hospitalization.tag', all.x=TRUE)
  # @ Add clone size mean across clusters.
  # For all of the cells from each cluster.
  cs.medians <- tmp.df[, .(all.clone.size.median=median(clon.size.tag, na.rm=TRUE)), by=.(cluster=get(clusters.tag))]
  colnames(cs.medians)[1] <- clusters.tag
  tmp.df <- merge(x=tmp.df, y=cs.medians, by=clusters.tag, all.x=TRUE)
  # For cells from each cluster-hospital status pair.
  cs.medians <- tmp.df[, .(hosp.clone.size.median=median(clon.size.tag, na.rm=TRUE)), by=.(cluster=cluster.hospitalization.tag)]
  colnames(cs.medians)[1] <- 'cluster.hospitalization.tag'
  tmp.df <- merge(x=tmp.df, y=cs.medians, by='cluster.hospitalization.tag', all.x=TRUE)

  # ---> Filter cells.
  # Get rid of cells with no TCR data.
  tmp.df <- tmp.df[!is.na(clon.size.tag)]
  # If necessary, get rid of cells coming from a specific set of clusters.
  if(!is.null(clusters.to.keep))  tmp.df <- tmp.df[get(clusters.tag) %in% clusters.to.keep]

  # ---> Plots
  # @ General violin plot.
  # Across clusters.
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersVlnPlot.pdf')
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill=clusters.tag)) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + theme(legend.position='none')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # Across clusters-disease severity pairs.
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseStatusVlnPlot.pdf')
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill=clusters.tag)) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + theme(legend.position='bottom')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()

  # @ Plot depicting Cell expansion.
  # Across clusters.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill='all.exp.cells.prop')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, all.exp.cells.prop], probs=c(0.1, 0.9)), colors=this.color.scale)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersVlnPlotWithExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # Across clusters-disease severity pairs.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.exp.cells.prop')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, hosp.exp.cells.prop], probs=c(0.1, 0.9)), colors=this.color.scale)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseSeverityVlnPlotWithExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()

  # @ Plot depicting Cell expansion.
  # Across clusters.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill='all.exp.clones.prop')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, all.exp.clones.prop], probs=c(0.1, 0.9)), colors=this.color.scale)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersVlnPlotWithExpandedClonesProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # Across clusters-disease severity pairs.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.exp.clones.prop')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, hosp.exp.clones.prop], probs=c(0.1, 0.9)), colors=this.color.scale)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseSeverityVlnPlotWithExpandedClonesProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()

  # @ Plot depicting mean clone size.
  max.y <- tmp.df[, max(clon.size.tag)]
  # - Across clusters.
  # Proportion of clonally expanded cells across cluster.
  cell.props <- unique(tmp.df[, .(cluster=get(clusters.tag), all.exp.cells.prop=round(all.exp.cells.prop, digits=3))] )
  colnames(cell.props)[1] <- clusters.tag
  # Plot.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill='all.clone.size.mean')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE, position = position_dodge(width = 0.9)) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, all.clone.size.mean], probs=c(0.1, 0.9)), colors=this.color.scale) + geom_text(data=cell.props, aes_string(x=clusters.tag, y=max.y, label='all.exp.cells.prop'), inherit.aes=FALSE, color='darkblue', fontface='bold', size=6) + geom_point(position = position_jitterdodge(seed=1, dodge.width=0.9), alpha=0.5, size=0.5)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersVlnPlotWithCloneSizeMean_And_ExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # - Across clusters-disease severity pairs.
  # Proportion of clonally expanded cells across cluster.
  cell.props <- unique(tmp.df[, .(cluster=cluster.hospitalization.tag, hosp.exp.cells.prop=round(hosp.exp.cells.prop, digits=3))] )
  colnames(cell.props)[1] <- 'cluster.hospitalization.tag'
  # Plot.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.clone.size.mean')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE, position = position_dodge(width = 0.9)) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10(breaks=seq(from=0, to=max(tmp.df$clon.size.tag), length.out=3)) + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, hosp.clone.size.mean], probs=c(0.1, 0.9)), colors=this.color.scale) + geom_text(data=cell.props, aes_string(x='cluster.hospitalization.tag', y=max.y, label='hosp.exp.cells.prop'), inherit.aes=FALSE, color='darkblue', fontface='bold', size=6) + geom_point(position = position_jitterdodge(seed=1, dodge.width=0.9), alpha=0.5, size=0.5)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseSeverityVlnPlotWithCloneSizeMean_And_ExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # Plot in blank version.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.clone.size.mean')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10(breaks=c(min(tmp.df$clon.size.tag), median(tmp.df$clon.size.tag), max(tmp.df$clon.size.tag))) + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + scale_fill_gradientn(limits=quantile(x=tmp.df[, hosp.clone.size.mean], probs=c(0.1, 0.9)), colors=this.color.scale) + theme(text=element_blank(), panel.background=element_blank())
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseSeverityVlnPlotWithCloneSizeMean_And_ExpandedCellsProp_Blank.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()

  # @ Plot depicting mean clone size in a clone-size-specific manner.
  # - Get data.
  max.y <- tmp.df[, max(clon.size.tag)]
  tmp.df.2 <- tmp.df[clon.size.tag>size.thold, .(clonotype.tag, clon.size.tag, cluster.hospitalization.tag, hosp.exp.cells.prop)]
  tmp.df.2 <- unique(tmp.df.2)
  # colnames(tmp.df.2)[3] <- clusters.tag
  # # @ Add clone size mean across clusters.
  # # For all of the cells from each cluster.
  # cs.means <- tmp.df.2[, .(all.clone.size.mean=mean(clon.size.tag, na.rm=TRUE)), by=.(cluster=get(clusters.tag))]
  # colnames(cs.means)[1] <- clusters.tag
  # tmp.df.2 <- merge(x=tmp.df.2, y=cs.means, by=clusters.tag, all.x=TRUE)
  # For cells from each cluster-hospital status pair.
  cs.means <- tmp.df.2[, .(hosp.clone.size.mean=mean(clon.size.tag, na.rm=TRUE)), by=.(cluster=cluster.hospitalization.tag)]
  colnames(cs.means)[1] <- 'cluster.hospitalization.tag'
  tmp.df.2 <- merge(x=tmp.df.2, y=cs.means, by='cluster.hospitalization.tag', all.x=TRUE)
  # - Across clusters-disease severity pairs.
  # Proportion of clonally expanded cells across cluster.
  cell.props <- unique(tmp.df.2[, .(cluster=cluster.hospitalization.tag, hosp.exp.cells.prop=round(hosp.exp.cells.prop, digits=3))] )
  colnames(cell.props)[1] <- 'cluster.hospitalization.tag'
  # Plot.
  tmp.ggplot <- ggplot(data=tmp.df.2, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.clone.size.mean')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE, position = position_dodge(width = 0.9)) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df.2[, hosp.clone.size.mean], probs=c(0.1, 0.9)), colors=this.color.scale) + geom_text(data=cell.props, aes_string(x='cluster.hospitalization.tag', y=max.y, label='hosp.exp.cells.prop'), inherit.aes=FALSE, color='darkblue', fontface='bold', size=6) + geom_point(position = position_jitterdodge(seed=1, dodge.width=0.9), alpha=0.5, size=0.5)
  tmp.file.name <- paste0(reports.path, '/CloneSizeCloneBasisAcrossClustersAndDiseaseSeverityVlnPlotWithCloneSizeMean_And_ExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()

  # @ Plot depicting median clone size.
  max.y <- tmp.df[, max(clon.size.tag)]
  # - Across clusters.
  # Proportion of clonally expanded cells across cluster.
  cell.props <- unique(tmp.df[, .(cluster=get(clusters.tag), all.exp.cells.prop=round(all.exp.cells.prop, digits=3))] )
  colnames(cell.props)[1] <- clusters.tag
  # Plot.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill='all.clone.size.median')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, all.clone.size.median], probs=c(0.1, 0.9)), colors=this.color.scale) + geom_text(data=cell.props, aes_string(x=clusters.tag, y=max.y, label='all.exp.cells.prop'), inherit.aes=FALSE, color='darkblue', fontface='bold', size=6)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersVlnPlotWithCloneSizeMedian_And_ExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # - Across clusters-disease severity pairs.
  # Proportion of clonally expanded cells across cluster.
  cell.props <- unique(tmp.df[, .(cluster=cluster.hospitalization.tag, hosp.exp.cells.prop=round(hosp.exp.cells.prop, digits=3))] )
  colnames(cell.props)[1] <- 'cluster.hospitalization.tag'
  # Plot.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.clone.size.median')) + geom_violin(width=1.5, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10(breaks=c(min(tmp.df$clon.size.tag), median(tmp.df$clon.size.tag), max(tmp.df$clon.size.tag))) + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme_minimal() + scale_fill_gradientn(limits=quantile(x=tmp.df[, hosp.clone.size.median], probs=c(0.1, 0.9)), colors=this.color.scale) + geom_text(data=cell.props, aes_string(x='cluster.hospitalization.tag', y=max.y, label='hosp.exp.cells.prop'), inherit.aes=FALSE, color='darkblue', fontface='bold', size=6)
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseSeverityVlnPlotWithCloneSizeMedian_And_ExpandedCellsProp.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()
  # Plot in blank version.
  tmp.ggplot <- ggplot(data=tmp.df, aes_string(x='cluster.hospitalization.tag', y='clon.size.tag', fill='hosp.clone.size.median')) + geom_violin(width=1, alpha=0.7, adjust=adjust.val, trim=FALSE) + geom_boxplot(width=0.04, alpha=0.7) + scale_y_log10(breaks=c(min(tmp.df$clon.size.tag), median(tmp.df$clon.size.tag), max(tmp.df$clon.size.tag))) + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + scale_fill_gradientn(limits=quantile(x=tmp.df[, hosp.clone.size.median], probs=c(0.1, 0.9)), colors=this.color.scale) + theme(text=element_blank(), panel.background=element_blank(), axis.line=element_line())
  tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersAndDiseaseSeverityVlnPlotWithCloneSizeMedian_And_ExpandedCellsProp_Blank.pdf')
  pdf(file=tmp.file.name, width=width.val)
  print(tmp.ggplot)
  dev.off()

  # # Violin plot with summary stats and output
  # tmp.file.name <- paste0(reports.path, '/CloneSizeAcrossClustersVlnPlotWithStatSumm.pdf')
  # tmp.ggplot <- ggplot(data=tmp.df, aes_string(x=clusters.tag, y='clon.size.tag', fill=clusters.tag)) + geom_violin(width=1.5, alpha=0.7) + stat_summary(fun.data=mean_sdl, fun.args=list(mult=1)) + stat_summary(geom='point', fun.y=mean, size=3, col='black', shape='X') + scale_y_log10() + labs(title='Clonse Size Across Clusters') + xlab('Clusters') + ylab('Clone size') + theme(legend.position='none', panel.background=element_blank())
  # pdf(file=tmp.file.name, width=width.val)
  # print(tmp.ggplot)
  # dev.off()

  # ---> Significance test: Mild vs Severe (hospital vs non-hospital).
  # @ Non-paired Wilcoxon Test.
  # Including each cell expansion value.
  mild.sizes <- tmp.df[hospitalization.tag=='No', clon.size.tag]
  severe.sizes <- tmp.df[hospitalization.tag=='Yes', clon.size.tag]
  cells.test.results <- wilcox.test(x=severe.sizes, y=mild.sizes, alternative='greater')
  # Including each unique clone expansion value.
  mild.sizes <- unique(tmp.df[hospitalization.tag=='No', .(clonotype.tag, clon.size.tag)])
  mild.sizes <- mild.sizes[, clon.size.tag]
  severe.sizes <- unique(tmp.df[hospitalization.tag=='Yes', .(clonotype.tag, clon.size.tag)])
  severe.sizes <- severe.sizes[, clon.size.tag]
  clones.test.results <- wilcox.test(x=severe.sizes, y=mild.sizes, alternative='greater')
  # Format and output.
  cells.significance <- ifelse(test=cells.test.results$p.value>0.05, yes=NA,
    no=ifelse(test=cells.test.results$p.value>0.01, yes='*',
    no=ifelse(test=cells.test.results$p.value>0.001, yes='**',
    no=ifelse(test=cells.test.results$p.value>0.0001, yes='***', no='****'))))
  clones.significance <- ifelse(test=clones.test.results$p.value>0.05, yes=NA,
    no=ifelse(test=clones.test.results$p.value>0.01, yes='*',
    no=ifelse(test=clones.test.results$p.value>0.001, yes='**',
    no=ifelse(test=clones.test.results$p.value>0.0001, yes='***', no='****'))))
  to.output <- data.frame(test=cells.test.results$method, cells.statistic=cells.test.results$statistic, cells.p.value=cells.test.results$p.value, cells.significance, clones.statistic=clones.test.results$statistic, clones.p.value=clones.test.results$p.value, clones.significance)
  tmp.file.name <- paste0(reports.path, '/WilcoxonTestSevereVsMildResults.csv')
  write.csv(file=tmp.file.name, x=to.output)
  # cat('Test p-value:\n', test.results$p.value, '\n')

  # ---> Significance test between groups (only when only two groups are provided.)
  if(length(clusters.to.keep)==2){
    # @ Non-paired Wilcoxon Test.
    # Including each cell expansion value.
    group.1.sizes <- tmp.df[get(clusters.tag)==clusters.to.keep[1], clon.size.tag]
    group.2.sizes <- tmp.df[get(clusters.tag)==clusters.to.keep[2], clon.size.tag]
    cells.test.results <- wilcox.test(x=group.1.sizes, y=group.2.sizes, alternative='two.sided')
    # Including each unique clone expansion value.
    group.1.sizes <- unique(tmp.df[get(clusters.tag)==clusters.to.keep[1], .(clonotype.tag, clon.size.tag)])
    group.1.sizes <- group.1.sizes[, clon.size.tag]
    group.2.sizes <- unique(tmp.df[get(clusters.tag)==clusters.to.keep[2], .(clonotype.tag, clon.size.tag)])
    group.2.sizes <- group.2.sizes[, clon.size.tag]
    clones.test.results <- wilcox.test(x=group.1.sizes, y=group.2.sizes, alternative='two.sided')
    # Format and output.
    cells.significance <- ifelse(test=cells.test.results$p.value>0.05, yes=NA,
      no=ifelse(test=cells.test.results$p.value>0.01, yes='*',
      no=ifelse(test=cells.test.results$p.value>0.001, yes='**',
      no=ifelse(test=cells.test.results$p.value>0.0001, yes='***', no='****'))))
    clones.significance <- ifelse(test=clones.test.results$p.value>0.05, yes=NA,
      no=ifelse(test=clones.test.results$p.value>0.01, yes='*',
      no=ifelse(test=clones.test.results$p.value>0.001, yes='**',
      no=ifelse(test=clones.test.results$p.value>0.0001, yes='***', no='****'))))
    to.output <- data.frame(test=cells.test.results$method, cells.statistic=cells.test.results$statistic, cells.p.value=cells.test.results$p.value, cells.significance, clones.statistic=clones.test.results$statistic, clones.p.value=clones.test.results$p.value, clones.significance)
    tmp.file.name <- paste0(reports.path, '/WilcoxonTestCluster', clusters.to.keep[1], 'VsCluster', clusters.to.keep[2], '.csv')
    write.csv(file=tmp.file.name, x=to.output)
  }
  # Mock return.
  return(NA)
}


# 5 -------------------------------------------------------------------->
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

# 6 -------------------------------------------------------------------->
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
  to.output <- meta.data[virus.tag=='CV', .('Mean clone size'=mean(clon.size.tag, na.rm=TRUE)), by=.(Cluster=RNA_snn_res.0.6)]
  tmp.file.name <- paste0(tmp.path, '/MeanCloneSizePerClusterCVOnly.csv')
  fwrite(file=tmp.file.name, x=to.output)
  # All cells.
  to.output <- meta.data[, .('Mean clone size'=mean(clon.size.tag, na.rm=TRUE)), by=.(Cluster=RNA_snn_res.0.6)]
  tmp.file.name <- paste0(tmp.path, '/MeanCloneSizePerClusterAll.csv')
  fwrite(file=tmp.file.name, x=to.output)
}

# 6 -------------------------------------------------------------------->
# Name: Get clone abundance info by tags, donor x cluster.
# Description:
#
# Arguments ------------------------->
# @ this.seurat.obj - Seurat object of interest.
# @ reports.path - Reports path where tag-specific analysis folder should already exist.
# Function:

get.abundance.by.tags <- function(this.seurat.obj, reports.path, clusters.tag){
  # Get table.
  meta.data <- as.data.table(this.seurat.obj@meta.data)
  meta.data[, expansion.tag:=ifelse(test=clon.size.tag>1, yes='abundant', no='non.abundant')]
  meta.data <- meta.data[!is.na(orig.donor) & !is.na(expansion.tag)]
  tmp.data <- meta.data[, .(cells=.N), by=.(donor=orig.donor, hospital=orig.hospital, cluster=get(clusters.tag), expansion.tag)]
  tmp.data[, expansion.tag:=paste0(expansion.tag, '.size')]
  tmp.data <- spread(data=tmp.data, key=expansion.tag, value=cells, fill=0)
  tmp.data[, total:=(abundant.size + non.abundant.size)]
  tmp.data[, abundant.prop:=(abundant.size/total)]
  tmp.data[, non.abundant.prop:=(non.abundant.size/total)]
  # Output results.
  tmp.path <- paste0(reports.path, '/tag_specific_analysis')
  if(!dir.exists(tmp.path)) stop('Tag-specific analysis does not exist.\n')
  tmp.file.name <- paste0(tmp.path, '/ExpandedCellsPerDonorByCluster.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na='NA', quote=FALSE)
}

############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############
cat('### ----------------------- General Arguments ----------------------- ###\n')
# ---> For tag-specific analysis.
# @ Tags to depict specifically.
gen.tags.of.int <- 'orig.hospital'
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
# ---> Reports path
gen.reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels'

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
reports.path.six <- paste0(gen.reports.path, '/expansion_6hrs')
if(!dir.exists(reports.path.six)) dir.create(reports.path.six)
# --> vdj files.
gen.vdj.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/06-12-2020/aggr_vdj/COVID_GenEx_Batches-1-2-3_CD4_STIM_6/COVID_GenEx_Batches-1-2-3_CD4_STIM_6_07-24-2020'
# @ csv file depicting relationships between clones and cells. Aggregation of filtered_contig_annotations.csv files output by cellranger vdj from each independent sample.
cells.clons.info.file <- paste0(gen.vdj.path, '/filtered_contig_annotations_aggr.csv')
file.exists(cells.clons.info.file)
# @ csv file describing clone's specific info. Same, file result of the aggregation of a file output by cellranger vdj from each independent sample, in this case 'clonotypes.csv'
clons.info.file <- paste0(gen.vdj.path, '/clonotypes_aggr.csv')
file.exists(clons.info.file)
# @ Seurat object file processed for all basic workflow steps (in brief, cell normalization (LogNormalize), cell normalization, picking of variable features, clustering and dimensionality reduction for UMAP embeddings).
seurat.obj.six.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2-3_CD4_STIM_6/COVID_GenEx_Batches-1-2-3_CD4_STIM_6_07-23-2020_qc-std_var-25_pc_38_hto-all_stim_stim_virus_all_stimtime_six_ciros/seurat_objects/SeuratObj.RDS'
file.exists(seurat.obj.six.file)
# ---> Others
# @ Clusters tag.
clusters.tag.six <- 'RNA_snn_res.0.6'
# @ Tags of interest.
tags.of.int.six <- c(clusters.tag.six, 'orig.donor', gen.tags.of.int)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj.six <- readRDS(file=seurat.obj.six.file)
# Check that for the seurat object, we have dim. reductions methods applied.
if(!typeof(seurat.obj.six@reductions$umap)=="S4") stop('There is not either UMAP or tSNE applied to the input seurat object.')
# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

### ---------------------- Data preprocessing ----------------------- ###

# ---> Meta data entries.
seurat.obj.six@meta.data[seurat.obj.six@meta.data=='void'] <- NA

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
# ---> TCR annotations.
seurat.obj.six <- do.tcr.annotations(reports.path=reports.path.six, this.seurat.obj=seurat.obj.six, clons.info=clons.info, cells.clons.info=cells.clons.info)

# ---> Recaclulate clone size based on the cells captured in the gene expression data.
tmp.clons.data <- as.data.table(seurat.obj.six@meta.data[, c('clonotype.tag', 'clon.size.tag')])
tmp.clons.data <- tmp.clons.data[!is.na(clonotype.tag), .(clon.size.tag=.N), by=clonotype.tag]
tmp.clons.data <- as.data.frame(tmp.clons.data)
tmp.names <- tmp.clons.data$clonotype.tag
tmp.clons.data <- tmp.clons.data$clon.size.tag
names(tmp.clons.data) <- tmp.names
seurat.obj.six@meta.data[, 'clon.size.tag'] <-tmp.clons.data[seurat.obj.six@meta.data[, 'clonotype.tag']]

# ---> General TCR and tag-specific analysis.
do.tcr.analysis(reports.path=reports.path.six, this.seurat.obj=seurat.obj.six, tags.of.int=tags.of.int.six)

# ---> Donors by clusters.
get.abundance.by.tags(this.seurat.obj=seurat.obj.six, reports.path=reports.path.six, clusters.tag=clusters.tag.six)

############    -----------------------------------------    ############
############    ------------   Expansion    -------------    ############
############    ---------   06 hrs. data set    ---------    ############
############    -----------------------------------------    ############
cta('############    ------------   Expansion    -------------    ############\n')
cat('############    ---------   06 hrs. data set    ---------    ############\n')

### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# ---> File and paths definitions.
# All should be absolute paths.
# @ Reports path (where results should be deposited).
reports.path.comb <- paste0(gen.reports.path, '/expansion_0-6hrs')
if(!dir.exists(reports.path.comb)) dir.create(reports.path.comb)
# --> vdj files.
gen.vdj.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/06-12-2020/aggr_vdj/COVID_GenEx_Batches-1-2-3_CD4_COMBINED_0-6/COVID_GenEx_Batches-1-2-3_CD4_COMBINED_0-6_07-24-2020'
# @ csv file depicting relationships between clones and cells. Aggregation of filtered_contig_annotations.csv files output by cellranger vdj from each independent sample.
cells.clons.info.file <- paste0(gen.vdj.path, '/filtered_contig_annotations_aggr.csv')
file.exists(cells.clons.info.file)
# @ csv file describing clone's specific info. Same, file result of the aggregation of a file output by cellranger vdj from each independent sample, in this case 'clonotypes.csv'
clons.info.file <- paste0(gen.vdj.path, '/clonotypes_aggr.csv')
file.exists(clons.info.file)
# @ Seurat object file processed for all basic workflow steps (in brief, cell normalization (LogNormalize), cell normalization, picking of variable features, clustering and dimensionality reduction for UMAP embeddings).
seurat.obj.file.comb <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2-3_CD4_COMBINED_0-6/COVID_GenEx_Batches-1-2-3_CD4_COMBINED_0-6_07-23-2020_qc-std_var-25_pc_30_hto-all_stim_combined_virus_all_stimtime_zero-six_ciros/seurat_objects/SeuratObj.RDS'
# @ Tags of interest.
tags.of.int.comb <- c('RNA_snn_res.0.6', gen.tags.of.int)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj.comb <- readRDS(file=seurat.obj.file.comb)
# Check that for the seurat object, we have dim. reductions methods applied.
if(!typeof(seurat.obj.comb@reductions$umap)=="S4") stop('There is not either UMAP or tSNE applied to the input seurat object.')
# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

### ---------------------- Data preprocessing ----------------------- ###
# ---> UMAP coordinates.
# Flip UMAP embbedings by changing UMAP 2 dimension sign.
# seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_2'] <- -seurat.obj@reductions$umap@cell.embeddings[, 'UMAP_2']

# ---> Meta data entries.
seurat.obj.comb@meta.data[seurat.obj.comb@meta.data=='void'] <- NA

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
# ---> General workflow.
# ---> TCR annotations.
seurat.obj.comb <- do.tcr.annotations(reports.path=reports.path.comb, this.seurat.obj=seurat.obj.comb, clons.info=clons.info, cells.clons.info=cells.clons.info)

# ---> Recaclulate clone size based on the cells captured in the gene expression data.
tmp.clons.data <- as.data.table(seurat.obj.comb@meta.data[, c('clonotype.tag', 'clon.size.tag')])
tmp.clons.data <- tmp.clons.data[!is.na(clonotype.tag), .(clon.size.tag=.N), by=clonotype.tag]
tmp.clons.data <- as.data.frame(tmp.clons.data)
tmp.names <- tmp.clons.data$clonotype.tag
tmp.clons.data <- tmp.clons.data$clon.size.tag
names(tmp.clons.data) <- tmp.names
seurat.obj.comb@meta.data[, 'clon.size.tag'] <-tmp.clons.data[seurat.obj.comb@meta.data[, 'clonotype.tag']]

# ---> General TCR and tag-specific analysis.
do.tcr.analysis(reports.path=reports.path.comb, this.seurat.obj=seurat.obj.comb, tags.of.int=tags.of.int.comb)

# ---> Clone size across clusters.
# @ For all clusters.
# vln.path <- paste0(reports.path.all, '/dim_reduction/violin_plots/all_clusters')
# if(!dir.exists(vln.path)) dir.create(vln.path, recursive=TRUE)
# # Only for COVID-reactive cells.
# tmp.path <- paste0(vln.path, '/covid_only')
# if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive=TRUE)
# get.clusts.expansion(this.seurat.obj=seurat.obj.all, clusters.tag=clusters.tag.all, reports.path=tmp.path, covid.only=TRUE)
# # For all virus-reactive cells.
# tmp.path <- paste0(vln.path, '/all_viruses')
# if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive=TRUE)
# get.clusts.expansion(this.seurat.obj=seurat.obj.all, clusters.tag=clusters.tag.all, reports.path=tmp.path, covid.only=FALSE)
#
# # @ For a set of clusters.
# vln.path <- paste0(reports.path.all, '/dim_reduction/violin_plots/clusters_subset')
# if(!dir.exists(vln.path)) dir.create(vln.path, recursive=TRUE)
# # Only for COVID-reactive cells.
# tmp.path <- paste0(vln.path, '/covid_only')
# if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive=TRUE)
# get.clusts.expansion(this.seurat.obj=seurat.obj.all, clusters.tag=clusters.tag, reports.path=tmp.path, covid.only=TRUE, clusters.to.keep=c(0, 4, 5, 9))
# # For all virus-reactive cells.
# tmp.path <- paste0(vln.path, '/all_viruses')
# if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive=TRUE)
# get.clusts.expansion(this.seurat.obj=seurat.obj.all, clusters.tag=clusters.tag, reports.path=tmp.path, covid.only=FALSE, clusters.to.keep=0:3)
#
# # @ For cluster 0 alone.
# vln.path <- paste0(reports.path.all, '/dim_reduction/violin_plots/clusters_subset_0')
# if(!dir.exists(vln.path)) dir.create(vln.path, recursive=TRUE)
# # Only for COVID-reactive cells.
# tmp.path <- paste0(vln.path, '/covid_only')
# if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive=TRUE)
# get.clusts.expansion(this.seurat.obj=seurat.obj.all, clusters.tag=clusters.tag, reports.path=tmp.path, covid.only=TRUE, clusters.to.keep=0, adjust.val=0.8, width.val=7)
# # For all virus-reactive cells.
# tmp.path <- paste0(vln.path, '/all_viruses')
# if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive=TRUE)
# get.clusts.expansion(this.seurat.obj=seurat.obj.all, clusters.tag=clusters.tag, reports.path=tmp.path, covid.only=FALSE, clusters.to.keep=0, adjust.val=0.8, width.val=7)


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
reports.path.tf <- paste0(gen.reports.path, '/expansion_24hrs')
if(!dir.exists(reports.path.tf)) dir.create(reports.path.tf)
# --> vdj files.
gen.vdj.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/06-12-2020/aggr_vdj/COVID_GenEx_Batches-1-2-3_CD4_STIM_24/COVID_GenEx_Batches-1-2-3_CD4_STIM_24_07-24-2020'
# @ csv file depicting relationships between clones and cells. Aggregation of filtered_contig_annotations.csv files output by cellranger vdj from each independent sample.
cells.clons.info.file <- paste0(gen.vdj.path, '/filtered_contig_annotations_aggr.csv')
file.exists(cells.clons.info.file)
# @ csv file describing clone's specific info. Same, file result of the aggregation of a file output by cellranger vdj from each independent sample, in this case 'clonotypes.csv'
clons.info.file <- paste0(gen.vdj.path, '/clonotypes_aggr.csv')
file.exists(clons.info.file)
# @ Seurat object file processed for all basic workflow steps (in brief, cell normalization (LogNormalize), cell normalization, picking of variable features, clustering and dimensionality reduction for UMAP embeddings).
seurat.obj.file.tf <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-2-3_CD4_STIM_24/COVID_GenEx_Batches-1-2-3_CD4_STIM_24_07-23-2020_qc-std_var-25_pc_16_hto-all_stim_stim_virus_all_stimtime_twentyfour_ciros/seurat_objects/SeuratObj.RDS'
# @ Tags of interest.
clusters.tag.tf <- 'RNA_snn_res.0.2'
tags.of.int.tf <- c('RNA_snn_res.0.2', gen.tags.of.int)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# Seurat object ------------------------------------------------------->
seurat.obj.tf <- readRDS(file=seurat.obj.file.tf)
# Check that for the seurat object, we have dim. reductions methods applied.
if(!typeof(seurat.obj.tf@reductions$umap)=="S4") stop('There is not either UMAP or tSNE applied to the input seurat object.')
# TCR data files ------------------------------------------------------>
# ---> Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# ---> Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

### ---------------------- Data preprocessing ----------------------- ###
# ---> Meta data entries.
seurat.obj.tf@meta.data[seurat.obj.tf@meta.data=='void'] <- NA

# ---> CLusters tag.
# Clusters tag. 0 nad 6 will be merged into a single cluster since they're both Tregs.
to.replace <- 1:length(Cells(seurat.obj.tf))
to.replace <- to.replace[seurat.obj.tf@meta.data[, clusters.tag.tf]=='6']
seurat.obj.tf@meta.data[, clusters.tag.tf] <- replace(x=seurat.obj.tf@meta.data[, clusters.tag.tf], list=to.replace, values='0')

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
# ---> General workflow.
# ---> TCR annotations.
seurat.obj.tf <- do.tcr.annotations(reports.path=reports.path.tf, this.seurat.obj=seurat.obj.tf, clons.info=clons.info, cells.clons.info=cells.clons.info)

# ---> Recaclulate clone size based on the cells captured in the gene expression data.
tmp.clons.data <- as.data.table(seurat.obj.tf@meta.data[, c('clonotype.tag', 'clon.size.tag')])
tmp.clons.data <- tmp.clons.data[!is.na(clonotype.tag), .(clon.size.tag=.N), by=clonotype.tag]
tmp.clons.data <- as.data.frame(tmp.clons.data)
tmp.names <- tmp.clons.data$clonotype.tag
tmp.clons.data <- tmp.clons.data$clon.size.tag
names(tmp.clons.data) <- tmp.names
seurat.obj.tf@meta.data[, 'clon.size.tag'] <-tmp.clons.data[seurat.obj.tf@meta.data[, 'clonotype.tag']]

# ---> General TCR and tag-specific analysis.
do.tcr.analysis(reports.path=reports.path.tf, this.seurat.obj=seurat.obj.tf, tags.of.int=tags.of.int.tf)

# ---> Donors by clusters.
get.abundance.by.tags(this.seurat.obj=seurat.obj.tf, reports.path=reports.path.tf, clusters.tag=clusters.tag.tf)

############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ----------   Clones' sharing    ---------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')
cat('############    ----------   Clones\' sharing    ---------    ############')
cat('############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############\n')

### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# ---> vdj data.
aggr.vdj.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/06-12-2020/aggr_vdj/COVID_GenEx_Batches-1-2-3_CD4_COMBINED_0-6-24/COVID_GenEx_Batches-1-2-3_CD4_COMBINED_0-6-24_07-24-2020'
cells.clons.info.file <- paste0(aggr.vdj.path, '/filtered_contig_annotations_aggr.csv')
clons.info.file <- paste0(aggr.vdj.path, '/clonotypes_aggr.csv')
file.exists(cells.clons.info.file)
file.exists(clons.info.file)
# ---> Seurat objects.
seurat.obj.comb.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/expansion_0-6hrs/seurat_objects/SeuratObj.RDS'
seurat.obj.six.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/expansion_6hrs/seurat_objects/SeuratObj.RDS'
seurat.obj.tf.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/expansion_24hrs/seurat_objects/SeuratObj.RDS'
file.exists(seurat.obj.comb.file)
file.exists(seurat.obj.six.file)
file.exists(seurat.obj.tf.file)
# ---> Translation dictionary.
# A file for a dictionary that allows us to translate the barcode suffixes between the 6hrs.-only aggregation and the 6-0 hrs. combination aggregation.
sffxs.dict.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/general_data/SufixxesDictionaryBetween0-6HrsAnd6Hrs.csv'
# ---> Clusters tags.
# Clusters column.
clusters.tag.six <- 'RNA_snn_res.0.6'
clusters.tag.tf <- 'RNA_snn_res.0.2'
clusters.tag.comb <- 'RNA_snn_res.0.6'
# Clusters name.
clusters.names.six <- 0:20
names(clusters.names.six) <- 0:(length(clusters.names.six)-1)
clusters.names.tf <- LETTERS
names(clusters.names.tf) <- 0:(length(clusters.names.tf)-1)
clusters.names.comb <- letters
names(clusters.names.comb) <- 0:(length(clusters.names.comb)-1)
# ---> Program-specific parameters.
cells.clons.file.suffix <- 21:26
names(cells.clons.file.suffix) <- 1:6
exp.thold <- 2
freq.tholds <- c(1, 2, 10)
# donor.id.tag <- 'donor.id.tag'
# Clusters to look at.
# six.cluster.sets <- paste0(paste0('cluster.', 0:12), '.six')
# tf.cluster.sets <- paste0(paste0('cluster.', 0:5), '.twentyfour')
cols.to.keep <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'orig.cell_type', 'orig.virus2', 'orig.peptide', 'orig.stim_time', 'orig.donor', 'ht_severity', 'orig.sex', 'orig.hospital', 'orig.severity', 'clusters.tag')
# ---> Path definitons.
sharing.reports.path <- paste0(gen.reports.path, '/clones_sharing')
if(!dir.exists(sharing.reports.path)) dir.create(sharing.reports.path)
#dir.create(reports.path)

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# ---> Seurat objects.
seurat.obj.comb <- readRDS(file=seurat.obj.comb.file)
seurat.obj.six <- readRDS(file=seurat.obj.six.file)
seurat.obj.tf <- readRDS(file=seurat.obj.tf.file)
# ---> Suffixes dictionary.
sffxs.dict <- read.csv(file=sffxs.dict.file, stringsAsFactors=FALSE)
# ---> TCR data files
# Clonotypes info
clons.info <- read.csv(file=clons.info.file, stringsAsFactors=FALSE)
cat('Both, TCR data and seurat object read to R objects. Check for warnings or errors if any.\n\n')
# Cells-clonotypes relationships info.
cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)

cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
cat('### ---------------------- Data Preprocessing ---------------------- ###\n')

data.path <- paste0(sharing.reports.path, '/source_data')
if(!dir.exists(data.path)) dir.create(data.path)

# ---> Seurat object.
tf.cell.names <- Cells(seurat.obj.tf)
for(suffix in names(cells.clons.file.suffix)) tf.cell.names <- str_replace(string=tf.cell.names, pattern=paste0('-', suffix), replacement=paste0('-', cells.clons.file.suffix[suffix]))
seurat.obj.tf <- RenameCells(object=seurat.obj.tf, new.names=tf.cell.names)

# ---> Clusters tag.
# 24 hrs.
new.levels <- mixedsort(unique(as.character(seurat.obj.tf@meta.data[, clusters.tag.tf])))
seurat.obj.tf@meta.data[, clusters.tag.tf] <-factor(x=as.character(seurat.obj.tf@meta.data[, clusters.tag.tf]), levels=new.levels)
seurat.obj.tf@meta.data[, 'clusters.tag'] <- seurat.obj.tf@meta.data[, clusters.tag.tf]
seurat.obj.tf@meta.data[, 'clusters.tag'] <- clusters.names.tf[seurat.obj.tf@meta.data[, 'clusters.tag']]
# 6 hrs.
new.levels <- mixedsort(unique(as.character(seurat.obj.six@meta.data[, clusters.tag.six])))
seurat.obj.six@meta.data[, clusters.tag.six] <-factor(x=as.character(seurat.obj.six@meta.data[, clusters.tag.six]), levels=new.levels)
seurat.obj.six@meta.data[, 'clusters.tag'] <- seurat.obj.six@meta.data[, clusters.tag.six]
seurat.obj.six@meta.data[, 'clusters.tag'] <- clusters.names.six[seurat.obj.six@meta.data[, 'clusters.tag']]
# 6-0 hrs. combinations.
new.levels <- mixedsort(unique(as.character(seurat.obj.comb@meta.data[, clusters.tag.comb])))
seurat.obj.comb@meta.data[, clusters.tag.comb] <-factor(x=as.character(seurat.obj.comb@meta.data[, clusters.tag.comb]), levels=new.levels)
seurat.obj.comb@meta.data[, 'clusters.tag'] <- seurat.obj.comb@meta.data[, clusters.tag.comb]
seurat.obj.comb@meta.data[, 'clusters.tag'] <- clusters.names.comb[seurat.obj.comb@meta.data[, 'clusters.tag']]

# ---> Keep columns of interest.
if(!all(all(cols.to.keep %in% colnames(seurat.obj.comb@meta.data)) & all(cols.to.keep %in% colnames(seurat.obj.six@meta.data))
 & all(cols.to.keep %in% colnames(seurat.obj.tf@meta.data)))) stop(paste0('Not all columns of interest defined appropriately.\n'))
seurat.obj.six@meta.data <- seurat.obj.six@meta.data[, cols.to.keep]
seurat.obj.tf@meta.data <- seurat.obj.tf@meta.data[, cols.to.keep]
seurat.obj.comb@meta.data <- seurat.obj.comb@meta.data[, cols.to.keep]

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

# ---> Combine TCR data from all data sets.
# Take the overall metadata
meta.data <- bind_rows(seurat.obj.comb@meta.data, seurat.obj.tf@meta.data, .id='data.set')
meta.data[, 'data.set'] <- ifelse(test=meta.data[, 'data.set']==1, yes='combination', no='twentyfour')
rownames(meta.data) <- c(rownames(seurat.obj.comb@meta.data), rownames(seurat.obj.tf@meta.data))
meta.data$barcode <- rownames(meta.data)
# Split clusters tags' appropriately.
meta.data <- as.data.table(meta.data)
meta.data[data.set=='twentyfour', clusters.tf.tag:=clusters.tag]
meta.data[data.set=='combination', clusters.comb.tag:=clusters.tag]
# @@ Add the six hours-only aggregation clusters.
clusters.info.six <- data.table(orig.cells=Cells(seurat.obj.six), clusters.tag=seurat.obj.six@meta.data[, 'clusters.tag'])
for(idx in 1:nrow(sffxs.dict)){
  six.sffx <- paste0("-", sffxs.dict[idx, 'suffix.six'], "$")
  comb.sffx <- paste0("-", sffxs.dict[idx, 'suffix.comb'])
  clusters.info.six[str_detect(string=orig.cells, pattern=six.sffx), alt.cells:=str_replace(string=orig.cells, pattern=six.sffx, replacement=comb.sffx)]
}
# Missing cells between aggregations.
missing.cells.from.six <- clusters.info.six[!alt.cells%chin%meta.data[, barcode], orig.cells]
missing.cells.info <- as.data.table(seurat.obj.six@meta.data[missing.cells.from.six, ])
# missing.cells.from.comb <- clusters.info.six[!alt.cells%chin%meta.data[, barcode], alt.cells]
# Making sure they're lost in the combination aggregation due to be low quality.
qc.tholds <- list('nCount_RNA'=c(1500+50, 20000-50), 'nFeature_RNA'=c(800+50, 4400-50), 'percent.mt'=c(0, 10-1))
missing.cells.evals <- lapply(X=names(qc.tholds), FUN=function(criteria){
  low.thold <- qc.tholds[[criteria]][1]
  up.thold <- qc.tholds[[criteria]][2]
  missing.cells.info[, get(criteria)>low.thold & get(criteria)<up.thold]
})
missing.cells.evals <- rlist::list.cbind(missing.cells.evals)
sum(apply(X=missing.cells.evals, MARGIN=1, FUN=all))
# Answer: All of these cells were removed from the 0-6 hrs. aggregation most likely due to their quality.
# Proceed to add the 6 hrs. info.
clusters.info.six <- clusters.info.six[, .(barcode=alt.cells, six.hrs.barcode=orig.cells, clusters.six.tag=clusters.tag)]
meta.data <- merge(x=meta.data, y=clusters.info.six, by='barcode', all.x=TRUE)

# ---> Merge GenEx and TCR data.
# Combine clones data with other meta data.
cells.clons.info <- as.data.table(merge(x=cells.clons.info, y=meta.data, by='barcode', all.x=TRUE, all.y=FALSE))

# ---> Clones info table.
# Filter out relationships with no clustering info available (no genex info).
cells.clons.info.tmp <- cells.clons.info[!is.na(clusters.six.tag) | !is.na(clusters.tf.tag) | !is.na(clusters.comb.tag)]

# ---> Save info.
tmp.file.name <- paste0(data.path, '/CellsClonsInfo.csv')
fwrite(file=tmp.file.name, x=cells.clons.info.tmp)

cat('\n\n')
############    ----------------------------------------    ############
############    ---------   Clones' sharing    ---------    ############
############    -------   for the 6 hrs. set    --------    ############
############    ----------------------------------------    ############
cat('############    ---------   Clones\' sharing    ---------    ############\n')
cat('############    -------   for the 6 hrs. set    --------    ############\n')

# Define reports path.
six.sharing.path <- paste0(sharing.reports.path, '/within_6hrs')
if(!dir.exists(six.sharing.path)) dir.create(six.sharing.path)

# Define the clusters for the data set.
set.names <- paste0('cluster.', cells.clons.info[!is.na(clusters.six.tag), unique(clusters.six.tag)])

############    ---------------------------------------    ############
############    ----------   UpSetR Plots    ----------    ############
############    ---------------------------------------    ############
upset.path <- paste0(six.sharing.path, '/upset_plots')
if(!dir.exists(upset.path)) dir.create(upset.path)

# ---> Load data.
tmp.file.name <- paste0(data.path, '/CellsClonsInfo.csv')
cells.clons.info <- fread(file=tmp.file.name, na.strings='')

############    ---------------------------------------    ############
############    ----------   UpSetR Plots    ----------    ############
############    --------   Based on Clones    ---------    ############
############    ---------------------------------------    ############

upset.clones.path <- paste0(upset.path, '/based_on_clones')
if(!dir.exists(upset.clones.path)) dir.create(upset.clones.path)

for(freq.thold in freq.tholds){
  # ---> Processing.
  clons.table <- cells.clons.info[orig.stim_time==6 & !is.na(clusters.six.tag) & orig.virus2=='CV', .SD[, .(cell.no=uniqueN(six.hrs.barcode)), by=.(cluster=clusters.six.tag)], by=.(clone.id=raw_clonotype_id)]
  clons.table <- spread(data=clons.table, key=cluster, value=cell.no, fill=0, sep='.')

  for(tmp.set in set.names){
    clons.table[, eval(tmp.set):=ifelse(test=get(tmp.set)>=freq.thold, yes=1, no=0)]
  }

  # ---> UpSetR Plot.
  # Normal version.
  tmp.file.name <- paste0(upset.clones.path, '/SharingAmongDataSetsWithThold_', freq.thold, '_UpSetR.pdf')
  pdf(file=tmp.file.name, 12)
  print(upset(data=clons.table,
    sets=set.names,
    nintersects=NA,
    query.legend="top",
    order.by=c("freq", "degree"),
    # keep.order=TRUE,
    mainbar.y.label='Number of intersected clones',
    sets.x.label=paste0("Total clones (freq. > ", freq.thold, ")"),
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    # scale.intersections=tmp.scale
  ))
  dev.off()
  # Blank version.
  tmp.file.name <- paste0(upset.clones.path, '/SharingAmongDataSetsWithThold_', freq.thold, '_UpSetR_Blank.pdf')
  pdf(file=tmp.file.name, 12)
  print(upset(data=clons.table,
    sets=set.names,
    nintersects=NA,
    query.legend="top",
    order.by=c("freq", "degree"),
    # keep.order=TRUE,
    mainbar.y.label='Number of intersected clones',
    sets.x.label=paste0("Total clones (freq. > ", freq.thold, ")"),
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    text.scale=c(0, 0, 0, 0, 0, 0)
  ))
  dev.off()
}

############    ---------------------------------------    ############
############    ----------   UpSetR Plots    ----------    ############
############    ---------   Based on Cells    ---------    ############
############    ---------------------------------------    ############

upset.cells.path <- paste0(upset.path, '/based_on_cells')
if(!dir.exists(upset.cells.path)) dir.create(upset.cells.path)

for(freq.thold in freq.tholds){
  # ---> Processing.
  clons.table <- cells.clons.info[orig.stim_time==6 & !is.na(clusters.six.tag) & orig.virus2=='CV', .SD[, .(cell.no=uniqueN(six.hrs.barcode)), by=.(cluster=clusters.six.tag)], by=.(clone.id=raw_clonotype_id)]
  clons.table <- spread(data=clons.table, key=cluster, value=cell.no, fill=0, sep='.')

  # ---> Create clone-basis structure.
  # Remove clones not seen in the 6H data set with a clone size larger than 1 for at least one of the clusters.
  clons.table <- as.data.frame(clons.table, stringsAsFactors=FALSE)
  clons.to.keep <- sapply(X=1:nrow(clons.table), FUN=function(idx){
    clone.info <- clons.table[idx, set.names]
    to.eval <- sum(clone.info>=freq.thold)
    return(to.eval>0)
  })
  clons.table <- clons.table[clons.to.keep, ]
  # Then, get cell-basis structure of the data. In brief, a clone row is repeated as many times as its clone size.
  clons.table$clone.size <- rowSums(clons.table[, set.names])
  clones.info.cb <- mclapply(X=1:nrow(clons.table), FUN=function(clone.row.id){
    clones.clusts.info <- clons.table[clone.row.id, set.names]
    times.to.rep <- clons.table[clone.row.id, 'clone.size']
    clones.clusts.info <- bind_rows(replicate(times.to.rep, clones.clusts.info, simplify = FALSE))
    return(clones.clusts.info)
  })
  clones.info.cb <- rbindlist(clones.info.cb)
  # Proof checking we got it right.
  # Total clone size should equal the total amount of rows we got in the processed data.
  sum(rowSums(clons.table[, set.names]))
  nrow(clones.info.cb)

  for(tmp.set in set.names){
    clones.info.cb[, eval(tmp.set):=ifelse(test=get(tmp.set)>=freq.thold, yes=1, no=0)]
  }

  # ---> UpSetR Plot.
  # Normal version.
  tmp.file.name <- paste0(upset.cells.path, '/SharingAmongDataSetsWithThold_', freq.thold, '_UpSetR.pdf')
  pdf(file=tmp.file.name, width=12)
  print(upset(data=clones.info.cb,
    sets=set.names,
    nintersects=NA,
    query.legend="top",
    order.by=c("freq", "degree"),
    # keep.order=TRUE,
    mainbar.y.label='Number of intersected cells',
    sets.x.label=paste0("Nonsense"),
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    # text.scale=c(0, 0, 0, 0, 0, 0)
  ))
  dev.off()
  # Blank version.
  tmp.file.name <- paste0(upset.cells.path, '/SharingAmongDataSetsWithThold_', freq.thold, '_UpSetR_Blank.pdf')
  pdf(file=tmp.file.name, width=12)
  print(upset(data=clones.info.cb,
    sets=set.names,
    nintersects=NA,
    query.legend="top",
    order.by=c("freq", "degree"),
    # keep.order=TRUE,
    mainbar.y.label='Number of intersected cells',
    sets.x.label=paste0("Nonsense"),
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    text.scale=c(0, 0, 0, 0, 0, 0)
  ))
  dev.off()
}

############    -----------------------------------------    ############
### ----------------------------- Cites ----------------------------- ###
############    -----------------------------------------    ############
# ---> UpSetR
# Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R Package for the Visualization of Intersecting Sets and their Properties doi: https://doi.org/10.1093/bioinformatics/btx364
