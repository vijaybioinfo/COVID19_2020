############    ----------   Clones' sharing    ---------    ############
############    ----------   among data sets    ---------    ############
cat('############    ----------   Clones\' sharing    ---------    ############')
cat('############    ----------   among data sets    ---------    ############\n')

### --------------------------- Libraries --------------------------- ###
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
source('/home/vfajardo/scripts/functions/R_handy_functions.R')


### ----------------------- General Arguments ----------------------- ###
# ---> Seurat objects.
seurat.obj.comb.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/expansion_0-6hrs/seurat_objects/SeuratObj.RDS'
seurat.obj.six.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/expansion_6hrs/seurat_objects/SeuratObj.RDS'
seurat.obj.tf.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/expansion_24hrs/seurat_objects/SeuratObj.RDS'
file.exists(seurat.obj.comb.file)
file.exists(seurat.obj.six.file)
file.exists(seurat.obj.tf.file)
# ----> Clonotypes data.
cells.clons.info.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/clones_sharing/source_data/CellsClonsInfo.csv'
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
# ---> Reports.
gen.reports.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/batches_1-2-3/CD4s_paper/final_figures/figure_3/figure_3_panels/clones_sharing'
reports.path <- paste0(gen.reports.path, '/among_datasets')
if(!dir.exists(reports.path)) dir.create(reports.path)
# ---> Program-specific parameters.
cols.to.keep <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'orig.cell_type', 'orig.virus2', 'orig.peptide', 'orig.stim_time', 'orig.donor', 'ht_severity', 'orig.sex', 'orig.hospital', 'orig.severity', 'clusters.tag', 'clon.size.tag')
freq.tholds <- c(1, 2, 10)
time.points <- c('0', '6', '24')
set.cols <- c('#8fbc8f', '#6495ed', '#cd5c5c')
set.names <- paste0('time.', time.points)
cells.clons.file.suffix <- 21:26
names(cells.clons.file.suffix) <- 1:6
tmp.cols <- binary.col.scale
names(tmp.cols) <- c('S', 'NS')
blank.complement <- theme(, line=element_blank(), text=element_blank(), legend.position='none', legend.background=element_blank())

cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# ---> Seurat objects.
seurat.obj.comb <- readRDS(file=seurat.obj.comb.file)
seurat.obj.six <- readRDS(file=seurat.obj.six.file)
seurat.obj.tf <- readRDS(file=seurat.obj.tf.file)
# ---> Cells-clonotypes info.
cells.clons.info <- fread(file=cells.clons.info.file, na.strings='')

cat('\n\n')
### ---------------------- Data Preprocessing ---------------------- ###
cat('### ---------------------- Data Preprocessing ---------------------- ###\n')
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

# ---> Set-specific barcodes.
# @ 6 hrs.
# Already defined.
# @ 24 hrs.
cells.clons.info[orig.stim_time==24, tf.hrs.barcode:=barcode]
for(sffx in names(cells.clons.file.suffix)){
  tf.sffx <- paste0("-", sffx)
  comb.sffx <- paste0("-", cells.clons.file.suffix[sffx], "$")
  cells.clons.info[orig.stim_time==24 & str_detect(string=barcode, pattern=comb.sffx), tf.hrs.barcode:=str_replace(string=barcode, pattern=comb.sffx, replacement=tf.sffx)]
}
# @ Combined data sets.
cells.clons.info[data.set=="combination" & orig.stim_time!=24, comb.hrs.barcode:=barcode]

# ---> Keep columns of interest.
if(!all(all(cols.to.keep %in% colnames(seurat.obj.comb@meta.data)) & all(cols.to.keep %in% colnames(seurat.obj.six@meta.data))
 & all(cols.to.keep %in% colnames(seurat.obj.tf@meta.data)))) stop(paste0('Not all columns of interest defined appropriately.\n'))
seurat.obj.six@meta.data <- seurat.obj.six@meta.data[, cols.to.keep]
seurat.obj.tf@meta.data <- seurat.obj.tf@meta.data[, cols.to.keep]
seurat.obj.comb@meta.data <- seurat.obj.comb@meta.data[, cols.to.keep]

# ---> 0-6 shared donors and subsample.
# Donors seen in six h.
donors.six <- as.data.table(seurat.obj.six@meta.data)
donors.six[, barcode:=Cells(seurat.obj.six)]
donors.six <- donors.six[!is.na(orig.donor) & !is.na(clon.size.tag), .(orig.donor, barcode)]
donors.cells.six <- donors.six[, .(cells=.N), by=.(donor=orig.donor)]
 # Donors seen in 0 h.
donors.zero <- as.data.table(seurat.obj.comb@meta.data)
donors.zero[, barcode:=Cells(seurat.obj.comb)]
donors.zero <- donors.zero[!is.na(orig.donor) & !is.na(clon.size.tag) & orig.stim_time==0, .(orig.donor, barcode)]
donors.cells.zero <- donors.zero[, .(cells=.N), by=.(donor=orig.donor)]
# Shared donors.
shared.donors <- merge(x=donors.cells.six, y=donors.cells.zero, by='donor', suffixes=c('.six', '.zero'))
shared.donors$cells.min <- apply(X=shared.donors[, .(cells.six, cells.zero)], MARGIN=1, FUN=min)
# Get whole sample per data set considering the minimum of cells per donor.
six.sample <- lapply(X=shared.donors$donor, FUN=function(tmp.donor){
  sample.size <- shared.donors[donor==tmp.donor, cells.min]
  universe <- donors.six[orig.donor==tmp.donor, barcode]
  tmp.sample <- sample(x=universe, size=sample.size, replace=FALSE)
  return(tmp.sample)
})
six.sample <- unlist(six.sample)
zero.sample <- lapply(X=shared.donors$donor, FUN=function(tmp.donor){
  sample.size <- shared.donors[donor==tmp.donor, cells.min]
  universe <- donors.zero[orig.donor==tmp.donor, barcode]
  tmp.sample <- sample(x=universe, size=sample.size, replace=FALSE)
  return(tmp.sample)
})
zero.sample <- unlist(zero.sample)

############    -----------------------------------------    ############
############    ----------   Clones' sharing    ---------    ############
############    -----------------------------------------    ############
cat('############    ----------   Clones\' sharing    ---------    ############')

############    -----------------------------------------    ############
############    ----------   Venn Diagrams    -----------    ############
############    -----------------------------------------    ############
diagrams.path <- paste0(reports.path, '/venn_diagrams')
if(!dir.exists(diagrams.path)) dir.create(diagrams.path)

# Raw data will be picked so that Venn diagrams can be created with the jvenn tool.

# ---> Based on clones
# Get data.
diagrams.info <- cells.clons.info[orig.virus2=='CV', .(time=orig.stim_time, cell.no=length(unique(barcode))), by=.(clone.id=raw_clonotype_id)]
diagrams.info <- unique(diagrams.info)
diagrams.info <- spread(data=diagrams.info, key=time, value=cell.no, fill=NA, sep='.')
venn.input <- lapply(X=set.names, FUN=function(tmp.set){
  to.output <- diagrams.info[!is.na(get(tmp.set)), clone.id]
  return(to.output)
})
names(venn.input) <- paste0(time.points, ' hrs.')
# Output Venn Diagram for all data sets.
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols, col=set.cols)
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_Blank.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols, col=set.cols, cex=0, cat.cex=0)
# Output Venn Diagram for 0 and 6 hrs. only.
venn.input <- venn.input[names(venn.input)!='24 hrs.']
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_0and6hrs.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2])
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_0and6hrs_Blank.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2], cex=0, cat.cex=0)
# Output data.
tmp.file.name <- paste0(diagrams.path, '/DataForSharingByVennDiagrams_BasedOnClones.csv')
fwrite(file=tmp.file.name, x=diagrams.info, na='NA')

# ---> Based on clones only for cells shared between 0 and 6 hrs, all cells.
# Get data.
diagrams.info <- cells.clons.info[orig.donor %chin% shared.donors$donor & orig.virus2=='CV', .(time=orig.stim_time, cell.no=length(unique(barcode))), by=.(clone.id=raw_clonotype_id)]
diagrams.info <- unique(diagrams.info)
diagrams.info <- spread(data=diagrams.info, key=time, value=cell.no, fill=NA, sep='.')
venn.input <- lapply(X=set.names, FUN=function(tmp.set){
  to.output <- diagrams.info[!is.na(get(tmp.set)), clone.id]
  return(to.output)
})
names(venn.input) <- paste0(time.points, ' hrs.')
# Output Venn Diagram for all data sets.
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_MatchedDonors_AllCells.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols, col=set.cols)
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_MatchedDonors_AllCells_Blank.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols, col=set.cols, cex=0, cat.cex=0)
# Output Venn Diagram for 0 and 6 hrs. only.
venn.input <- venn.input[names(venn.input)!='24 hrs.']
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_0and6hrs_MatchedDonors_AllCells.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2])
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_0and6hrs_MatchedDonors_AllCells_Blank.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2], cex=0, cat.cex=0)
# Output data.
tmp.file.name <- paste0(diagrams.path, '/DataForSharingByVennDiagrams_BasedOnClones_MatchedDonors_AllCells.csv')
fwrite(file=tmp.file.name, x=diagrams.info, na='NA')

# ---> Based on clones only for cells shared between 0 and 6 hrs, subsample.
# Get data.
tmp.1 <- cells.clons.info[orig.stim_time==0 & comb.hrs.barcode %chin% zero.sample, .(time=orig.stim_time, cell.no=length(unique(barcode))), by=.(clone.id=raw_clonotype_id)]
tmp.2 <- cells.clons.info[orig.stim_time==6 & six.hrs.barcode %chin% six.sample, .(time=orig.stim_time, cell.no=length(unique(barcode))), by=.(clone.id=raw_clonotype_id)]
# diagrams.info <- cells.clons.info[orig.donor %chin% shared.donors$donor & orig.virus2=='CV', .(time=orig.stim_time, cell.no=length(unique(barcode))), by=.(clone.id=raw_clonotype_id)]
diagrams.info <- rbind(tmp.1, tmp.2)
diagrams.info <- unique(diagrams.info)
diagrams.info <- spread(data=diagrams.info, key=time, value=cell.no, fill=NA, sep='.')
venn.input <- lapply(X=c('time.0', 'time.6'), FUN=function(tmp.set){
  to.output <- diagrams.info[!is.na(get(tmp.set)), clone.id]
  return(to.output)
})
names(venn.input) <- paste0(c('time.0', 'time.6'), ' hrs.')
# Output Venn Diagram for 0 and 6 hrs. only.
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_0and6hrs_MatchedDonors_Subsample.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2])
tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnClones_0and6hrs_MatchedDonors_Subsample_Blank.png')
venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2], cex=0, cat.cex=0)
# Output data.
tmp.file.name <- paste0(diagrams.path, '/DataForSharingByVennDiagrams_BasedOnClones_MatchedDonors_Subsample.csv')
fwrite(file=tmp.file.name, x=diagrams.info, na='NA')

# # # ---> Based on cells for all sets.
# # # Get data.
# # diagrams.info <- as.data.frame(diagrams.info, stringsAsFactors=FALSE)
# # diagrams.info$clone.size <- rowSums(diagrams.info[, set.names], na.rm=TRUE)
# # diagrams.info <- mclapply(X=1:nrow(diagrams.info), FUN=function(clone.row.id){
# #   clones.clusts.info <- diagrams.info[clone.row.id, set.names]
# #   times.to.rep <- diagrams.info[clone.row.id, 'clone.size']
# #   clones.clusts.info <- bind_rows(replicate(times.to.rep, clones.clusts.info, simplify = FALSE))
# #   clones.clusts.info$clone.id <- paste(diagrams.info[clone.row.id, 'clone.id'], 1:times.to.rep, sep='.')
# #   return(clones.clusts.info)
# # })
# # diagrams.info <- rbindlist(diagrams.info)
# # venn.input <- lapply(X=set.names, FUN=function(tmp.set){
# #   to.output <- diagrams.info[!is.na(get(tmp.set)), clone.id]
# #   return(to.output)
# # })
# # names(venn.input) <- paste0(time.points, ' hrs.')
# # # Output Venn Diagram.
# # tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnCells.png')
# # venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols, col=set.cols)
# # tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnCells_Blank.png')
# # venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols, col=set.cols, cex=0, cat.cex=0)
# # # Output data.
# # tmp.file.name <- paste0(diagrams.path, '/DataForSharingByVennDiagrams_BasedOnCellsForAllSets.csv')
# # fwrite(file=tmp.file.name, x=diagrams.info, na='NA')
#
# # ---> Based on cells for 0 hrs. and 6 hrs. sets.
# # Get data.
# diagrams.info <- cells.clons.info[, .(time=orig.stim_time, cell.no=length(unique(barcode))), by=.(clone.id=raw_clonotype_id)]
# diagrams.info <- unique(diagrams.info)
# diagrams.info <- diagrams.info[time!=24]
# diagrams.info <- spread(data=diagrams.info, key=time, value=cell.no, fill=NA, sep='.')
# diagrams.info <- as.data.frame(diagrams.info, stringsAsFactors=FALSE)
# alt.set.names <- set.names[set.names!='time.24']
# diagrams.info$clone.size <- rowSums(diagrams.info[, alt.set.names], na.rm=TRUE)
# diagrams.info <- mclapply(X=1:nrow(diagrams.info), FUN=function(clone.row.id){
#   clones.clusts.info <- diagrams.info[clone.row.id, alt.set.names]
#   times.to.rep <- diagrams.info[clone.row.id, 'clone.size']
#   clones.clusts.info <- bind_rows(replicate(times.to.rep, clones.clusts.info, simplify = FALSE))
#   clones.clusts.info$clone.id <- paste(diagrams.info[clone.row.id, 'clone.id'], 1:times.to.rep, sep='.')
#   return(clones.clusts.info)
# })
# diagrams.info <- rbindlist(diagrams.info)
# venn.input <- lapply(X=alt.set.names, FUN=function(tmp.set){
#   to.output <- diagrams.info[!is.na(get(tmp.set)), clone.id]
#   return(to.output)
# })
# names(venn.input) <- paste0(time.points[time.points!=24], ' hrs.')
# # Output Venn Diagram.
# tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnCells_0and6hrs.png')
# venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2])
# tmp.file.name <- paste0(diagrams.path, '/SharingByVennDiagrams_BasedOnCells_0and6hrs_Blank.png')
# venn.diagram(x=venn.input, filename=tmp.file.name, fill=set.cols[1:2], col=set.cols[1:2], cex=0, cat.cex=0)
# # Output data.
# tmp.file.name <- paste0(diagrams.path, '/DataForSharingByVennDiagrams_BasedOnCellsForSets0and6hrs.csv')
# fwrite(file=tmp.file.name, x=diagrams.info, na='NA')

# To get rid of all logs produced by Venn.
tmp.command <- paste0('rm ', diagrams.path, '/*.log')
system(command=tmp.command)

# ---> Venn diagrams with jvenn.
# jvenn tool. Saved to the same folder.

############    -----------------------------------------    ############
############    ----------   UpSetR Plots    ------------    ############
############    -----------------------------------------    ############
upset.path <- paste0(reports.path, '/upset_plots')
if(!dir.exists(upset.path)) dir.create(upset.path)

############    -----------   UpSetR Plot    ------------    ############
############    ---------   Based on Clones    ----------    ############

upset.clones.path <- paste0(upset.path, '/based_on_clones')
if(!dir.exists(upset.clones.path)) dir.create(upset.clones.path)

for(freq.thold in freq.tholds){
  # ---> Format.
  time.sets.info <- cells.clons.info[, .SD[, .(cell.no=length(unique(barcode))), by=.(time=orig.stim_time)], by=.(clone.id=raw_clonotype_id)]
  time.sets.info <- spread(data=time.sets.info, key=time, value=cell.no, fill=0, sep='.')

  for(tmp.set in set.names){
    time.sets.info[, eval(tmp.set):=ifelse(test=get(tmp.set)>=freq.thold, yes=1, no=0)]
  }

  # ---> UpSetR Plot.
  tmp.file.name <- paste0(upset.clones.path, '/SharingAmongDataSetsWithThold_', freq.thold, '_UpSetR.pdf')
  pdf(file=tmp.file.name)
  print(upset(data=time.sets.info,
    sets=set.names,
    nintersects=NA,
    query.legend="top",
    order.by=c("freq", "degree"),
    # keep.order=TRUE,
    mainbar.y.label='Number of intersected clones',
    sets.x.label=paste0("Total clones (freq. > ", freq.thold, ")"),
    main.bar.color='#0080ff',
    sets.bar.color='#0d0d0d',
    # text.scale=c(0, 0, 0, 0, 0, 0), # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
    # scale.intersections=tmp.scale
  ))
  dev.off()
}

############    -----------   UpSetR Plot    ------------    ############
############    ----------   Based on Cells    ----------    ############

upset.cells.path <- paste0(upset.path, '/based_on_cells')
if(!dir.exists(upset.cells.path)) dir.create(upset.cells.path)

for(freq.thold in freq.tholds){
  # ---> Format.
  time.sets.info <- cells.clons.info[, .SD[, .(cell.no=length(unique(barcode))), by=.(time=orig.stim_time)], by=.(clone.id=raw_clonotype_id)]
  time.sets.info <- spread(data=time.sets.info, key=time, value=cell.no, fill=0, sep='.')

  # ---> Create clone-basis structure.
  # Remove clones not seen in the 6H data set with a clone size larger than 1 for at least one of the clusters.
  time.sets.info <- as.data.frame(time.sets.info, stringsAsFactors=FALSE)
  clons.to.keep <- sapply(X=1:nrow(time.sets.info), FUN=function(idx){
    clone.info <- time.sets.info[idx, set.names]
    to.eval <- sum(clone.info>=freq.thold)
    return(to.eval>0)
  })
  time.sets.info <- time.sets.info[clons.to.keep, ]
  # Then, get cell-basis structure of the data. In brief, a clone row is repeated as many times as its clone size.
  time.sets.info$clone.size <- rowSums(time.sets.info[, set.names])
  clones.info.cb <- mclapply(X=1:nrow(time.sets.info), FUN=function(clone.row.id){
    clones.clusts.info <- time.sets.info[clone.row.id, set.names]
    times.to.rep <- time.sets.info[clone.row.id, 'clone.size']
    clones.clusts.info <- bind_rows(replicate(times.to.rep, clones.clusts.info, simplify = FALSE))
    return(clones.clusts.info)
  })
  clones.info.cb <- rbindlist(clones.info.cb)
  # Proof checking we got it right.
  # Total clone size should equal the total amount of rows we got in the processed data.
  sum(rowSums(time.sets.info[, set.names]))
  nrow(clones.info.cb)

  for(tmp.set in set.names){
    clones.info.cb[, eval(tmp.set):=ifelse(test=get(tmp.set)>=freq.thold, yes=1, no=0)]
  }

  # ---> UpSetR Plot.
  tmp.file.name <- paste0(upset.cells.path, '/SharingAmongDataSetsWithThold_', freq.thold, '_UpSetR.pdf')
  pdf(file=tmp.file.name)
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
    # text.scale=c(0, 0, 0, 0, 0, 0), # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
    # scale.intersections=tmp.scale
  ))
  dev.off()
}


############    -----------------------------------------    ############
############    -----   Dimensionality reduction    -----    ############
############    -----------------------------------------    ############

dim.red.path <- paste0(reports.path, '/dim_reduction')
if(!dir.exists(dim.red.path)) dir.create(dim.red.path)

# ---> Main program as a function.
# Given a seurat object, this function will add the annotations regarding clone sharing among data sets per cell.
# Arguments:
#   reference.time: '0', '6' or '24'
#   barcodes.tag: Specific tag that indicates the barcodes from the set of interest.
#   time.sets.info as calculated below.
# Value: Annotated seurat object.

annotate.seurat.obj.ss <- function(seurat.obj, reference.time, barcodes.tag, time.sets.info){
  # ---> First part: Get clones sharing-related info.
  # To data frame.
  time.sets.info <- as.data.frame(time.sets.info, stringsAsFactors=FALSE)
  # Data sets different than the reference. Define accordingly.
  other.sets <- setdiff(x=time.points, y=reference.time)
  # Add the clone sharing status.
  time.sets.info$status <- ifelse(test=time.sets.info[, other.sets[1]]>0 & time.sets.info[, other.sets[2]]>0, yes='Both', no=
    ifelse(test=time.sets.info[, other.sets[1]]>0, yes=other.sets[1], no=
      ifelse(test=time.sets.info[, other.sets[2]]>0, yes=other.sets[2], no='NS')))
  time.sets.info[time.sets.info[, reference.time]==0, 'status'] <- NA
  # Subset to keep the info of interest.
  sharing.info <- time.sets.info[, c('clone.id', 'status')]
  sharing.info <- sharing.info[!is.na(sharing.info[, 'status']), ]
  # ---> Second part: Get clones-barcodes relationships.
  barcodes.info <- cells.clons.info[!is.na(get(barcodes.tag)), .(barcode=get(barcodes.tag), clone.id=raw_clonotype_id)]
  barcodes.info <- unique(barcodes.info)
  # ---> Third part: Merge all clones info according to cells order in seurat object.
  to.annotate <- merge(x=barcodes.info, y=sharing.info, by='clone.id')
  seurat.obj.info <- data.table(barcode=Cells(seurat.obj))
  to.annotate <- merge(x=seurat.obj.info, y=to.annotate, by='barcode', all.x=TRUE, all.y=FALSE, sort=FALSE)
  # ---> Four part: Annotate the seurat object accordingly.
  seurat.obj@meta.data[, 'ams.clonotype.id.tag'] <- to.annotate[, clone.id]
  seurat.obj@meta.data[, 'clone.status.tag'] <- to.annotate[, status]
  # Return.
  return(seurat.obj)
}


# ---> Annotations in seurat objects: Cells shared with a different data set.
time.sets.info <- cells.clons.info[, .SD[, .(cell.no=length(unique(barcode))), by=.(time=orig.stim_time)], by=.(clone.id=raw_clonotype_id)]
time.sets.info <- spread(data=time.sets.info, key=time, value=cell.no, fill=0)

# Combination aggr.
seurat.obj.comb <- annotate.seurat.obj.ss(seurat.obj=seurat.obj.comb, reference.time=0, barcodes.tag='comb.hrs.barcode', time.sets.info=time.sets.info)
# 6 hrs.
seurat.obj.six <- annotate.seurat.obj.ss(seurat.obj=seurat.obj.six, reference.time='6', barcodes.tag="six.hrs.barcode", time.sets.info=time.sets.info)
# 24 hrs.
seurat.obj.tf <- annotate.seurat.obj.ss(seurat.obj=seurat.obj.tf, reference.time='24', barcodes.tag="tf.hrs.barcode", time.sets.info=time.sets.info)

# ---> Dimensionality reduction.
# Function.
depict.sharing <- function(seurat.obj, description, reports.path){
  # All cells.
  tmp.file.name <- paste0(reports.path, '/SharingFor', description, 'AggregationOnUMAP_AllCells.pdf')
  pdf(file=tmp.file.name)
  print(DimPlot(object=seurat.obj, reduction='umap', group.by='clone.status.tag'))
  dev.off()
  # Only shared cells.
  to.depict <- Cells(seurat.obj)[seurat.obj@meta.data[, 'clone.status.tag']!='NS' & !is.na(seurat.obj@meta.data[, 'clone.status.tag'])]
  tmp.file.name <- paste0(reports.path, '/SharingFor', description, 'AggregationOnUMAP_SharedOnly.pdf')
  pdf(file=tmp.file.name)
  print(DimPlot(object=seurat.obj, reduction='umap', group.by='clone.status.tag', cells=to.depict))
  dev.off()
}

# Combination aggr.
depict.sharing(seurat.obj=seurat.obj.comb, description='Combination', reports.path=dim.red.path)
# 6 hrs.
depict.sharing(seurat.obj=seurat.obj.six, description='Six', reports.path=dim.red.path)
# 24 hrs.
depict.sharing(seurat.obj=seurat.obj.tf, description='Twentyfour', reports.path=dim.red.path)

############    -----------------------------------------    ############
############    -   Clones shared between 0 and 6hrs    -    ############
############    -----------------------------------------    ############

saz.sharing.path <- paste0(reports.path, '/six_and_zero_sharing')
if(!dir.exists(saz.sharing.path)) dir.create(saz.sharing.path)

tmp.tag <- 'clone.status.tag'

# ---> Function to combine preprocessed data and provide outputs.
sets.cell.sharing <- function(six.hrs.clones, zero.hrs.clones, file.desc){
  # Update the tag.
  six.hrs.clones[, status:=ifelse(test=status=='NS', yes='NS', no=ifelse(status=='0', yes='S', no=ifelse(test=status=='Both', yes='S', no='NS')))]
  zero.hrs.clones[, status:=ifelse(test=status=='NS', yes='NS', no=ifelse(status=='6', yes='S', no=ifelse(test=status=='Both', yes='S', no='NS')))]
  # Combine data from both sets.
  tmp.data <- rbindlist(l=list(zero.hrs.clones, six.hrs.clones), idcol='Set')
  tmp.data[, Set:=ifelse(Set==1, yes='Zero', no='Six')]

  # ---> Shared vs non-shared CELLS per data set.

  # @ Proportion of Cells.
  # Calculate proportion of cells to use it as a label.
  tmp.text <- data.frame(Set=c('Zero', 'Six'), text=c(paste0(tmp.data[status=='S' & Set=='Zero', .N], '/', tmp.data[Set=='Zero', .N]), paste0(tmp.data[status=='S' & Set=='Six', .N], '/', tmp.data[Set=='Six', .N])), stringsAsFactors=FALSE)
  # Plot and output.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=Set, fill=status)) + geom_bar(position='fill', width=0.6) + geom_text(data=tmp.text, aes(x=Set, y=0.9, label=text), inherit.aes=FALSE) + scale_fill_manual(values=tmp.cols) + ylab('Proportion') + theme_minimal()
  tmp.file.name <- paste0(saz.sharing.path, '/SharingPerSetPerCellBarPlot_', file.desc, '.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
  # Blank version.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=Set, fill=status)) + geom_bar(position='fill', width=0.6) + scale_fill_manual(values=tmp.cols) + ylab('Proportion') + theme_minimal() + blank.complement
  tmp.file.name <- paste0(saz.sharing.path, '/SharingPerSetPerCellBarPlot_', file.desc, '_Blank.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()

  # @ Clonal expansion of Cells.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=Set, fill=status, y=size)) + geom_violin() + scale_y_log10() + theme_minimal() + ylab('Clone size') + scale_fill_manual(values=tmp.cols)
  tmp.file.name <- paste0(saz.sharing.path, '/SizeBySharingAndSetBarPlot_', file.desc, '.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()

  # Mock return.
  return(NA)
}

# ---> All cells from all donors.
# Pick important data for both.
# Consider only the cells coming from donors shared by both datasets.
# @ Six hrs. data set.
six.hrs.clones <- as.data.table(seurat.obj.six@meta.data)
six.hrs.clones[, cell:=Cells(seurat.obj.six)]
six.hrs.clones <- six.hrs.clones[!is.na(clone.status.tag) & orig.virus2=='CV', .(cell=cell, status=clone.status.tag, cluster=clusters.tag, size=clon.size.tag)]
# @ Combination aggr. data set.
zero.hrs.clones <- as.data.table(seurat.obj.comb@meta.data)
zero.hrs.clones[, cell:=Cells(seurat.obj.comb)]
zero.hrs.clones <- zero.hrs.clones[!is.na(clone.status.tag) & orig.stim_time==0, .(cell=cell, status=clone.status.tag, cluster=clusters.tag, size=clon.size.tag)]
# Get output.
sets.cell.sharing(six.hrs.clones=six.hrs.clones, zero.hrs.clones=zero.hrs.clones, file.desc='AllDonors_AllCells')

# ---> All cells from matched donors.
# Pick important data for both.
# Consider only the cells coming from donors shared by both datasets.
# @ Six hrs. data set.
six.hrs.clones <- as.data.table(seurat.obj.six@meta.data)
six.hrs.clones[, cell:=Cells(seurat.obj.six)]
six.hrs.clones <- six.hrs.clones[!is.na(clone.status.tag) & orig.virus2=='CV' & orig.donor %chin% shared.donors$donor, .(cell=cell, status=clone.status.tag, cluster=clusters.tag, size=clon.size.tag)]
# @ Combination aggr. data set.
zero.hrs.clones <- as.data.table(seurat.obj.comb@meta.data)
zero.hrs.clones[, cell:=Cells(seurat.obj.comb)]
zero.hrs.clones <- zero.hrs.clones[!is.na(clone.status.tag) & orig.stim_time==0 & orig.donor %chin% shared.donors$donor, .(cell=cell, status=clone.status.tag, cluster=clusters.tag, size=clon.size.tag)]
# Get output.
sets.cell.sharing(six.hrs.clones=six.hrs.clones, zero.hrs.clones=zero.hrs.clones, file.desc='MatchedDonors_AllCells')

# ---> All cells from matched donors.
# Pick important data for both.
# Consider only the cells coming from donors shared by both datasets.
# @ Six hrs. data set.
six.hrs.clones <- as.data.table(seurat.obj.six@meta.data)
six.hrs.clones[, cell:=Cells(seurat.obj.six)]
six.hrs.clones <- six.hrs.clones[cell %chin% six.sample, .(cell=cell, status=clone.status.tag, cluster=clusters.tag, size=clon.size.tag)]
# @ Combination aggr. data set.
zero.hrs.clones <- as.data.table(seurat.obj.comb@meta.data)
zero.hrs.clones[, cell:=Cells(seurat.obj.comb)]
zero.hrs.clones <- zero.hrs.clones[cell %chin% zero.sample, .(cell=cell, status=clone.status.tag, cluster=clusters.tag, size=clon.size.tag)]
# Get output.
sets.cell.sharing(six.hrs.clones=six.hrs.clones, zero.hrs.clones=zero.hrs.clones, file.desc='MatchedDonors_Subsample')


# ---> Cluster info for shared cells.

# Pick data (only shared cells) and give appropriate format.
# @ Six hrs.
tmp.six <- six.hrs.clones[status=='S']
tmp.six <- as.data.frame(tmp.six[, .(cell.no=.N), by=.(cluster, status)], stringsAsFactors=FALSE)
tmp.six$cluster <- factor(tmp.six$cluster, levels=mixedsort(unique(tmp.six$cluster)))
# @ Six hrs.
tmp.zero <- zero.hrs.clones[status=='S']
tmp.zero <- as.data.frame(tmp.zero[, .(cell.no=.N), by=.(cluster, status)], stringsAsFactors=FALSE)
tmp.zero$cluster <- factor(tmp.zero$cluster, levels=mixedsort(unique(tmp.zero$cluster)))
# @ Preprocessing.
# Compute percentages
tmp.six$fraction <- tmp.six$cell.no/sum(tmp.six$cell.no)
tmp.zero$fraction <- tmp.zero$cell.no/sum(tmp.zero$cell.no)
# Compute the cumulative percentages (top of each rectangle)
tmp.six$ymax <- cumsum(tmp.six$fraction)
tmp.zero$ymax <- cumsum(tmp.zero$fraction)
# Compute the bottom of each rectangle
tmp.six$ymin <- c(0, head(tmp.six$ymax, n=-1))
tmp.zero$ymin <- c(0, head(tmp.zero$ymax, n=-1))
# Plot and output, 6 hrs.
tmp.ggplot <- ggplot(data=tmp.six, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=cluster)) + geom_rect() + coord_polar(theta="y") + xlim(c(2, 4)) + theme_minimal() + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line=element_blank())
# tmp.ggplot <- ggplot(data=tmp.six, aes(x=status, fill=cluster)) + geom_bar(position='fill') + coord_polar(theta="y") + xlim(c(2, 4))
tmp.file.name <- paste0(saz.sharing.path, '/6hrsSharedCellsClusterInfo.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()
# Plot and output, 0 hrs.
tmp.ggplot <- ggplot(data=tmp.zero, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=cluster)) + geom_rect() + coord_polar(theta="y") + xlim(c(2, 4)) + theme_minimal() + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.line=element_blank())
# tmp.ggplot <- ggplot(data=tmp.six, aes(x=status, fill=cluster)) + geom_bar(position='fill') + coord_polar(theta="y") + xlim(c(2, 4))
tmp.file.name <- paste0(saz.sharing.path, '/0hrsSharedCellsClusterInfo.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()


############    -----------------------------------------    ############
### ----------------------------- Cites ----------------------------- ###
############    -----------------------------------------    ############
# ---> jvenn
# Cite: Philippe Bardou, Jérôme Mariette, Frédéric Escudié, Christophe Djemiel and Christophe Klopp. jvenn: an interactive Venn diagram viewer. BMC Bioinformatics 2014, 15:293 doi:10.1186/1471-2105-15-293
# Web page: http://jvenn.toulouse.inra.fr/app/example.html
# ---> UpSetR
# Jake R Conway, Alexander Lex, Nils Gehlenborg UpSetR: An R Package for the Visualization of Intersecting Sets and their Properties doi: https://doi.org/10.1093/bioinformatics/btx364
