#!/bin/bash
######################
# MONOCLE 3 WRAPPER ##
######################

## DESCRIPTION ##
# Wrapper for Monocle 3 analysis


#### parameters #### -----------------------------------------------------------
MYCELLSF=/home/ciro/large/covid19/results/clustering/CD4T6_R1n2n3_sng_25p/clustering/zetInfo/clustCells38PCs_30Ks_0.06667JD.RData # metadata, CSV or Seurat object
# RAWDATA # 10X folder PATH (before outs), CSV or Seurat object (takes MYCELLSF if not set)
OUTDIR=/home/ciro/large/covid19/results/trajectory/m3alpha_`basename ${MYCELLSF/clustering\/zet*/}` # =${OUTDIR:-/mnt/BioScratch/cramirez/} # dont put "/" at the end
CATG=RNA_snn_res.0.6~orig.virus2~orig.hospital~orig.severity # =${CATG:-res.06} # colunm name in meta data or annotation file
CATGFILT=orig.virus2~CV # =${CATGFILT:-mycolumn~myclass} # a single char of mycolumn~myclass or a "list(c('column1','class2'), c('column14','-class1'))"
# SMPCELLS # =${SMPCELLS:-FALSE} # if cells should be sampled to the smaller group
# CPMCONV # =${CPMCONV:-TRUE} # CPM conversion for plots
GCOLOUR=/home/ciro/covid19/info/global_colours.csv # =${GCOLOUR:-ningun_file} # header: group,colour
# TOPGENES # =${TOPGENES:-16} # top genes to show in plots
# MYSEED # =${MYSEED:-27} # seed for determinism
# specific paramaters #
REDTYPE=DDRTree # =${REDTYPE:-tSNE} # type of dimReduction, options: "DDRTree", "ICA", "tSNE", "UMAP"
SEUDISPS=TRUE # =${SEUDISPS:-TRUE} # use seurat HVGs
PCS_COMP=38 # =${PCS_COMP:-2} # number of PCs to use
# RESMF # =${RESMF:-nUMI~percent.mito} # regressing variables
# GENESETF # =${GENESETF:-no_file} # file with genes or a cutoff for monocle mean_expression
# CTSTYPE #=${CTSTYPE:-raw} # type of data
FILTCELLS=FALSE #=${FILTCELLS:-TRUE} filter cells according to distribution
#### parameters #### -----------------------------------------------------------

MEM=80gb
OUTDIR=/home/ciro/large/covid19/results/trajectory/m3alpha_`basename ${MYCELLSF/clustering\/zet*/}`_tfh0n5_200
CATGFILT="list(c('orig.virus2', 'CV'), c('RNA_snn_res.0.6', '0', '5'))"
GENESETF="c('1', '200')"

source monocle_pipeline.sh
