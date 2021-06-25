#!/bin/bash
######################
# MONOCLE 3 WRAPPER ##
######################

## DESCRIPTION ##
# Wrapper for Monocle 3 analysis


#### Default parameters #### ---------------------------------------------------
MYCELLSF=${MYCELLSF}
RAWDATA=${RAWDATA:-${MYCELLSF}}
OUTDIR=${OUTDIR:-/mnt/BioScratch/cramirez}
CATG=${CATG:-res.0.6}
CATGFILT=${CATGFILT:-mycolumn~myclass}
SMPCELLS=${SMPCELLS:-FALSE}
CPMCONV=${CPMCONV:-FALSE}
GCOLOUR=${GCOLOUR:-ningun_file}
TOPGENES=${TOPGENES:-16}
MYSEED=${MYSEED:-27}
# specific paramaters #
REDTYPE=${REDTYPE:-tSNE} # #c("DDRTree", "ICA", "tSNE", "UMAP")
SEUDISPS=${SEUDISPS:-TRUE}
PCS_COMP=${PCS_COMP:-50}
RESMF=${RESMF:-NULL}
GENESETF=${GENESETF:-no_file}
CTSTYPE=${CTSTYPE:-raw}
FILTCELLS=${FILTCELLS:-TRUE}
# RESMF=${RESMF:-nUMI~percent.mito}
MEM=${MEM:-50gb}

echo; echo "**** Ciro R-Suastegui - LJI 2018 (v1)"
echo "------------------------------- PRESENTING PARAMETERS -------------------------------"
echo "Annotation file: ${MYCELLSF}"
echo "Expression matrix: ${RAWDATA}"
echo "Output directory: ${OUTDIR}"
echo "Columns for groups: ${CATG}"
echo "Filter: ${CATGFILT}"
echo "Sample cell for test: ${SMPCELLS}"
echo "CPM convert: ${CPMCONV}"
echo "No. of markers: ${TOPGENES}"
echo "Colours file: ${GCOLOUR}"
echo "Seed: ${MYSEED}"
# specific paramaters #
echo "Reduction types: ${REDTYPE}"
echo "Use Seurat's HVGs: ${SEUDISPS}"
echo "Number of PCs: ${PCS_COMP}"
echo "Regression Fromula: ${RESMF}"
echo "Selected genes: ${GENESETF}"
echo "Counts type: ${CTSTYPE}"
echo "Filter cells: ${FILTCELLS}"
echo "------------------------------- --------------------- -------------------------------"
read -n 1 -s -r -p "Press any key to continue"
echo

# # Clean and create directory
if [ ! -d ${OUTDIR} ]; then mkdir -p ${OUTDIR}; fi
clea=`ls ${OUTDIR} | wc -l`
if [ ${clea} -gt 0 ]; then
  delim="=========== ${OUTDIR}: ==========="
  echo $delim
  echo -e "Contains:\n`ls ${OUTDIR}`"
  printf "%0.s=" $(seq 1 ${#delim}); echo
  while true; do
    read -p "Do you wish to erase previous information and continue with the anaylisis? (y/n/s): " yn
    case $yn in
      [Yy]* ) read -p "Are you completely sure? " yn
        case $yn in
          [Yy]* ) echo "Cleaning..."; rm -r ${OUTDIR}/*; break;;
          [Nn]* ) exit;;
          [Ss]* ) break;;
          * ) echo "Please answer yes (Y/y) or no (N/n) or skip (S/s).";;
        esac
        ;;
      [Nn]* ) echo "Aborting..."; exit;;
      [Ss]* ) break;;
      * ) echo "Please answer yes (Y/y) or no (N/n) or skip (S/s).";;
    esac
  done
fi
if [ ! -d ${OUTDIR}/jobs ]; then mkdir ${OUTDIR}/jobs; fi
cd ${OUTDIR}
echo "Working at: `pwd`";echo
echo "Running analysis:"
SUBNAME=${REDTYPE}_${CATG:0:7}_${PCS_COMP}PC
if [ "${SEUDISPS}" == "TRUE" ]; then SUBNAME=${SUBNAME}_SEU; else SUBNAME=${SUBNAME}_MONO; fi
FNAME=${OUTDIR}/jobs/${SUBNAME}
if [ -s ${FNAME}.out.txt ]; then rm ${FNAME}.out.txt; fi
cat <<EOT > ${FNAME}.sh
#PBS -N ${SUBNAME}
#PBS -o ${FNAME}.out.txt
#PBS -e ${OUTDIR}/jobs/a1.err.txt
#PBS -m a
#PBS -M ciro@lji.org
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -l mem=${MEM}
#PBS -l walltime=12:00:00

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat ${PBS_NODEFILE}
echo; echo ------------------------------------------------------
echo "PBS: qsub is running on ${PBS_O_HOST}"
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = ${PBS_O_PATH}
echo ------------------------------------------------------
echo; echo

cd ${OUTDIR}
/share/apps/R/3.5/bin/Rscript /mnt/BioHome/ciro/scripts/monocle/monocle3_v1.R \\
          --mycellsf=${MYCELLSF} \\
          --rawdata=${RAWDATA} \\
          --outdir=${OUTDIR} \\
          --catg=${CATG} \\
          --catgfilt="${CATGFILT}" \\
          --smpcells=${SMPCELLS} \\
          --cpmconv=${CPMCONV} \\
          --gcolour=${GCOLOUR} \\
          --topgenes=${TOPGENES} \\
          --myseed=${MYSEED} \\
          --redtype=${REDTYPE} \\
          --seudisps=${SEUDISPS} \\
          --pcs_comp=${PCS_COMP} \\
          --resMF=${RESMF} \\
          --genesetf="${GENESETF}" \\
          --ctstype=${CTSTYPE} \\
          --filtcells=${FILTCELLS} \\
          --myseed=${MYSEED}
echo
EOT
qsub ${FNAME}.sh
echo
echo "Job file: ${FNAME}.sh"
echo "Good luck!"
echo
