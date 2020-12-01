#!/bin/sh

if [[ $# -lt 5 ]]; then
  echo "Usage: ./extractOPfeatures.sh <ssdb.csv> <cluster name> <gap limit> <read list> <mate list> [KIDs]"
  echo "ssdb file         - from Kegg SSDB gene cluster tool in CSV"
  echo "gene cluster name - arbitrary label"
  echo "gap size          - limit on gaps between genes in junction profiling"
  echo "read list         - comma separated list of fastq files"
  echo "mate list         - same for reverse reads in pe sequencing"
  echo "KIDs              - comma separated list KEGG Genes"
  exit
fi

## SCRIPT ARGUMENTS

SSDB=$1
ONAME=$2
X=$3
FOR_READS=$4
REV_READS=$5
KIDS=$6

## PROJECT DIRECTORY SETUP

SCRIPT=$(readlink -f "$0")
SCR_HOME=$(dirname "$SCRIPT")

cd $SCR_HOME

PDIR=${SCR_HOME}/Projects/${ONAME}

if [ ! -d "${PDIR}" ]; then mkdir $PDIR; fi
mkdir ${PDIR}/inputs
mkdir ${PDIR}/outputs
mkdir ${PDIR}/int
mkdir ${PDIR}/int/bed
mkdir ${PDIR}/int/count
mkdir ${PDIR}/int/tmp
mkdir ${PDIR}/int/refdata
mkdir ${PDIR}/int/refdata/btind
mkdir ${PDIR}/int/refdata/genomes

GENOMEDIR=${PDIR}/int/refdata/genomes


### Run the makeBED.R script

 # Uses the REST API to search KEGG GENOME for genes from SSDB Results
 # Produces bed file and strain IDs to be used in extracting FASTAs of genes and junctions
module load R
Rscript ${SCR_HOME}/bin/makeBED.R $SSDB ${PDIR}/int/bed/${ONAME} $KIDS

### Exclude BGCs not physically close on genome

awk '$3 - $2 > 40000 || $2 == "NA" || $3 == "NA" {print $1}' ${PDIR}/int/bed/${ONAME}bgc.bed > ${PDIR}/int/tmp/exclude.txt

for b in ${PDIR}/int/bed/*.bed; do
  base=${b::-4}
  grep -vwFf ${PDIR}/int/tmp/exclude.txt $b > ${PDIR}/int/bed/${base}_inc.bed
done

grep -vwFf ${PDIR}/int/tmp/exclude.txt ${PDIR}/int/bed/${ONAME}Kstrains.txt > ${PDIR}/int/bed/${ONAME}Kstrain_inc.txt

### Convert BED files to FASTA

BEDSARR=()
for i in ${PDIR}/int/bed/*_inc.bed; do BEDSARR+=("$i"); done
delim=""
BEDS=""
for i in "${BEDSARR[@]}"; do
  BEDS="${BEDS}${delim}${i}"
  delim=","
done

${SCR_HOME}/bin/BED-to-FASTA.sh ${PDIR}/int/bed/${ONAME}Kstrain_inc.txt $GENOMEDIR $BEDS ${PDIR}/int/refdata ${PDIR}/int/tmp

### Align reads to each FASTA Reference

for f in ${PDIR}/int/refdata/*.fna; do
  
  ### Create subdirs for each FASTA reference
  base=${f::-4}
  mkdir ${PDIR}/int/count/${base}
  mkdir ${PDIR}/int/count/${base}/counts
  mkdir ${PDIR}/int/count/${base}/tot.counts
  
  ### Get Lengths
  
  ${SCR_HOME}/bin/getLens.sh $f > ${PDIR}/outputs/${base}_lens.txt
  
  ### Create SLURM Scripts Splitting Sample Reads Across Four Nodes
  # Count files per sample written to counts directory
  
  ${SCR_HOME}/bin/bowtie2-wrapper.sh ${PDIR}/int/count $PDIR $SCR_HOME $f ${PDIR}/int/count/${base}/counts ${PDIR}/int/count/${base}/tot.counts
  
  ### Merge cts into count table

  ${SCR_HOME}/bin/makeCountTable.sh    ${SCR_HOME}/data/fqlist/approach_1_1_fq.gz.names ${PDIR}/int/count/${base}/counts     > ${PDIR}/outputs/${base}_cts.txt
  ${SCR_HOME}/bin/makeTotCountTable.sh ${SCR_HOME}/data/fqlist/approach_1_1_fq.gz.names ${PDIR}/int/count/${base}/tot.counts > ${PDIR}/outputs/${base}_totcts.txt
  
done


