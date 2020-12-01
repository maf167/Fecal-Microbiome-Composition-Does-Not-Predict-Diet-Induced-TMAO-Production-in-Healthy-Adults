#!/bin/sh
# BOWTIE WRAPPER WITH OPTION TO SPLIT OVER 4 COMPUTE NODES

module load bowtie2/2.3.1

FASTA=$1
FOR=$2
REV=$3
INDDIR=$4
CTS=$5
TCTS=$6

basename=$(echo $bed | awk -F/ '{print $NF}' | cut -d "." -f 1)
bowtie2-build $FASTA ${INDDIR}/basename

bowtie2 -x ${INDDIR}/basename -1 $FOR -2 $REV --phred64 --very-sensitive-local --no-unal | samtools view -bS | samtools sort - | tee >(samtools idxstats - > ${CTS}/basename.cts) | samtools view -cF 0x100 > ${TCTS}/basename.totcts
