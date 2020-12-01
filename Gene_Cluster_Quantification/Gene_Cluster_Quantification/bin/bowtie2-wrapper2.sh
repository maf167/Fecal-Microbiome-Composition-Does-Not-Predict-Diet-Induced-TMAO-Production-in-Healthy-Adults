#!/bin/sh

SLDIR=$1
PDIR=$2
SCR_HOME=$3
FASTA=$4
CTS=$5
TCTS=$6

FBASE=$(echo $FASTA | awk -F/ '{print $NF}' | cut -d "." -f 1)

for i in `seq 1  4`; do

  S="""
#!/bin/bash\n
#SBATCH --mail-type=ALL\n
#SBATCH --mail-user=ferrelm@ccf.org\n
#SBATCH --time=600:00:00\n
#SBATCH --job-name=Job${i}\n
\n
echo "BEGIN"\n
date\n
\n
module load samtools \n
module load bowtie2/2.3.1 \n
\n
FASTA=${FASTA}\n
PDIR=${PDIR}\n
\n
FOR_READS1=${SCR_HOME}/data/fqlist/approach_1_1_fq.gz.names\n
FOR_READS2=${SCR_HOME}/data/fqlist/approach_1_2_fq.gz.names\n
FOR_READS3=${SCR_HOME}/data/fqlist/approach_1_3_fq.gz.names\n
FOR_READS4=${SCR_HOME}/data/fqlist/approach_1_4_fq.gz.names\n
REV_READS1=${SCR_HOME}/data/fqlist/approach_2_1_fq.gz.names\n
REV_READS2=${SCR_HOME}/data/fqlist/approach_2_2_fq.gz.names\n
REV_READS3=${SCR_HOME}/data/fqlist/approach_2_3_fq.gz.names\n
REV_READS4=${SCR_HOME}/data/fqlist/approach_2_4_fq.gz.names\n
\n
\n
l=\$(awk -F, '{print NF}' \${REV_READS${i}})\n
\n
for i in \`seq 1 \$l\`; do \n
  REV=\$(awk -v i=\${i} -F, '{print \$i}' \$REV_READS${i})\n
  FOR=\$(awk -v i=\${i} -F, '{print \$i}' \$FOR_READS${i})\n
  SAMP=\$(echo \$FOR | awk -F/ '{print \$9}')\n
\n
  basename=\$(echo \$FASTA | awk -F/ '{print \$NF}' | cut -d "." -f 1)
  bowtie2-build \$FASTA \${PDIR}/int/refdata/btind/\${basename}
  bowtie2 -x \${PDIR}/int/refdata/btind/\${basename} -1 \$FOR -2 \$REV --phred64 --very-sensitive-local --no-unal | samtools view -bS | samtools sort - | tee >(samtools idxstats - > \${CTS}/\${basename}.cts) | samtools view -cF 0x100 > ${TCTS}/\${basename}.totcts
done\n

echo "END"\n
date\n
  """
  echo -e $S > ${SLDIR}/${FBASE}_${i}_job.sh
  sbatch ${SLDIR}/${FBASE}_${i}_job.sh
done