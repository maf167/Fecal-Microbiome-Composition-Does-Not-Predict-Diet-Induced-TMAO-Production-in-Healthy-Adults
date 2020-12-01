#!/bin/sh
# CREATE GENOMEDIR WITH GENOME FASTA FILES AND EXTRACT WITH GIVEN BED FILES


STRAINLIST=$1
GENOMEDIR=$2
IFS=',' read -r -a BEDS <<< "$3"
FNADIR=$4
TMPDIR=$5
HDIR=$PWD

# Download KEGG GENOME Catalog webpage that links to GenBank

wget -O ${TMPDIR}/tmpFTPs genome.jp/kegg/catalog/org_list.html

>FTPtab
while read line;do

  # Extract GenBank FTP directory from webpage
  SITE=$(grep -A5 ">${line}</a>" ${TMPDIR}/tmpFTPs | grep "ftp" | sed -e 's/href=\x27/\n/' | tail -n1 |  sed -e 's/\x27>/\n/' | head -n1)
  
  # Specify the genomic FASTA file within GenBank directory
  EXT=$(echo $SITE | awk -F/ '{print $NF}')
  LINK="${SITE}/${EXT}_genomic.fna.gz"
  
  # Write to FTPtab
  echo -e "${line}\t${LINK}" >> ${TMPDIR}/FTPtab

done<$STRAINLIST

# Download FASTA Genomes
while read line; do

  link=$(grep "^${line}" ${TMPDIR}/FTPtab | cut -f 2)                      # Lookup download link
  fname=$(echo $link | awk -F/ '{print $NF}')
  fname2=$(echo $fname | awk -F. 'OFS="." {$NF=""; print$0}' | sed 's/.$//')
  cd $GENOMEDIR
  wget  $link                                                              # Perform download
  gzip -d ${fname}                                                         # Decompress file
  nl=$(wc -l ${fname2} | awk '{print $1}')                                 # Get line count
  nll=$(echo "${nl} - 1" | bc)
  tail -n ${nll} ${fname2} > ${TMPDIR}/tmp1                                # tmp1 is missing first line
  echo ">${line}" | cat - ${TMPDIR}/tmp1 > ${TMPDIR}/tmp2
  mv ${TMPDIR}/tmp2 ${line}_wgs.fna                                        # new first line with label matching BED file
  rm ${fname2}

done<$STRAINLIST

rm GCA*
rm GCF*

cat ./* > all.fna

## Perform FASTA extraction with bedtools

module load bedtools/2.26.0

for bed in "${BEDS[@]}"; do
  basename=$(echo $bed | awk -F/ '{print $NF}' | cut -d "." -f 1)
  bedtools getfasta -fi all.fna -bed $bed -name -fo ${FNADIR}/${basename}.fna
done

cd $HDIR

## Clean Up Temporary Files

rm ${TMPDIR}/tmp1
rm ${TMPDIR}/tmp2
rm ${TMPDIR}/tmpFTPs
rm ${TMPDIR}/FTPtab