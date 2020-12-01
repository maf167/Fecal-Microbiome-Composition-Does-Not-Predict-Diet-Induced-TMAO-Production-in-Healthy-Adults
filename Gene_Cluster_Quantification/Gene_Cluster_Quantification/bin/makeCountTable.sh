#!/bin/sh

## Get sample names from list of fastq files

FOR_FILES=$1
l=$(awk -F, '{print NF}' $FOR_FILES)

CTDIR=$2
HDIR=$PWD

cd CTDIR

>tmpSamp
for i in `seq 1 $l`; do

  FOR=$(awk -v i=${i} -F, '{print $i}' $FOR_FILES)
  SAMP=$(echo $FOR | awk -F/ '{print $9}')
  echo $SAMP >> tmpSamp
done

## Take header from an arbitrary counts file
HEAD=$(awk 'BEGIN{ ORS = "\t"}{print $1}' $(ls | grep "1000A"))
echo -e "\t${HEAD}"

## Make the second column of each counts file a row in the table

while read line; do
  f=$(ls | grep "${line}")
  nl=$(wc -l $f | cut -f 1 -d " ")
  # Skip counts files with no data
  if [[ $nl -gt 1 ]]; then
    l=$(awk 'BEGIN{ ORS = "\t"}{print $3}' $f)
    echo -e "${line}\t${l}"
  fi
done<tmpSamp

rm tmpSamp

cd $HDIR