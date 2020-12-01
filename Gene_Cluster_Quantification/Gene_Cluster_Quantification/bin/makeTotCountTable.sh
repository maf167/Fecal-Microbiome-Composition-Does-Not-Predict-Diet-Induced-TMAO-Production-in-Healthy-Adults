#!/bin/sh

FOR_FILES=$1
l=$(awk -F, '{print NF}' $FOR_FILES)

CTDIR=$2
HDIR=$PWD

cd $CTDIR

>tmpSamp
for i in `seq 1 $l`; do

  FOR=$(awk -v i=${i} -F, '{print $i}' $FOR_FILES)
  SAMP=$(echo $FOR | awk -F/ '{print $9}')
  echo $SAMP >> tmpSamp
done

for i in $(ls | grep "tot.counts"); do

  SAMP=$(echo ${i} | awk -F_ '{print $1}')
  CT=$(cat $i)

  echo -e "${SAMP}\t${CT}"

done

while read line; do
  f=$(ls | grep "${line}")
  nl=$(wc -l $f | cut -f 1 -d " ")
  # Skip counts files with no data
  if [[ $nl -gt 1 ]]; then
    SAMP=$(echo ${f} | awk -F_ '{print $1}')
    CT=$(cat $f)
    echo -e "${SAMP}\t${CT}"
  fi
done<tmpSamp

rm tmpSamp

cd$HDIR