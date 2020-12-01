#!/bin/sh
# FASTA MUST BE FORMATTED WITH SINGLE LINE SEQS
# BEDTOOLS-created

# PRINTS TO stdout

FASTA=$1

# Chomp '>' from Sequence IDs 

seqs=($(awk -F% '$1 ~ /^>/ {print substr($1,2)}' $FASTA))

for i in "${seqs[@]}"
do
  # Sequence is the line after the ID Line, all whitespace removed 
  ISEQ=$(grep -A1 "^>${i}$" $FASTA | tail -n1 | tr -d '[:space:]')
  echo -e "${i}\t${#ISEQ}"
done