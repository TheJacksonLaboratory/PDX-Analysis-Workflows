#!/bin/bash

# A script to Add Caller column

export file=$1
export output1=$2




awk -F$'\t'  '/^[^#]/ { print $1,$2,$3,$4,$5,$6,$7,$8";CALLER=Pindel",$9,$10;next } {print $0}' $file|tr ' ' '\t' > $output1

  sed -i '1 a\##INFO=<ID=CALLER,Number=1,Type=String,Description="The variant caller Used to call the variant GATK UG or Pindel">\ ' $output1