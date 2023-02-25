#!/bin/bash

############################################################
# Global settings                                          #
############################################################

# Abort the script at the first error, when a command exits with non-zero status
set -e

# Default parameters
adapterlength=44
pctmatch=97
threads=8
outdir=$(pwd)
mem=8G

Help()
{
   # Display Help
   echo "Filter adapters in HiFi reads."
   echo "Requires: 1) bam file to be filtered; 2) blast database of adapters"
   echo "Modified from HiFiAdapterFilt (https://github.com/sheinasim/HiFiAdapterFilt)."
   echo "Cite them: Sim, S.B., Corpuz, R.L., Simmonds, T.J. et al. \
    HiFiAdapterFilt, a memory efficient read processing pipeline, \
    prevents occurrence of adapter sequence in PacBio HiFi reads and \
    their negative impacts on genome assembly. \
    BMC Genomics 23, 157 (2022). \
    https://doi.org/10.1186/s12864-022-08375-1"
   echo "Usage:"
   echo $0
   echo "-b bam file"
   echo "-d database path"
   echo "[ -l minimum match Length to filter. Default=44 ]"
   echo "[ -m minimum Match percentage to filter. Default=97]  "
   echo "[ -t number of Threads for blastn. Default=8 ] "
   echo "[ -o output directory Default=./]"
   echo "[ -M memory Default=8G]"
   exit
}

while getopts ":b:d:l:m:t:o:M:h" option; do
   case $option in
      h) # display Help
         Help;;
      b)
        bam=${OPTARG};;
      d)
        db=${OPTARG};;
      l)
        adapterlength=${OPTARG};;
      m)
        pctmatch=${OPTARG};;
      t)
        threads=${OPTARG} ;;
      o)
        outdir=${OPTARG} ;;
      M)
        mem=${OPTARG} ;;
   esac
done

############################################################
# Checkers                                                 #
############################################################

# 1. if no option provided, print help
if [ $# -eq 0 ]; then
    Help
fi

# 2. if bam or db not found, print error message

if [ -z $bam ]; then
    echo "Input bam file is not set"
    exit
fi
if [ ! -f $bam ]; then
	echo $bam" is not found"
	exit
fi
if [ -z $db ]; then
    echo "Database directory is not set"
    exit
fi
if [ ! -d $db ]; then
	echo $db" is not found"
	exit
fi

## Create out directory if necessary
if [ ! -d ${outdir} ]
then 
	mkdir -p ${outdir}
fi

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

bam_prefix="$(basename $bam .bam)"
fastq=$outdir/$bam_prefix.fastq
fasta=$outdir/$bam_prefix.fasta
filt_fastq=$outdir/$bam_prefix.filt.fastq.gz

bamtools convert -format fastq -in $bam -out $fastq &
bamtools convert -format fasta -in $bam -out $fasta &
wait

blastn -db $db/pacbio_vectors_db \
  -query $fasta \
  -num_threads ${threads} \
  -task blastn \
  -reward 1 -penalty -5 -gapopen 3 -gapextend 3 \
  -dust no \
  -soft_masking true \
  -evalue 700 -searchsp 1750000000000 \
  -outfmt 6 \
  > ${outdir}/$bam_prefix.contaminant.blastout

cat ${outdir}/$bam_prefix.contaminant.blastout | \
  grep 'NGB0097' | \
    awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" \
    '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | \
    sort -u \
    > ${outdir}/$bam_prefix.blocklist &

filterbyname.sh -Xmx$mem ignorebadquality in=$fastq out=$filt_fastq names=${outdir}/$bam_prefix.blocklist

rm $fasta
rm $fastq
