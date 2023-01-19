#!/bin/bash
#SBATCH --job-name=ATACseq    ## Name of the job.
#SBATCH -A etom2       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --cpus-per-task=1  ## number of cores the job needs

SourceDir="/data/class/ecoevo283/public/RAWDATA/ATACseq"
DestDir="/pub/etom2/EE283/EE283/Week2/ATACseq/raw"
File="ATACseq.labels.txt"
while read line
do
   echo "$line"
   barcode=$(echo $line | cut -f1 -d" ")
   genotype=$(echo $line | cut -f2 -d" ")
   tissue=$(echo $line | cut -f3 -d" ")
   bioRep=$(echo $line | cut -f4 -d" ")
   READ1=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R1.fq.gz")
   READ2=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R2.fq.gz") 
   ln -s $READ1 $DestDir/${genotype}_${tissue}_${bioRep}_R1.fq.gz
   ln -s $READ2 $DestDir/${genotype}_${tissue}_${bioRep}_R2.fq.gz

done < $File

