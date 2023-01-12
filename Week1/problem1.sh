#!/bin/bash
#SBATCH --job-name=problem1    ## Name of the job.
#SBATCH -A etom2       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --cpus-per-task=1  ## number of cores the job needs

wget https://wfitch.bio.uci.edu/~tdlong/problem1.tar.gz
tar -xvf problem1.tar.gz
rm problem1.tar.gz

cd problem1
touch answers.txt

head -10 p.txt | tail -1 > answers.txt  
head -10 f.txt | tail -1 >> answers.txt

