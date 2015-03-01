#!/bin/bash 

#----------------------------------------------------
# Generic SLURM script 
#----------------------------------------------------

#SBATCH -J scan-motif # Job name
#SBATCH -o gls.%j.out # stdout; %j expands to jobid
#SBATCH -e gls.%j.err # stderr; skip to combine stdout and stderr
#SBATCH -p normal # queue
#SBATCH -N 1 -n 8 # one node and one task
#SBATCH -t 48:00:00 # max time
#SBATCH --mail-user=zw355@cornell.edu
#SBATCH --mail-type=ALL

$R --vanilla --quiet  < gls-test-plink.r > gls-test-plink.out 
