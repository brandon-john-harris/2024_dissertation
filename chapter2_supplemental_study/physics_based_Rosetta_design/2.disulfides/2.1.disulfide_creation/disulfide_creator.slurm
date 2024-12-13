#!/bin/bash                                                                                                                                                      
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=DsfCreate # Job name
#SBATCH --array=1 # for this script adjust to match number of array directories
#SBATCH --nodes=1 # single node run
#SBATCH --mem=1210
#SBATCH --time=0-01:00:00 #time to run in days/hrs/min/sec
##SBATCH --output=slurm-out-%A_%a.out # output file
##SBATCH --error=slurm-err-%A_%a.err # error file

### Updated 2021-04-02 by bjh. Adjusted execcutable options. Not Time/Memory optimized ###

export line=$line
chmod +x disulfide_creator.sh
./disulfide_creator.sh