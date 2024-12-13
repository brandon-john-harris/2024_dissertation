#!/bin/bash                                                                                                                                                      
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=TopoBuild # Job name
#SBATCH --array=1-500 # for this script adjust to match number of array directories
#SBATCH --nodes=1 # single node run
#SBATCH --mem=10000
#SBATCH --time=1-00:00:00 # time to run in days/hrs/min/sec
##SBATCH --output=slurm-out-%A_%a.out # output file
##SBATCH --error=slurm-err-%A_%a.err # error file

echo "My SLURM_ARRAY_TASK_ID:  $SLURM_ARRAY_TASK_ID"
number="$(printf "%04d" $SLURM_ARRAY_TASK_ID)"
mkdir $number
cp -r ./design_topo.sh /${PWD}/${number}/
cd /${PWD}/${number}/
array_suffix=$number
export array=$array_suffix
chmod +x design_topo.sh
./design_topo.sh