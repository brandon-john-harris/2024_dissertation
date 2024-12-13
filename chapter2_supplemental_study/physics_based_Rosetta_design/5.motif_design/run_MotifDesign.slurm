#!/bin/bash
#SBATCH --partition=production # partition to submit to 
#SBATCH --job-name=MotifDesign # Job name                
#SBATCH --array=1-30 # for this script adjust to match number of array directories
#SBATCH --nodes=1 # single node run                                                                       
#SBATCH --time=0-12:00:00 # time to run in days/hrs/min/sec       
#SBATCH --mem=10000                                        
#SBATCH --output=slurm-out-%A_%a.out # output file                                                        
#SBATCH --error=slurm-err-%A_%a.err # error file

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
array_suffix=$SLURM_ARRAY_TASK_ID
export array_suffix=$array_suffix
chmod +x MotifDesign.sh