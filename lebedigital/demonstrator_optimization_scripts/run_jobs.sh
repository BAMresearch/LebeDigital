#!/bin/bash
#SBATCH --job-name=LBD_optimization
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --partition=batch_SKL,batch_SNB
#SBATCH --array=1-100
#SBATCH --output=slurm-%A_%a.out

# 1-$1, the $1 is for the first argument of the sbatch run_jobs 100, where 100 samples total
# Load Python Module
source /home/atul/.bashrc

. /home/atul/miniconda3/etc/profile.d/conda.sh

conda activate lebedigital

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
#echo "This is array task ${SLURM_ARRAY_TASK_ID}" >> output.txt

# creates multiple folders for each sample
#python farm_workflow.py

# Run the workflow in folder parallely
python parallel_compute_workflow.py $SLURM_ARRAY_TASK_ID

echo "!!!! Job array completed.!!!!"