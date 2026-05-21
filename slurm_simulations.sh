#!/bin/bash
#SBATCH --partition=main
#SBATCH --no-requeue
#SBATCH --job-name=acic-test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=00:30:00
#SBATCH --array=1-400                 
#SBATCH --output=logs/slurm.%N.%j.%a.out
#SBATCH --error=logs/slurm.%N.%j.%a.err


cd /scratch/am3923/second_order_correction

mkdir -p logs results

export PATH=/home/am3923/R-4.4.2/bin:$PATH

PARAM_NUM=$(( (SLURM_ARRAY_TASK_ID - 1) / 50 + 17 ))
SIM_NUM=$(( (SLURM_ARRAY_TASK_ID - 1) % 50 + 1 ))
RUN_ID=${1:?Error: RUN_ID required. Usage: sbatch script.sh <run_id>}
echo "Task $SLURM_ARRAY_TASK_ID → param=$PARAM_NUM sim=$SIM_NUM run_id=$RUN_ID"

module load gcc/5.4

Rscript design_5_acic_data_1.R $PARAM_NUM $SIM_NUM $RUN_ID
