#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=def-mzhen
#SBATCH --mem=15G
#SBATCH --array=1-10
#SBATCH --cpus-per-task=12
#SBATCH --job-name=2
#SBATCH --output=../../../results/log/consensus_ds%x-k%a-%j.out

module load matlab
export NUM_THREAD=$SLURM_CPUS_PER_TASK
export DS_IDX=$SLURM_JOB_NAME
export NUM_CLUSTERS=$SLURM_ARRAY_TASK_ID
export FOLDER_NAME='B-all'
echo bash prepare to run ds${DS_IDX} k ${NUM_CLUSTERS} using ${NUM_THREAD} cores
matlab -batch 'run("run_fitConsensus.m");exit;'
echo bash finished running matlab script

