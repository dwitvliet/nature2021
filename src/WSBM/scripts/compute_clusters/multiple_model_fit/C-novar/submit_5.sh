#!/bin/bash
#SBATCH -n 64
#SBATCH -t 6-00:00
#SBATCH -p serial_requeue
#SBATCH --mem=60000
#SBATCH --array=10
#SBATCH --job-name=C_ds5
#SBATCH -o ../../../results/log/%x-k%a-%j.out
#SBATCH -e ../../../results/log/%x-k%a-%j.err

module load matlab/R2019a-fasrc01
export NUM_THREAD=64
export DS_IDX=5
export NUM_CLUSTERS=$SLURM_ARRAY_TASK_ID
export FOLDER_NAME='C-no-variable'
echo bash prepare to run ds${DS_IDX} k ${NUM_CLUSTERS} using ${NUM_THREAD} cores
matlab -batch 'run("wsbm_consensus_batch.m");exit;'
echo bash finished running matlab script

