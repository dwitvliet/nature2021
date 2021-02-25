#!/bin/bash
#SBATCH -n 64
#SBATCH -t 3-00:00
#SBATCH -p serial_requeue
#SBATCH --mem=60000
#SBATCH --array=9
#SBATCH --job-name=D_ds2
#SBATCH -o ../../../results/log/%x-k%a-%j.out
#SBATCH -e ../../../results/log/%x-k%a-%j.err

module load matlab/R2019a-fasrc01
export NUM_THREAD=64
export DS_IDX=2
export NUM_CLUSTERS=$SLURM_ARRAY_TASK_ID
export FOLDER_NAME='D-only_stable'
echo bash prepare to run ds${DS_IDX} k ${NUM_CLUSTERS} using ${NUM_THREAD} cores
matlab -batch 'run("wsbm_consensus_batch.m");exit;'
echo bash finished running matlab script

