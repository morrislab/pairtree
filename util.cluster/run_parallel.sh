#!/bin/bash
# Usage: sbatch -N <num_nodes> run_parallel.sh jobs.txt jobs.log

#SBATCH --ntasks-per-node=40
#SBATCH --time=23:58:00
#SBATCH --job-name lolparallel
#SBATCH --output /home/q/qmorris/jawinter/scratch/tmp/slurm-%j.stdout
#SBATCH --error /home/q/qmorris/jawinter/scratch/tmp/slurm-%j.stderr

PARALLEL=$SLURM_NTASKS_PER_NODE
PARALLEL=10
jobsfn=$1
logfn=$2

module load gnu-parallel
  
# DIRECTORY TO RUN - $SLURM_SUBMIT_DIR is the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR
HOSTS=$(scontrol show hostnames $SLURM_NODELIST | tr '\n' ,)
cat $jobsfn | parallel --halt 1 --resume --joblog $logfn -j $PARALLEL -S $HOSTS --wd $PWD
