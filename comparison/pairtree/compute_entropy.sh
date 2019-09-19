#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree

BATCH=sims.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.lol69
TRUTH_DIR=$BASEDIR/scratch/results/sims.truth

PARALLEL=80

source $SCRIPTDIR/util.sh

function compute_entropy {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    is_run_complete $resultsfn || continue

    cmd="cd $outdir && "
    cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/compute_entropy.py "
    if ! is_big_K $runid; then
      cmd+="--truth $TRUTH_DIR/$runid.results.npz "
    else
      continue
    fi
    cmd+="$resultsfn "
    cmd+="> $outdir/$runid.entropy.txt"
    echo $cmd
  done 
}

function main {
  compute_entropy | parallel -j$PARALLEL --halt 1 --eta
}

main
