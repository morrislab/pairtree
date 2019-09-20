#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree

BATCH=sims.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.lol69
TRUTH_DIR=$BASEDIR/scratch/results/sims.truth
SCORESDIR=$BASEDIR/scratch/scores

PARALLEL=10

source $SCRIPTDIR/util.sh

function compute_entropy {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    is_run_complete $resultsfn || continue

    cmd="cd $outdir && "
    cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/compute_entropy.py "
    cmd+="--truth $TRUTH_DIR/$runid.results.npz "
    cmd+="$resultsfn "
    cmd+="> $outdir/$runid.entropy.txt"
    echo $cmd
  done 
}

function combine {
  header=$(head -n1 $(ls $PAIRTREE_RESULTS_DIR/*/*.entropy.txt | head -n1))
  (
    echo 'runid',$header
    for entfn in $PAIRTREE_RESULTS_DIR/*/*.entropy.txt; do
      runid=$(basename $entfn | cut -d. -f1)
      if [[ $(head -n1 $entfn) != $header ]]; then
        echo "Header in $entfn ($(head -n1 $entfn)) doesn't match expected $header" >&2
        exit 1
      fi
      echo $runid,$(tail -n+2 $entfn)
    done
  ) > $SCORESDIR/entropy.txt
}

function plot {
  cmd="python3 $SCRIPTDIR/plot_entropy.py "
  cmd+="$SCORESDIR/entropy.txt "
  cmd+="$SCORESDIR/entropy.html"
  echo $cmd
}

function main {
  #compute_entropy | parallel -j$PARALLEL --halt 1 --eta
  combine
  plot
}

main
