#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree

BATCH=sims.smallalpha.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.multichain
TRUTH_DIR=$BASEDIR/scratch/results/sims.smallalpha.truth
SCORESDIR=$BASEDIR/scratch/scores

source $SCRIPTDIR/util.sh

function compute_entropy {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    is_run_complete $resultsfn || continue

    cmd="cd $outdir && "
    cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/compute_entropy.py "
    cmd+="--truth $TRUTH_DIR/$runid/$runid.results.npz "
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
  echo $cmd | bash
}

function main {
  compute_entropy | grep -v -e K30_ -e K100_ | parallel -j80 --halt 1 --eta
  compute_entropy | grep    -e K30_ -e K100_ | parallel -j10 --halt 1 --eta
  combine
  plot
}

main
