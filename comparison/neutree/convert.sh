#!/bin/bash
set -euo pipefail
shopt -s nullglob

BASEDIR=~/work/pairtree
RESULTSDIR=$BASEDIR/scratch/results
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.pairtree
TRUTH_DIR=$RESULTSDIR/sims.smallalpha.truth
NEUTREEDIR=$BASEDIR/comparison/neutree
PARALLEL=40

function create_evals {
  results=$1
  nresults=$(ls $RESULTSDIR/$results/*/*.neutree.pickle | wc -l)
  echo -e "$results\t$nresults" >&2

  for neutreefn in $RESULTSDIR/$results/*/*.neutree.pickle; do
    outdir=$(dirname $neutreefn)
    runid=$(basename $neutreefn | cut -d. -f1)
    basepath="${outdir}/${runid}"

    (
      cmd="cd $outdir && "
      cmd+="OMP_NUM_THREADS=1 python3 $NEUTREEDIR/make_mutphis.py "
      cmd+="$neutreefn "
      cmd+="${PAIRTREE_INPUTS_DIR}/${runid}.ssm "
      cmd+="${basepath}.mutphi.npz "
      echo $cmd

      cmd="cd $outdir && "
      cmd+="OMP_NUM_THREADS=1 python3 $NEUTREEDIR/make_mutdists.py "
      cmd+="$neutreefn "
      cmd+="${TRUTH_DIR}/${runid}/${runid}.phi.npz "
      cmd+="${basepath}.mutdist.npz "
      echo $cmd

      cmd="cd $outdir && "
      cmd+="OMP_NUM_THREADS=1 python3 $NEUTREEDIR/make_mutrels.py "
      cmd+="$neutreefn "
      cmd+="${basepath}.mutrel.npz "
      echo $cmd
    )
  done
}

function main {
  results=""
  results+="sims.smallalpha.citup.rawvars.qip "
  results+="sims.smallalpha.lichee "
  results+="sims.smallalpha.pairtree.multichain "
  results+="sims.smallalpha.pairtree.rprop "
  results+="sims.smallalpha.pastri "
  results+="sims.smallalpha.pwgs.supervars "

  cmds=$(for result in $results; do
    create_evals $result
  done)

  echo "$cmds" | grep mutphi                   | sort --random-sort | parallel -j$PARALLEL --halt 1 --eta
  echo "$cmds" | grep mutdist | grep -v lichee | sort --random-sort | parallel -j$PARALLEL --halt 1 --eta
  echo "$cmds" | grep mutphi                   | sort --random-sort | parallel -j5         --halt 1 --eta
}

main
