#!/bin/bash
set -euo pipefail
module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
BATCH=sims
RESULTSDIR=$BASEDIR/scratch/results
TRUTH_DIR=$RESULTSDIR/$BATCH.truth
PARALLEL=40

function make_truth_mutrels {
  mkdir -p $TRUTH_DIR
  for paramsfn in $PAIRTREE_INPUTS_DIR/*.params.json; do
    runid=$(basename $paramsfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/make_truth_mutrel.py" \
      "$PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$TRUTH_DIR/$runid.mutrel.npz"
  done | parallel -j$PARALLEL --halt 1
}

function run_eval {
  cd $RESULTSDIR
  prefix=steph.xeno.nogarb
  SCORESDIR=$RESULTSDIR/$prefix/scores
  mkdir -p $SCORESDIR

  for mutrelfn in ${prefix}.truth/*.mutrel.npz; do
    runid=$(basename $mutrelfn | cut -d. -f1)
    mutrels="truth=${prefix}.truth/$runid.mutrel.npz "
    mutrels+="pwgs=${prefix}.pwgs/$runid/$runid.mutrel.npz "
    mutrels+="pairtree_trees=${prefix}.pairtree/$runid.trees.mutrel.npz "
    mutrels+="pairtree_pairs=${prefix}.pairtree/$runid.pairs.mutrel.npz"

    for M in $(echo $mutrels | tr ' ' '\n' | cut -d= -f2); do
      [[ ! -f $M ]] && continue 2
    done

    echo "cd $RESULTSDIR && python3 $SCRIPTDIR/eval.py $mutrels > $SCORESDIR/$runid.score.txt"
  done | parallel -j$PARALLEL --halt 0

  (
    cd $SCORESDIR
    echo 'runid,'$(head -n1 $(ls *.score.txt | head -n1))
    for foo in *.score.txt; do
      echo $(echo $foo | cut -d. -f1),$(tail -n+2 $foo) 
    done
  ) | grep -v ',$' | curl -F c=@- https://ptpb.pw
}

function main {
  make_truth_mutrels
  run_eval
}

main
