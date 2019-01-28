#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
BATCH=sims
RESULTSDIR=$BASEDIR/scratch/results
TRUTH_DIR=$RESULTSDIR/$BATCH.truth
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PREFIX=sims
SCORESDIR=$RESULTSDIR/$PREFIX/scores
PARALLEL=20

function make_truth_mutrels {
  mkdir -p $TRUTH_DIR
  for datafn in $PAIRTREE_INPUTS_DIR/*.data.pickle; do
    runid=$(basename $datafn | cut -d. -f1)
    echo "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutrel.py" \
      "$datafn" \
      "$TRUTH_DIR/$runid.mutrel.npz" \
      "$TRUTH_DIR/$runid.mutphi.npz"
  done #| parallel -j$PARALLEL --halt 1
}

function run_eval {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR
  #rm -f $SCORESDIR/*.score.txt

  for mutrelfn in $(ls $TRUTH_DIR/*.mutrel.npz | sort --random-sort); do
    runid=$(basename $mutrelfn | cut -d. -f1)
    mutrels="truth=${PREFIX}.truth/$runid.mutrel.npz "
    mutrels+="pairtree_trees_llh=${PREFIX}.pairtree.fixedclusters/$runid.pairtree_trees_llh.mutrel.npz "
    mutrels+="pairtree_trees_uniform=${PREFIX}.pairtree.fixedclusters/$runid.pairtree_trees_uniform.mutrel.npz "
    mutrels+="pairtree_clustrel=${PREFIX}.pairtree.fixedclusters/$runid.pairtree_clustrel.mutrel.npz "
    #mutrels+="pwgs_trees_single_llh=${PREFIX}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.mutrel.npz "
    #mutrels+="pwgs_trees_single_uniform=${PREFIX}.pwgs.supervars/$runid/$runid.pwgs_trees_single_uniform.mutrel.npz "
    #mutrels+="pwgs_trees_multi_llh=${PREFIX}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.mutrel.npz "
    #mutrels+="pwgs_trees_multi_uniform=${PREFIX}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_uniform.mutrel.npz "
    #mutrels+="pastri_trees_llh=${PREFIX}.pastri.informative/$runid.pastri_trees_llh.mutrel.npz "
    #mutrels+="pastri_trees_uniform=${PREFIX}.pastri.informative/$runid.pastri_trees_uniform.mutrel.npz"

    #for M in $(echo $mutrels | tr ' ' '\n' | cut -d= -f2); do
    #  [[ -f $M ]] || continue 2
    #done

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval.py " \
      "--discard-garbage" \
      "--params $PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$mutrels " \
      "> $SCORESDIR/$runid.score.txt"
  done #| parallel -j$PARALLEL --halt 1
}

function compile_scores {
  outfn=$RESULTSDIR/$PREFIX.txt
  (
    cd $SCORESDIR
    echo 'runid,'$(head -n1 $(ls *.score.txt | head -n1))
    for foo in *.score.txt; do
      echo $(echo $foo | cut -d. -f1),$(tail -n+2 $foo) 
    done
  ) > $outfn
  cat $outfn | curl -F c=@- https://ptpb.pw >&2
}

function main {
  #make_truth_mutrels
  #run_eval
  compile_scores
}

main
