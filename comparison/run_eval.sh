#!/bin/bash
set -euo pipefail
#module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
BATCH=sims
RESULTSDIR=$BASEDIR/scratch/results
TRUTH_DIR=$RESULTSDIR/$BATCH.truth
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PARALLEL=40

function make_truth_mutrels {
  mkdir -p $TRUTH_DIR
  for datafn in $PAIRTREE_INPUTS_DIR/*.data.pickle; do
    runid=$(basename $datafn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/make_truth_mutrel.py" \
      "$datafn" \
      "$TRUTH_DIR/$runid.mutrel.npz"
  done | parallel -j$PARALLEL --halt 1
}

function run_eval {
  cd $RESULTSDIR
  prefix=sims
  outfn=$RESULTSDIR/$prefix.txt
  SCORESDIR=$RESULTSDIR/$prefix/scores
  mkdir -p $SCORESDIR
  rm -f $SCORESDIR/*.score.txt

  for mutrelfn in $TRUTH_DIR/*.mutrel.npz; do
    runid=$(basename $mutrelfn | cut -d. -f1)
    mutrels="truth=${prefix}.truth/$runid.mutrel.npz "
    mutrels+="pwgs_trees_single_llh=${prefix}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.mutrel.npz "
    mutrels+="pwgs_trees_single_uniform=${prefix}.pwgs.supervars/$runid/$runid.pwgs_trees_single_uniform.mutrel.npz "
    mutrels+="pwgs_trees_multi_llh=${prefix}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.mutrel.npz "
    mutrels+="pwgs_trees_multi_uniform=${prefix}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_uniform.mutrel.npz "
    mutrels+="pairtree_trees_llh=${prefix}.pairtree.fixedclusters/$runid.pairtree_trees_llh.mutrel.npz "
    mutrels+="pairtree_trees_uniform=${prefix}.pairtree.fixedclusters/$runid.pairtree_trees_uniform.mutrel.npz "
    mutrels+="pairtree_clustrel=${prefix}.pairtree.fixedclusters/$runid.pairtree_clustrel.mutrel.npz "
    mutrels+="pastri_trees_llh=${prefix}.pastri.informative/$runid.pastri_trees_llh.mutrel.npz "
    mutrels+="pastri_trees_uniform=${prefix}.pastri.informative/$runid.pastri_trees_uniform.mutrel.npz"

    #for M in $(echo $mutrels | tr ' ' '\n' | cut -d= -f2); do
    #  [[ ! -f $M ]] && continue 2
    #done

    echo "cd $RESULTSDIR && python3 $SCRIPTDIR/eval.py $mutrels > $SCORESDIR/$runid.score.txt"
  done | parallel -j$PARALLEL

  (
    cd $SCORESDIR
    echo 'runid,'$(head -n1 $(ls *.score.txt | head -n1))
    for foo in *.score.txt; do
      echo $(echo $foo | cut -d. -f1),$(tail -n+2 $foo) 
    done
  ) > $outfn
  cat $outfn | curl -F c=@- https://ptpb.pw
}

function main {
  #make_truth_mutrels
  run_eval
}

main
