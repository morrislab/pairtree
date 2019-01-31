#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
BATCH=sims
RESULTSDIR=$BASEDIR/scratch/results
SCORESDIR=$BASEDIR/scratch/scores
TRUTH_DIR=$RESULTSDIR/$BATCH.truth
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
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

function make_results_paths {
  runid=$1
  result_type=$2

  paths="truth=${BATCH}.truth/$runid.$result_type.npz "
  paths+="pairtree_trees_llh=${BATCH}.pairtree.fixedclusters/$runid.pairtree_trees_llh.$result_type.npz "
  paths+="pairtree_trees_uniform=${BATCH}.pairtree.fixedclusters/$runid.pairtree_trees_uniform.$result_type.npz "
  paths+="pairtree_clustrel=${BATCH}.pairtree.fixedclusters/$runid.pairtree_clustrel.$result_type.npz "
  paths+="pwgs_trees_single_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
  paths+="pwgs_trees_single_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_uniform.$result_type.npz "
  paths+="pwgs_trees_multi_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
  paths+="pwgs_trees_multi_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_uniform.$result_type.npz "
  paths+="pastri_trees_llh=${BATCH}.pastri.informative/$runid.pastri_trees_llh.$result_type.npz "
  paths+="pastri_trees_uniform=${BATCH}.pastri.informative/$runid.pastri_trees_uniform.$result_type.npz"

  echo $paths
}

function eval_mutrels {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for mutrelfn in $(ls $TRUTH_DIR/*.mutrel.npz | sort --random-sort); do
    runid=$(basename $mutrelfn | cut -d. -f1)
    mutrels=$(make_results_paths $runid mutrel)

    #for M in $(echo $mutrels | tr ' ' '\n' | cut -d= -f2); do
    #  [[ -f $M ]] || continue 2
    #done

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutrels.py " \
      "--discard-garbage" \
      "--params $PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$mutrels " \
      "> $SCORESDIR/$BATCH/$runid.mutrel_score.txt"
  done #| parallel -j$PARALLEL --halt 1
}

function eval_mutphis {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for mutphifn in $(ls $TRUTH_DIR/*.mutphi.npz | sort --random-sort); do
    runid=$(basename $mutphifn | cut -d. -f1)
    mutphis=$(make_results_paths $runid mutphi)

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutphis.py " \
      "--params $PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$mutphis " \
      "> $SCORESDIR/$BATCH/$runid.mutphi_score.txt"
  done #| parallel -j$PARALLEL --halt 1
}

function compile_scores {
  score_type=$1
  suffix=${score_type}_score.txt
  outfn=$SCORESDIR/$BATCH.$score_type.txt
  (
    cd $SCORESDIR/$BATCH
    methods=$(head -n1 $(ls *.$suffix | head -n1))
    echo 'runid,'$methods
    for foo in *.$suffix; do
      if [[ $(head -n1 $foo) != $methods ]]; then
        echo "Methods in $foo don't match expected $methods" >&2
        exit 1
      fi
      echo $(echo $foo | cut -d. -f1),$(tail -n+2 $foo) 
    done
  ) > $outfn
  cat $outfn | curl -F c=@- https://ptpb.pw >&2
}

function run_server {
  cd $SCRIPTDIR
  MUTREL_RESULTS=$SCORESDIR/sims.mutrel.txt MUTPHI_RESULTS=$SCORESDIR/sims.mutphi.txt gunicorn -w 4 -b 0.0.0.0:8000 visualize:server
}

function main {
  #make_truth_mutrels
  #eval_mutrels
  #compile_scores mutrel
  #eval_mutphis
  #compile_scores mutphi
  run_server
}

main
