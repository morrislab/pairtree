#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
BATCH=sims
RESULTSDIR=$BASEDIR/scratch/results
SCORESDIR=$BASEDIR/scratch/scores
TRUTH_DIR=$RESULTSDIR/${BATCH}.truth
MLE_MUTPHIS_DIR=$RESULTSDIR/${BATCH}.mle_unconstrained
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PARALLEL=20

function make_truth {
  mkdir -p $TRUTH_DIR

  for datafn in $PAIRTREE_INPUTS_DIR/*.data.pickle; do
    runid=$(basename $datafn | cut -d. -f1)
      #"--enumerate-trees" \
    echo "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutrel.py" \
      "$datafn" \
      "$TRUTH_DIR/$runid.mutrel.npz"
    echo "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutphi.py" \
      "$datafn" \
      "$TRUTH_DIR/$runid.mutphi.npz"
  done #| parallel -j$PARALLEL --halt 1
}

function make_mle_mutphis {
  mkdir -p $MLE_MUTPHIS_DIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    [[ $runid =~ K30 || $runid =~ K100 ]] && continue

    echo "python3 $SCRIPTDIR/calc_mle_mutphis.py" \
      "$ssmfn" \
      "$MLE_MUTPHIS_DIR/$runid.mutphi.npz"
  done #| parallel -j$PARALLEL --halt 1
}

function make_results_paths {
  runid=$1
  result_type=$2
  paths=""

  if [[ $result_type == mutphi ]]; then
    paths+="mle_unconstrained=${BATCH}.mle_unconstrained/$runid.$result_type.npz "
  fi

  if [[ $BATCH == steph ]]; then
    paths+="truth=${BATCH}.truth/$runid.$result_type.npz "
    paths+="pairtree_trees_llh=${BATCH}.xeno.withgarb.pairtree/$runid.pairtree_trees_llh.$result_type.npz "
    paths+="pairtree_trees_uniform=${BATCH}.xeno.withgarb.pairtree/$runid.pairtree_trees_uniform.$result_type.npz "
    paths+="pairtree_handbuilt=${BATCH}.pairtree.hbstruct/$runid.$result_type.npz "
    paths+="pwgs_allvars_single_uniform=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_single_uniform.$result_type.npz "
    paths+="pwgs_allvars_single_llh=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    paths+="pwgs_allvars_multi_uniform=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_multi_uniform.$result_type.npz "
    paths+="pwgs_allvars_multi_llh=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
    paths+="pwgs_supervars_single_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_uniform.$result_type.npz "
    paths+="pwgs_supervars_single_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    paths+="pwgs_supervars_multi_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_uniform.$result_type.npz "
    paths+="pwgs_supervars_multi_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz"
  elif [[ $BATCH == sims ]]; then
    paths+="truth=${TRUTH_DIR}/$runid.$result_type.npz "
    paths+="pairtree_trees_llh=${BATCH}.pairtree.fixedclusters/$runid.pairtree_trees_llh.$result_type.npz "
    paths+="pairtree_trees_uniform=${BATCH}.pairtree.fixedclusters/$runid.pairtree_trees_uniform.$result_type.npz "
    paths+="pairtree_clustrel=${BATCH}.pairtree.fixedclusters/$runid.pairtree_clustrel.$result_type.npz "
    paths+="pwgs_trees_single_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    paths+="pwgs_trees_single_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_uniform.$result_type.npz "
    paths+="pwgs_trees_multi_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
    paths+="pwgs_trees_multi_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_uniform.$result_type.npz "
    paths+="pastri_trees_llh=${BATCH}.pastri.informative/$runid.pastri_trees_llh.$result_type.npz "
    paths+="pastri_trees_uniform=${BATCH}.pastri.informative/$runid.pastri_trees_uniform.$result_type.npz"
  fi

  echo $paths
}

function eval_mutrels {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for mutrelfn in $(ls $TRUTH_DIR/*.mutrel.npz | sort --random-sort); do
    runid=$(basename $mutrelfn | cut -d. -f1)
    mutrels=$(make_results_paths $runid mutrel)

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

  for mutphifn in $(ls $MLE_MUTPHIS_DIR/*.mutphi.npz | sort --random-sort); do
    runid=$(basename $mutphifn | cut -d. -f1)
    mutphis=$(make_results_paths $runid mutphi)

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutphis.py " \
      "--params $PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "--ssms $PAIRTREE_INPUTS_DIR/$runid.ssm" \
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
      S=$(echo $foo | cut -d. -f1)

      # These are the runs we did not use in the paper.
      for bad_sampid in SJETV010{,nohypermut,stephR1,stephR2} SJBALL022610; do
        if [[ $S == $bad_sampid ]]; then
          continue 2
        fi
      done
      echo $S,$(tail -n+2 $foo)
    done
  ) > $outfn
  cat $outfn | curl -F c=@- https://ptpb.pw >&2
}

function plot_comparison {
  # To redirect port 80 to 8000:
  #   iptables -A PREROUTING -t nat -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 8000
  cd $SCRIPTDIR

  export MUTREL_RESULTS=$SCORESDIR/$BATCH.mutrel.txt
  export MUTPHI_RESULTS=$SCORESDIR/$BATCH.mutphi.txt
  export SIM_PARAMS="GKMST"
  production=false

  if [[ production == true ]]; then
    gunicorn -w 4 -b 0.0.0.0:8000 plot_comparison:server
  else
    python3 plot_comparison.py
  fi
}

function plot_individual {
  cd $SCORESDIR
  for batch in steph sims; do
    for ptype in mutphi mutrel; do
      echo "python3 ~/work/pairtree/comparison/plot_individual.py" \
        "--template plotly_white" \
        "$( [[ $batch == sims && $ptype == mutphi ]] && echo --max-y 10)" \
        "$( [[ $batch == sims ]] && echo --partition-by-samples)" \
        "$batch.$ptype.txt" \
        "$batch.$ptype.html"
    done
  done | parallel -j40 --halt 1
}

function main {
  #make_truth
  #make_mle_mutphis

  eval_mutrels
  eval_mutphis
  #compile_scores mutrel
  #compile_scores mutphi

  #plot_individual
  #plot_comparison
}

main
