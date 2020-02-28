#!/bin/bash
function para {
  parallel -j$PARALLEL --halt 1 --eta
}
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
EXPTSDIR=~/work/pairtree-experiments
PARALLEL=80

function para {
  parallel -j$PARALLEL --halt 1 --eta
}

function make_mle_mutphis {
  mkdir -p $MLE_MUTPHIS_DIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)

    echo "python3 $SCRIPTDIR/make_mle_mutphis.py" \
      "$ssmfn" \
      "$MLE_MUTPHIS_DIR/$runid.mutphi.npz"
  done | parallel -j80 --halt 1 --eta
}

function make_baseline_phi {
  for resultfn in $BASELINE_DIR/*/*.results.npz; do
    runid=$(basename $resultfn | cut -d. -f1)
    resultdir=$(dirname $resultfn)

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_baseline_mutdist_from_pairtree_run.py "
    cmd+="$resultfn "
    cmd+="$resultdir/${runid}.phi.npz"
    echo $cmd
  done
}

function make_results_paths {
  runid=$1
  result_type=$2
  paths=""

  if [[ $result_type == mutphi ]]; then
    paths+="mle_unconstrained=${BATCH}.mle_unconstrained/$runid.$result_type.npz "
  fi

  paths+="baseline=${BASELINE_DIR}/$runid/${runid}.${result_type}.npz "
  paths+="pairtree=${BATCH}.xeno.pairtree.multichain.rprop/${runid}/${runid}.${result_type}.npz "
  paths+="pairtree_proj=${BATCH}.xeno.pairtree.multichain.projection/${runid}/${runid}.${result_type}.npz "
  paths+="pairtree_hbstruct=${BASELINE_DIR}/${runid}/${runid}.${result_type}.npz "
  paths+="pwgs_supervars=${BATCH}.pwgs.supervars/$runid/${runid}.${result_type}.npz "
  paths+="lichee=${BATCH}.xeno.lichee/$runid/$runid.$result_type.npz "

  echo $paths
}

function eval_mutphis {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for mutphifn in $(ls $BASELINE_DIR/*/*.mutphi.npz | sort --random-sort); do
    runid=$(basename $mutphifn | cut -d. -f1)
    mutphis=$(make_results_paths $runid mutphi)

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutphis.py " \
      "--params $PAIRTREE_INPUTS_DIR/${runid}.nostruct.params.json" \
      "$mutphis " \
      "> $SCORESDIR/$BATCH/$runid.mutphi_score.txt"
  done
}

function eval_mutdists {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for phifn in $(ls $BASELINE_DIR/*/*.phi.npz | sort --random-sort); do
    # We don't actually use `phifn`, but that's okay -- we just want a
    # reference for what runs exist.
    runid=$(basename $phifn | cut -d. -f1)
    mutdists=$(make_results_paths $runid mutdist)

    for P in 1 2; do
      cmd="cd $RESULTSDIR && "
      cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutdists.py "
      cmd+="-p $P "
      cmd+="--params $PAIRTREE_INPUTS_DIR/${runid}.nostruct.params.json "
      cmd+="$mutdists "
      cmd+="> $SCORESDIR/$BATCH/$runid.mutdistl${P}_score.txt "
      echo $cmd
    done
  done
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
        echo "Methods in $foo ($(head -n1 $foo)) don't match expected $methods" >&2
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
  #cat $outfn | curl -F c=@- https://ptpb.pw >&2
}

function plot_single_vs_others {
  single_method=$1
  score_type=$2
  infn=$3
  outfn=$4

  cmd="python3 $SCRIPTDIR/plotter/plot_single_vs_others.py "
  cmd+="--template seaborn "
  cmd+="--score-type $score_type "
  if [[ $score_type == mutphi ]]; then
    cmd+="--baseline baseline "
  fi
  cmd+="$infn $method $outfn"
  echo $cmd
}

function plot_results_steph {
  (
    eval_mutphis
    eval_mutdists
  ) | para
  for type in mutphi mutdistl1 mutdistl2; do
    compile_scores $type
  done

  for score_type in mutphi mutdistl1 mutdistl2; do
    cmd="python3 $SCRIPTDIR/plotter/plot_steph.py "
    cmd+="--template seaborn "
    cmd+="--score-type $score_type "
    cmd+="--hide-methods pairtree_hbstruct,pairtree_proj "
    if [[ $score_type == mutphi ]]; then
      cmd+="--baseline baseline "
    fi
    cmd+="$SCORESDIR/${BATCH}.${score_type}.{txt,html}"
    echo $cmd
  done | para
}

function make_index {
  cd $SCORESDIR
  (
    for anal in mutdistl1 mutdistl2 mutphi; do
      echo "<h3>$anal</h3><ul>"
      for foo in $BATCH.*$anal*.html; do
        echo "<li><a href=$foo>$foo</a></li>"
      done
      echo "</ul>"
    done
  ) > index.${BATCH}.html
}

function main {
  export BATCH=steph
  export RESULTSDIR=$BASEDIR/scratch/results
  export SCORESDIR=$BASEDIR/scratch/scores
  export PAIRTREE_INPUTS_DIR=$EXPTSDIR/inputs/steph.xeno.pairtree
  export BASELINE_DIR=$RESULTSDIR/steph.xeno.pairtree.hbstruct.projection
  export MLE_MUTPHIS_DIR=$RESULTSDIR/${BATCH}.mle_unconstrained

  #make_mle_mutphis | para
  #make_baseline_phi | para
  plot_results_steph
  make_index
}

main
