#!/bin/sh
set -euo pipefail
shopt -s nullglob

BASEDIR=~/work/pairtree
PEARSIMDIR=~/work/pearsim
JOBDIR=~/jobs
PYTHON=python3

K=30
S=30
T=1000
M_per_cluster=20
G_per_cluster=2
ALPHA=0.1
PARA=80
INDIR=$BASEDIR/scratch/inputs/missedcna
RESULTSDIR=$BASEDIR/scratch/results/missedcna


function commafy {
  echo $(echo $1 | sed 's/\./,/g')
}

function make_inputs {
  mkdir -p $INDIR && cd $INDIR

  for K in 10 30; do
  for MIN_GARB_PHI_DELTA in 0.0005 0.001 0.01 0.05 0.1; do
  for run in $(seq 100); do
    M=$(echo "$K * $M_per_cluster" | bc)
    G=$(echo "$K * $G_per_cluster" | bc)
    jobname="sim_missedcna_K${K}_mindelta$(commafy $MIN_GARB_PHI_DELTA)_run${run}"

    echo "$PYTHON $PEARSIMDIR/make_simulated_data.py" \
      "--seed $run" \
      "--write-clusters" \
      "-K $K" \
      "-S $S" \
      "-T $T" \
      "-M $M" \
      "-G $G" \
      "--alpha $ALPHA" \
      "--garbage-type obvs_missed_cna" \
      "--min-garb-phi-delta $MIN_GARB_PHI_DELTA" \
      "--min-garb-pairs 3" \
      "--min-garb-samps 3" \
      "$INDIR/$jobname.truth.pickle" \
      "$INDIR/$jobname.params.json" \
      "$INDIR/$jobname.ssm" \
      "> $INDIR/$jobname.stdout" \
      "2>$INDIR/$jobname.stderr"
  done
  done
  done | parallel -j$PARA --halt 2 --eta
}

function detect_missedcna {
  for ssmfn in $INDIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    results=$RESULTSDIR/$runid
    mkdir -p $results && cd $results

    cmd=""
    cmd+="$PYTHON $BASEDIR/util/remove_bad_var_read_prob.py"
    cmd+=" --keep-existing-garbage"
    cmd+=" --verbose"
    cmd+=" --logbf-threshold 10"
    cmd+=" $INDIR/$runid.ssm"
    cmd+=" $INDIR/$runid.params.json"
    cmd+=" $results/$runid.params.json"
    cmd+=" > $results/$runid.stdout"
    cmd+=" 2>$results/$runid.stderr"
    echo $cmd
  done | parallel -j$PARA --halt 2 --eta
}

function eval_missedcna {
  for results in $RESULTSDIR/sim*; do
    runid_base=$(basename $results)
    for jsonfn in $results/*.params.json; do
      runid=$(basename $jsonfn | cut -d. -f1)
      cmd="$PYTHON $BASEDIR/comparison/pairtree/eval_garbage.py"
      cmd+=" $INDIR/$runid_base.params.json"
      cmd+=" $results/$runid.params.json"
      cmd+=" >$results/$runid.eval.json"
      echo $cmd
    done
  done | parallel -j$PARA --halt 2 --eta
}

function plot_missedcna {
  cd $RESULTSDIR
  echo "$PYTHON $BASEDIR/comparison/pairtree/plot_garbage.py" \
    "--plot-fn $RESULTSDIR/missedcna.html" \
    "--plot-bars" \
    "--filter-mindelta 0.01" \
    "sim*/*.eval.json" \
    "> $RESULTSDIR/missedcna.json" | parallel -j$PARA --halt 2 --eta
}

function main {
  make_inputs
  detect_missedcna
  eval_missedcna
  plot_missedcna
}

main
