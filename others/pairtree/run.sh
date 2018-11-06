#!/bin/bash
set -euo pipefail
module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
JOBDIR=$SCRATCH/jobs
HANDBUILTDIR=~/work/steph/data/handbuilt_trees
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.nogarb.pairtree
PHYLOSTEPH_INPUTS_DIR=~/work/steph/data/inputs/steph.xenos.nocns
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/steph.xeno.nogarb.pairtree
TRUTH_DIR=$BASEDIR/scratch/results/steph.xeno.nogarb.truth
PARALLEL=40

function convert_phylosteph_inputs {
  mkdir -p $PAIRTREE_INPUTS_DIR
  cd $PAIRTREE_INPUTS_DIR

  for ssmfn in $PHYLOSTEPH_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/convert_phylosteph_inputs_to_pairtree.py " \
      "--handbuilt $HANDBUILTDIR/$runid.json" \
      "$PHYLOSTEPH_INPUTS_DIR/$runid.{sampled.ssm,params.json}" \
      "$PAIRTREE_INPUTS_DIR/$runid.{ssm,params.json}"
  done | parallel -j$PARALLEL --halt 1
}

function make_truth_mutrels {
  mkdir -p $TRUTH_DIR
  for paramsfn in $PAIRTREE_INPUTS_DIR/*.params.json; do
    runid=$(basename $paramsfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/make_truth_mutrel.py" \
      "$PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$TRUTH_DIR/$runid.mutrel.npz"
  done | parallel -j$PARALLEL --halt 1
}

function create_mutrels {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    echo "cd $PAIRTREE_RESULTS_DIR &&" \
      "python3 $SCRIPTDIR/convert_outputs.py" \
      "$PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$PAIRTREE_RESULTS_DIR/$runid.results.npz" \
      "$PAIRTREE_RESULTS_DIR/$runid.{pairs,trees}.mutrel.npz"
  done | parallel -j$PARALLEL --halt 1
}

function main {
  #convert_phylosteph_inputs
  #make_truth_mutrels
  create_mutrels
}

main
