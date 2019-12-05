#!/bin/bash
set -euo pipefail
module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.patient.withgarb.pairtree
HANDBUILTDIR=~/work/steph/data/handbuilt_trees
PHYLOSTEPH_INPUTS_DIR=~/work/steph/data/inputs/steph.xenos.nocns
PARALLEL=40

function convert_phylosteph_inputs {
  mkdir -p $PAIRTREE_INPUTS_DIR
  cd $PAIRTREE_INPUTS_DIR

  for ssmfn in $PHYLOSTEPH_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
      #"--discard-garbage" \
    echo "python3 $SCRIPTDIR/convert_phylosteph_inputs_to_pairtree.py " \
      "--handbuilt $HANDBUILTDIR/$runid.json" \
      "$PHYLOSTEPH_INPUTS_DIR/$runid.{sampled.ssm,params.json}" \
      "$PAIRTREE_INPUTS_DIR/$runid.{ssm,params.json}"
  done #| parallel -j$PARALLEL --halt 1
}

function main {
  convert_phylosteph_inputs
}

main
