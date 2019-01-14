#!/bin/bash
set -euo pipefail
module load gnu-parallel

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
JOBDIR=$SCRATCH/jobs
INBASE=$BASEDIR/scratch/inputs/sims.sciclone
OUTDIR=$BASEDIR/scratch/results/sims.sciclone
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PARALLEL=40

function convert_inputs {
  mkdir -p $INBASE

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    indir=$INBASE/$sampid

    echo "mkdir -p $indir &&" \
      "python3 $SCRIPTDIR/convert_inputs.py" \
      "$ssmfn" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$indir"
  done | parallel -j$PARALLEL --halt 1
}

function run_sciclone {
  mkdir -p $OUTDIR
  cd $INBASE

  for sampid in *; do
    echo "module load r &&" \
      "Rscript --vanilla $SCRIPTDIR/run_sciclone.R" \
      "$INBASE/$sampid/*.dat" \
      "$OUTDIR/$sampid.results.txt" \
      "$OUTDIR/$sampid.sampnames.json" \
      ">  $OUTDIR/$sampid.stdout" \
      "2> $OUTDIR/$sampid.stderr"
  done | parallel -j$PARALLEL
}

function convert_outputs {
  for resultsfn in $OUTDIR/*.results.txt; do
    sampid=$(basename $resultsfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/convert_outputs.py" \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$resultsfn" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$OUTDIR/$sampid.params.json"
  done | parallel -j$PARALLEL --halt 1
}

function main {
  #convert_inputs
  #run_sciclone
  convert_outputs
}

main
