#!/bin/bash
set -euo pipefail
module load gnu-parallel

BASEDIR=~/work/pairtree
JOBDIR=$SCRATCH/jobs
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
#BATCH=sims.pastri.informative
#BATCH=sims.pastri.uniform
BATCH=$1
INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH
PARALLEL=80
NUM_ITERS=10000

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $BASEDIR/comparison/pastri/convert_inputs.py " \
      "--uniform-proposal" \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.counts" \
      "$INDIR/$sampid.proposal"
  done | parallel -j$PARALLEL --halt 1
}

function run {
  mkdir -p $OUTDIR

  for countsfn in $INDIR/*.counts; do
    runid=$(basename $countsfn | cut -d. -f1)
    jobname="steph_pastri_${runid}"

    (
      # Must source ~/.bash_host_specific to get PATH set properly for
      # Miniconda.
      echo "source $HOME/.bash_host_specific && " \
        "cd $OUTDIR && " \
        "python2 $HOME/.apps/pastri/src/RunPASTRI.py" \
        "--output_prefix $runid" \
        "--num_iters $NUM_ITERS" \
        "$INDIR/${runid}.counts" \
        "$INDIR/${runid}.proposal" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) 
  done | parallel -j$PARALLEL --joblog $SCRATCH/tmp/$BATCH.log
}

function main {
  #convert_inputs
  run
}

main