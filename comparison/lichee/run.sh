#!/bin/bash
command -v java > /dev/null || module load java

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
LICHEE_DIR=$HOME/.apps/lichee/LICHeE/release
NUM_TREES=3000

BATCH=sims.lichee
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/convert_inputs.py " \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.snv" \
      "$INDIR/$sampid.cluster"
  done
}

function run_lichee {
  mkdir -p $OUTDIR
  cd $OUTDIR

  for snvfn in $INDIR/*.snv; do
    runid=$(basename $snvfn | cut -d. -f1)
    echo "module load java && " \
      "cd $OUTDIR &&" \
      "TIMEFORMAT='%R %U %S'; time (java -jar $LICHEE_DIR/lichee.jar" \
      "-build" \
      "-i $INDIR/$runid.snv" \
      "-o $OUTDIR/$runid.trees" \
      "-clustersFile $INDIR/$runid.cluster" \
      "-maxVAFAbsent 0.005" \
      "-minVAFPresent 0.005" \
      "-n 0" \
      "-s $NUM_TREES" \
      ">  $OUTDIR/$runid.stdout" \
      "2> $OUTDIR/$runid.stderr) 2>$runid.time"
  done
}

function convert_outputs {
  cd $OUTDIR
  for treesfn in $OUTDIR/*.trees; do
    runid=$(basename $treesfn | cut -d. -f1)
    cmd="python3 $SCRIPTDIR/convert_outputs.py "
    cmd+="--weight-trees-by llh "
    if ! [[ $runid =~ K30 || $runid =~ K100 ]]; then
      cmd+="--mutrel $OUTDIR/$runid.mutrel.npz "
    fi
    cmd+="--structures $OUTDIR/$runid.params.json "
    cmd+="$OUTDIR/$runid.trees "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "
    echo $cmd
  done
}

function main {
  #convert_inputs
  #run_lichee
  convert_outputs
}

main
