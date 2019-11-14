#!/bin/bash
command -v java > /dev/null || module load java

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
LICHEE_DIR=$HOME/.apps/lichee/LICHeE/release
NUM_TREES=3000

BATCH=sims.smallalpha.lichee
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.pairtree
#BATCH=steph.xeno.lichee
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTBASE=$BASEDIR/scratch/results/$BATCH

export OMP_NUM_THREADS=1

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
  mkdir -p $OUTBASE
  cd $OUTBASE

  for snvfn in $INDIR/*.snv; do
    runid=$(basename $snvfn | cut -d. -f1)
    outd="$OUTBASE/$runid"
    echo \
      "mkdir -p $outd && cd $outd &&" \
      "TIMEFORMAT='%R %U %S'; time (java -jar $LICHEE_DIR/lichee.jar" \
      "-build" \
      "-i $INDIR/$runid.snv" \
      "-o $outd/$runid.trees" \
      "-clustersFile $INDIR/$runid.cluster" \
      "-maxVAFAbsent 0.005" \
      "-minVAFPresent 0.005" \
      "-n 0" \
      "-s $NUM_TREES" \
      ">  $outd/$runid.stdout" \
      "2> $outd/$runid.stderr) 2>$outd/$runid.time"
  done
}

function convert_outputs {
  cd $OUTBASE
  for treesfn in $OUTBASE/*/*.trees; do
    runid=$(basename $treesfn | cut -d. -f1)
    outd=$(dirname $treesfn)

    cmd="python3 $SCRIPTDIR/convert_outputs.py "
    cmd+="--mutrel $outd/$runid.mutrel.npz "
    cmd+="--structures $outd/$runid.params.json "
    cmd+="$INDIR/$runid.snv "
    cmd+="$outd/$runid.trees "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "

    echo $cmd
  done
}

function compute_phis {
  cd $OUTBASE
  for structfn in $OUTBASE/*/*.params.json; do
    runid=$(basename $structfn | cut -d. -f1)
    outd=$(dirname $structfn)
    resultsfn=$outd/$runid.results.npz

    cmd="python3 $BASEDIR/bin/pairtree "
    cmd+="--seed 1 "
    cmd+="--phi-fitter projection "
    cmd+="--parallel 0 "
    cmd+="--params $structfn "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$resultsfn "

    cmd+="&& python3 $SCRIPTDIR/replace_llh.py "
    cmd+="$structfn "
    cmd+="$resultsfn "

    echo $cmd
  done
}

function compute_mutphis {
  cd $OUTBASE

  for resultsfn in $OUTBASE/*/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    outd=$(dirname $resultsfn)
    ssmfn=${PAIRTREE_INPUTS_DIR}/${runid}.ssm 
    mutphifn=$outd/$runid.mutphi.npz

    cmd="python3 $BASEDIR/comparison/pairtree/make_mutphis.py "
    cmd+="$resultsfn "
    cmd+="$ssmfn "
    cmd+="$mutphifn "

    cmd+="&& python3 $BASEDIR/comparison/impute_missing_mutphis.py "
    cmd+="$ssmfn "
    cmd+="$PAIRTREE_INPUTS_DIR/${runid}.params.json "
    cmd+="$mutphifn "

    echo $cmd
  done
}

function main {
  convert_inputs
  #run_lichee
  #convert_outputs | sort --random-sort | parallel -j8 --halt 1 --eta
  #compute_phis | parallel -j40 --halt 1 --eta
  #compute_mutphis | parallel -j80 --halt 1 --eta
}

main
