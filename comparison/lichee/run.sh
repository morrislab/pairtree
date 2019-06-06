#!/bin/bash
command -v java > /dev/null || module load java

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
LICHEE_DIR=$HOME/.apps/lichee/LICHeE/release
NUM_TREES=3000

#BATCH=sims.lichee
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
BATCH=steph.xeno.lichee
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH

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
  mkdir -p $OUTDIR
  cd $OUTDIR

  for snvfn in $INDIR/*.snv; do
    runid=$(basename $snvfn | cut -d. -f1)

    echo "(command -v java > /dev/null || module load java) && " \
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
    for weight_trees_by in llh; do
      cmd="python3 $SCRIPTDIR/convert_outputs.py "
      cmd+="--weight-trees-by $weight_trees_by "
      cmd+="--mutrel $OUTDIR/$runid.$weight_trees_by.mutrel.npz "
      cmd+="$INDIR/$runid.snv "
      cmd+="$OUTDIR/$runid.trees "
      cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "

      # Only output mutrels for `K <= 10`.
      if ! [[ $runid =~ K30 || $runid =~ K100 ]]; then
        echo $cmd
      fi
    done

    cmd="python3 $SCRIPTDIR/convert_outputs.py "
    cmd+="--structures $OUTDIR/$runid.params.json "
    cmd+="$INDIR/$runid.snv "
    cmd+="$OUTDIR/$runid.trees "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "
    echo $cmd
  done
}

function compute_phis {
  cd $OUTDIR
  for structfn in $OUTDIR/*.params.json; do
    runid=$(basename $structfn | cut -d. -f1)
    resultsfn=$OUTDIR/$runid.results.npz

    cmd="LD_LIBRARY_PATH=$HOME/tmp/jose/bin:$LD_LIBRARY_PATH "
    cmd+="python3 $BASEDIR/basic.py "
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
  cd $OUTDIR

  for resultsfn in $OUTDIR/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    ssmfn=${PAIRTREE_INPUTS_DIR}/${runid}.ssm 

    for weight_trees_by in llh; do
      mutphifn=$OUTDIR/$runid.$weight_trees_by.mutphi.npz

      cmd="python3 $BASEDIR/comparison/pairtree/make_mutphis.py "
      cmd+="--weight-trees-by $weight_trees_by "
      cmd+="--ssms $ssmfn "
      cmd+="$resultsfn "
      cmd+="$mutphifn "

      cmd+="&& python3 $BASEDIR/comparison/impute_missing_mutphis.py "
      cmd+="$ssmfn "
      cmd+="$PAIRTREE_INPUTS_DIR/${runid}.params.json "
      cmd+="$mutphifn "

      echo $cmd
    done
  done
}

function main {
  #convert_inputs
  #run_lichee
  #convert_outputs
  #compute_phis
  compute_mutphis
}

main
