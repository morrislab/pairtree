#!/bin/bash
BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
CITUP_DIR=$HOME/.apps/anaconda2/bin

#CITUP_MODE=qip
CITUP_MODE=iter
PARALLEL=80

BATCH=sims.citup
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
#BATCH=steph.xeno.citup
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH.$CITUP_MODE

export OMP_NUM_THREADS=1

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/convert_inputs.py " \
      "--use-supervars" \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.snv" \
      "$INDIR/$sampid.vid" \
      "$INDIR/$sampid.cluster"
  done
}

function run_citup {
  for snvfn in $INDIR/*.snv; do
    runid=$(basename $snvfn | cut -d. -f1)
    outdir="$OUTDIR/$runid"
    clusterfn+="$INDIR/${runid}.cluster "
    # Does this need to be `clusters + 1`?
    let num_clusters=$(cat $clusterfn | sort | uniq | wc -l)+1

    cmd="mkdir -p $outdir && cd $outdir &&"
    cmd+="TIMEFORMAT='%R %U %S'; time ("
    cmd+="PATH=$HOME/.apps/anaconda2/bin:$PATH run_citup_${CITUP_MODE}.py " 
    cmd+="--submit local "
    cmd+="--maxjobs $PARALLEL "
    cmd+="--tmpdir \$(mktemp -d --suffix .citup.$CITUP_MODE.$runid) "
    cmd+="--min_nodes $num_clusters "
    cmd+="--max_nodes $num_clusters "
    cmd+="$INDIR/${runid}.snv "
    if [[ "$CITUP_MODE" == "qip" ]]; then
      cmd+="$clusterfn "
    fi
    cmd+="$outdir/${runid}.results.hdf5 "
    cmd+=">$runid.stdout "
    cmd+="2>$runid.stderr "
    cmd+=") 2>$runid.time"

    echo $cmd
  done
}

function convert_outputs {
  for outdir in $OUTDIR/*; do
    runid=$(basename $outdir)
    resultfn="$outdir/$runid.results.hdf5"
    [[ -f $resultfn ]] || continue

    cmd="cd $outdir && "
    cmd+="python3 $SCRIPTDIR/convert_outputs.py "
    cmd+="$resultfn "
    cmd+="$INDIR/$runid.vid "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "

    echo $cmd
  done
}

function main {
  #convert_inputs
  #run_citup
  convert_outputs
}

main
