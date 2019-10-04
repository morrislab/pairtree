#!/bin/bash
export OMP_NUM_THREADS=1

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
CITUP_DIR=$HOME/.apps/anaconda2/bin

#CITUP_MODE=iter
CITUP_MODE=qip
USE_SUPERVARS=false
PARALLEL=40

#BATCH=sims.citup
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
BATCH=steph.xeno.citup
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

if [[ "$USE_SUPERVARS" == "true" ]]; then
  suffix="supervars"
else
  suffix="rawvars"
fi

INDIR=$BASEDIR/scratch/inputs/$BATCH.$suffix
OUTDIR=$BASEDIR/scratch/results/$BATCH.$suffix.$CITUP_MODE

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    cmd="python3 $SCRIPTDIR/convert_inputs.py "
    if [[ "$USE_SUPERVARS" == "true" ]]; then
      cmd+="--use-supervars "
    fi
    cmd+="$PAIRTREE_INPUTS_DIR/$sampid.ssm "
    cmd+="$PAIRTREE_INPUTS_DIR/$sampid.params.json "
    cmd+="$INDIR/$sampid.snv "
    cmd+="$INDIR/$sampid.vid "
    cmd+="$INDIR/$sampid.cluster"
    echo $cmd
  done
}

function run_citup {
  #for snvfn in $INDIR/*.snv; do
  for runid in SJBALL022610 SJBALL031 SJETV047 SJETV010nohypermut SJBALL022614 SJMLL026 SJBALL022611 SJERG009 SJETV010 SJETV043 SJBALL022613 SJMLL039 SJBALL022612 SJBALL036 SJBALL022609 SJBALL022610steph SJETV010stephR1 SJETV010stephR1R2 SJETV010stephR2; do
    snvfn=$INDIR/$runid.snv
    [[ -f $snvfn ]] || continue
    #runid=$(basename $snvfn | cut -d. -f1)
    outdir="$OUTDIR/$runid"
    clusterfn="$INDIR/${runid}.cluster"
    let num_clusters=$(cat $clusterfn | sort | uniq | wc -l)+1
    [[ -f $outdir/${runid}.results.hdf5 ]] && continue

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
    cmd+="--weight-trees-by llh "
    if [[ "$USE_SUPERVARS" == "true" ]]; then
      cmd+="--use-supervars "
    fi
    if [[ "$CITUP_MODE" == "qip" ]]; then
      cmd+="--citup-clusters $INDIR/$runid.cluster "
    fi
    cmd+="--mutrel $outdir/$runid.mutrel.npz "
    cmd+="--mutphi $outdir/$runid.mutphi.npz "
    cmd+="$resultfn "
    cmd+="$INDIR/$runid.vid "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "

    echo $cmd
  done
}

function main {
  #convert_inputs
  run_citup
  #convert_outputs
}

main
