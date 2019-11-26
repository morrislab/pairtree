#!/bin/bash
export OMP_NUM_THREADS=1

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
CITUP_DIR=$HOME/.apps/miniconda2/bin
JOB_DIR=$HOME/jobs

#CITUP_MODE=iter
CITUP_MODE=qip
USE_SUPERVARS=false
PARALLEL=1

BATCH=sims.smallalpha.citup
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.pairtree
#BATCH=steph.xeno.citup
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

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

function is_big_K {
  runid=$1
  if [[ $runid =~ K30 || $runid =~ K100 ]]; then
    return 0
  else
    return 1
  fi
}

function run_citup {
  for snvfn in $INDIR/*.snv; do
  #for runid in SJBALL022610 SJBALL031 SJETV047 SJETV010nohypermut SJBALL022614 SJMLL026 SJBALL022611 SJERG009 SJETV010 SJETV043 SJBALL022613 SJMLL039 SJBALL022612 SJBALL036 SJBALL022609 SJBALL022610steph SJETV010stephR1 SJETV010stephR1R2 SJETV010stephR2; do
    runid=$(basename $snvfn | cut -d. -f1)
    snvfn=$INDIR/$runid.snv
    outfn="$outdir/${runid}.results.hdf5"

    [[ -f $snvfn ]] || continue
    is_big_K $runid || continue
    [[ -f $outfn ]] && continue

    outdir="$OUTDIR/$runid"
    clusterfn="$INDIR/${runid}.cluster"
    let num_clusters=$(cat $clusterfn | sort | uniq | wc -l)+1

    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=$PARALLEL\n"
    cmd+="#SBATCH --mem=20GB\n"
    cmd+="#SBATCH --time=23:59:00\n"
    cmd+="#SBATCH --job-name citup_$runid\n"
    cmd+="#SBATCH --output=$JOB_DIR/slurm_citup_${runid}_%j.txt\n"
    cmd+="#SBATCH --partition=cpu\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="mkdir -p $outdir && cd $outdir && "
    cmd+="CITUPTMP=\$(mktemp -d --suffix .citup.$CITUP_MODE.$runid) "
    cmd+="TIMEFORMAT='%R %U %S'; time ("
    cmd+="PATH=$CITUP_DIR:$PATH run_citup_${CITUP_MODE}.py "
    cmd+="--submit local "
    cmd+="--maxjobs $PARALLEL "
    cmd+="--tmpdir \$CITUPTMP "
    cmd+="--min_nodes $num_clusters "
    cmd+="--max_nodes $num_clusters "
    cmd+="$INDIR/${runid}.snv "
    if [[ "$CITUP_MODE" == "qip" ]]; then
      cmd+="$clusterfn "
    fi
    cmd+="$outfn "
    cmd+=">$runid.stdout "
    cmd+="2>$runid.stderr "
    cmd+=") 2>$runid.time"
    cmd+="; cp -a \$CITUPTMP/log/latest/pipeline.log $outdir/"

    jobfn=$(mktemp)
    echo -e $cmd > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function convert_outputs {
  for outdir in $OUTDIR/*; do
    runid=$(basename $outdir)
    resultfn="$outdir/$runid.results.hdf5"
    [[ -f $resultfn ]] || continue

    cmd="cd $outdir && "
    cmd+="python3 $SCRIPTDIR/convert_outputs.py "
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
  done | parallel -j40 --halt 1 --eta
}

function main {
  #convert_inputs
  #run_citup
  convert_outputs
}

main
