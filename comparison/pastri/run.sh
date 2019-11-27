#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
PARALLEL=40
NUM_ITERS=10000
PASTRI_DIR=$HOME/.apps/pastri
JOB_DIR=$HOME/jobs
PYTHON2=$HOME/.apps/miniconda2/bin/python2
PYTHON3=$HOME/.apps/miniconda3/bin/python3

BATCH=sims.smallalpha.pastri
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.pairtree
#BATCH=steph.xeno.pastri.hbclusters
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
      #"--uniform-proposal" \
    echo "$PYTHON3 $SCRIPTDIR/convert_inputs.py " \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.counts" \
      "$INDIR/$sampid.proposal"
  done #| parallel -j$PARALLEL --halt 1
}

function run_pastri {
  for countsfn in $INDIR/*.counts; do
    runid=$(basename $countsfn | cut -d. -f1)

    jobname="pastri_${runid}"
    outpath=$OUTDIR/$runid
    mkdir -p $outpath

    [[ -f $OUTDIR/$runid.trees ]] && continue

    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=1\n"
    cmd+="#SBATCH --mem=20GB\n"
    cmd+="#SBATCH --time=23:59:00\n"
    cmd+="#SBATCH --job-name pastri_$runid\n"
    cmd+="#SBATCH --output=$JOB_DIR/slurm_pastri_${runid}_%j.txt\n"
    cmd+="#SBATCH --partition=cpu\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="cd $outpath && "
    cmd+="TIMEFORMAT='%R %U %S'; time (OMP_NUM_THREADS=1 "
    cmd+="$PYTHON2 $PASTRI_DIR/src/RunPASTRI.py "
    cmd+="--output_prefix $runid "
    cmd+="--num_iters $NUM_ITERS "
    cmd+="$INDIR/${runid}.counts "
    cmd+="$INDIR/${runid}.proposal "
    cmd+=">$runid.stdout "
    cmd+="2>$runid.stderr) 2>$runid.time"

    jobfn=$(mktemp)
    echo -e $cmd > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function get_F_and_C {
  for treesfn in $OUTDIR/*/*.trees; do
    runid=$(basename $treesfn | cut -d. -f1)
    outpath=$(dirname $treesfn)
    # Trees in $treesfn are ranked by likelihood. `get_F_and_C.py` takes tree
    # rank as a parameter. Thus, if we count N trees with LH > 0, we know their
    # ranks are [0, 1, ..., N-1].
    valid_count=$(cat $treesfn | grep '^>' | grep -v -- '-inf$' | wc -l)
    for idx in $(seq $valid_count); do
      echo "cd $outpath && " \
        "$PYTHON2 $PASTRI_DIR/src/get_F_and_C.py" \
        "-i $idx" \
        "-o $outpath/$runid" \
        "$INDIR/${runid}.counts" \
        "$treesfn" \
        "$outpath/${runid}.fsamples" \
        "> $outpath/${runid}.${idx}.output_conversion.stdout" \
        "2>$outpath/${runid}.${idx}.output_conversion.stderr"
    done
  done | parallel -j$PARALLEL
}

function convert_outputs {
  for treesfn in $OUTDIR/*/*.trees; do
    runid=$(basename $treesfn | cut -d. -f1)
    outpath=$(dirname $treesfn)
    echo "cd $outpath && " \
      "OMP_NUM_THREADS=1 $PYTHON3 $SCRIPTDIR/convert_outputs.py" \
      "$runid" \
      "$PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$treesfn" \
      "${outpath}/${runid}.neutree.pickle"
  done | parallel -j$PARALLEL --halt 1 --eta
}

function main {
  #convert_inputs
  #run_pastri
  #get_F_and_C
  convert_outputs
}

main
