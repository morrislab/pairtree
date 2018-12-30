#!/bin/bash
set -euo pipefail
module load gnu-parallel

BASEDIR=~/work/pairtree
JOBDIR=$SCRATCH/jobs
INDIR=$BASEDIR/scratch/inputs/steph.xeno.unclust.pwgs
OUTBASE=$BASEDIR/scratch/results/steph.xeno.unclust.pwgs
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree
PARALLEL=40

function convert_inputs {
  mkdir -p $INDIR
  cp -a $BASEDIR/comparison/pwgs/empty.cnv $INDIR/

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $BASEDIR/comparison/pwgs/convert_inputs.py " \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.ssm" \
      "$INDIR/$sampid.params.json"
  done | parallel -j$PARALLEL --halt 1
}

function run_steph {
  for ssmfn in $INDIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    jobname="steph_pwgs_${runid}"
    OUTDIR=$OUTBASE/$runid

    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=24:00:00"
      echo "#SBATCH --job-name $jobname"
      echo "#SBATCH --output=$JOBDIR/slurm_${jobname}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      # Must source ~/.bash_host_specific to get PATH set properly for
      # Miniconda.
      echo "source $HOME/.bash_host_specific && " \
        "module load gsl &&" \
        "mkdir -p $OUTDIR &&" \
        "cd $OUTDIR && " \
        "python2 $HOME/.apps/phylowgs/multievolve.py" \
        "--num-chains $PARALLEL" \
        "--ssms $INDIR/${runid}.ssm" \
        "--cnvs $INDIR/empty.cnv" \
        "--params $INDIR/${runid}.params.json" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function create_mutrels {
  cd $OUTBASE
  for runid in *; do
    OUTDIR=$OUTBASE/$runid
    [[ -f $OUTDIR/$runid.summ.json.gz ]] || continue
    echo "cd $OUTDIR &&" \
      "python3 $BASEDIR/comparison/pwgs/convert_outputs.py" \
      "$PAIRTREE_INPUTS_DIR/$runid.{ssm,params.json}" \
      "$OUTDIR/$runid.{summ.json.gz,muts.json.gz,mutass.zip}" \
      "$OUTDIR/$runid.mutrel.npz"
  done | parallel -j$PARALLEL --halt 1
}

function main {
  #convert_inputs
  run_steph
  #create_mutrels
}

main
