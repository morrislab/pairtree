#!/bin/bash
set -euo pipefail
module load gnu-parallel

BASEDIR=~/work/pairtree
JOBDIR=$SCRATCH/jobs
INDIR=$BASEDIR/scratch/inputs/sims.pwgs.supervars
OUTBASE=$BASEDIR/scratch/results/sims.pwgs.supervars
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PWGS_PATH=~/.apps/phylowgs
PARALLEL=40

function convert_inputs {
  mkdir -p $INDIR
  cp -a $BASEDIR/comparison/pwgs/empty.cnv $INDIR/

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $BASEDIR/comparison/pwgs/convert_inputs.py " \
      "--use-supervars" \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.ssm" \
      "$INDIR/$sampid.params.json"
  done | parallel -j$PARALLEL --halt 1
}

function run_pwgs {
  for ssmfn in $INDIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    jobname="pwgs_${runid}"
    OUTDIR=$OUTBASE/$runid

    treefn=$OUTDIR/chains/trees.zip
    [[ -f $treefn ]] && continue

    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $jobname"
      echo "#SBATCH --output=$JOBDIR/slurm_${jobname}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      # Must source ~/.bash_host_specific to get PATH set properly for
      # Miniconda.
      echo "source $HOME/.bash_host_specific && " \
        "module load gsl &&" \
        "mkdir -p $OUTDIR &&" \
        "cd $OUTDIR && " \
        "python2 $PWGS_PATH/multievolve.py" \
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

function convert_outputs {
  cd $OUTBASE
  for runid in *; do
    OUTDIR=$OUTBASE/$runid
    treefn=$OUTDIR/chains/trees.zip
    [[ -f $treefn ]] || continue
    [[ -f "$OUTDIR/$runid.mutrel.npz" ]] && continue

    echo "cd $OUTDIR &&" \
      "python2 $PWGS_PATH/write_results.py" \
      "--include-multiprimary" \
      "$runid" \
      "$treefn" \
      "$OUTDIR/$runid.{summ.json.gz,muts.json.gz,mutass.zip}" \
      ">  $OUTDIR/$runid.results.stdout" \
      "2> $OUTDIR/$runid.results.stderr &&" \
    "OMP_NUM_THREADS=1 python3 $BASEDIR/comparison/pwgs/convert_outputs.py" \
      "--weight-trees-by uniform" \
      "$PAIRTREE_INPUTS_DIR/$runid.{ssm,params.json}" \
      "$OUTDIR/$runid.{summ.json.gz,muts.json.gz,mutass.zip}" \
      "$OUTDIR/$runid.mutrel.npz" \
      ">>  $OUTDIR/$runid.results.stdout" \
      "2>> $OUTDIR/$runid.results.stderr"
  done | parallel -j$PARALLEL --halt 1
}

function main {
  #convert_inputs
  #run_pwgs
  convert_outputs
}

main
