#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel
PYTHON2=$HOME/.apps/miniconda2/bin/python2

BASEDIR=~/work/pairtree
JOBDIR=/tmp
PWGS_PATH=~/.apps/phylowgs
PARALLEL=40
NUM_CHAINS=1

BATCH=sims.smallalpha.pwgs.supervars
INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTBASE=$BASEDIR/scratch/results/$BATCH
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.pairtree

#BATCH=steph.pwgs.supervars
#BATCH=steph.pwgs.allvars
#INDIR=$BASEDIR/scratch/inputs/$BATCH
#OUTBASE=$BASEDIR/scratch/results/$BATCH
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

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
        "TIMEFORMAT='%R %U %S' && " \
        "time (OMP_NUM_THREADS=1 $PYTHON2 $PWGS_PATH/multievolve.py" \
        "--num-chains $NUM_CHAINS" \
        "--ssms $INDIR/${runid}.ssm" \
        "--cnvs $INDIR/empty.cnv" \
        "--params $INDIR/${runid}.params.json" \
        ">$runid.stdout" \
        "2>$runid.stderr) 2>$runid.time"
    ) #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done
}

function convert_outputs {
  USE_MULTICHAIN=$1

  cd $OUTBASE
  for runid in *; do
    OUTDIR=$OUTBASE/$runid

    if [ "$USE_MULTICHAIN" = true ]; then
      multichain_treefn=$OUTDIR/chains/trees.zip
      json_base=$runid.multichain
      treefn=$multichain_treefn
      results_base=$runid.pwgs_trees_multi
      [[ -f $multichain_treefn ]] || continue
    else
      chainpath=$(ls -d $OUTDIR/chains/chain* | sort --random-sort | head -n1)
      chain=$(basename $chainpath)
      json_base=$runid.single_$chain
      treefn=$chainpath/trees.zip
      results_base=$runid.pwgs_trees_single
      [[ $(unzip -l $treefn | grep tree_ | wc -l) > 0 ]] || continue
    fi

    [[ -f $OUTDIR/${results_base}.mutrel.npz && -f $OUTDIR/${results_base}.mutphi.npz ]] && continue

    cmd="cd $OUTDIR && export OMP_NUM_THREADS=1 "

    cmd+="&& $PYTHON2 $PWGS_PATH/write_results.py "
    cmd+="--include-multiprimary "
    cmd+="--keep-superclones "
    cmd+="$runid "
    cmd+="$treefn "
    cmd+="$OUTDIR/$json_base.{summ.json.gz,muts.json.gz,mutass.zip} "
    cmd+=">  $OUTDIR/$runid.results.stdout "
    cmd+="2> $OUTDIR/$runid.results.stderr "

    cmd+="&& OMP_NUM_THREADS=1 python3 $BASEDIR/comparison/pwgs/make_mutrel_and_mutphi.py "
    if [[ $OUTDIR =~ supervars ]]; then
      cmd+="--use-supervars "
    fi
    cmd+="--trees-mutrel $OUTDIR/${results_base}.mutrel.npz "
    cmd+="--phi $OUTDIR/${results_base}.mutphi.npz "
    cmd+="--pairtree-params $PAIRTREE_INPUTS_DIR/$runid.params.json "
    cmd+="--pairtree-ssms $PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$OUTDIR/$json_base.{summ.json.gz,muts.json.gz,mutass.zip} "
    cmd+=">>  $OUTDIR/$runid.results.stdout "
    cmd+="2>> $OUTDIR/$runid.results.stderr"

    echo $cmd
  done
}

function main {
  convert_inputs
  #run_pwgs

  #convert_outputs true
  # Don't run `parallel` with `--halt 1`, as some jobs will fail --
  # `write_results.py` will exit with non-zero when burn-in hasn't finished.
  #convert_outputs false | grep -v -e K30_ -e K100_ | parallel -j80 --eta
  #convert_outputs false #| grep    -e K30_ -e K100_ | parallel -j10 --eta
}

main
