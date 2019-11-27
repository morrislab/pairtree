#!/bin/bash
#command -v parallel > /dev/null || module load gnu-parallel
PYTHON2=$HOME/.apps/miniconda2/bin/python2

BASEDIR=~/work/pairtree
JOB_DIR=$HOME/jobs
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

    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=1\n"
    cmd+="#SBATCH --mem=20GB\n"
    cmd+="#SBATCH --time=23:59:00\n"
    cmd+="#SBATCH --job-name pwgs_$runid\n"
    cmd+="#SBATCH --output=$JOB_DIR/slurm_pwgs_${runid}_%j.txt\n"
    cmd+="#SBATCH --partition=cpu\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="mkdir -p $OUTDIR && "
    cmd+="cd $OUTDIR && "
    cmd+="TIMEFORMAT='%R %U %S' && "
    cmd+="time (LD_LIBRARY_PATH=$HOME/.apps/gsl/lib:$LD_LIBRARY_PATH OMP_NUM_THREADS=1 $PYTHON2 $PWGS_PATH/multievolve.py "
    cmd+="--num-chains $NUM_CHAINS "
    cmd+="--ssms $INDIR/${runid}.ssm "
    cmd+="--cnvs $INDIR/empty.cnv "
    cmd+="--params $INDIR/${runid}.params.json "
    cmd+=">$runid.stdout "
    cmd+="2>$runid.stderr) 2>$runid.time"

    jobfn=$(mktemp)
    echo -e $cmd > $jobfn
    sbatch $jobfn
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

    cmd="cd $OUTDIR && export OMP_NUM_THREADS=1 "

    cmd+="&& $PYTHON2 $PWGS_PATH/write_results.py "
    cmd+="--include-multiprimary "
    cmd+="--keep-superclones "
    cmd+="$runid "
    cmd+="$treefn "
    cmd+="$OUTDIR/$json_base.{summ.json.gz,muts.json.gz,mutass.zip} "

    cmd+="&& OMP_NUM_THREADS=1 python3 $BASEDIR/comparison/pwgs/convert_outputs.py "
    if [[ $OUTDIR =~ supervars ]]; then
      cmd+="--use-supervars "
    fi
    cmd+="$OUTDIR/$json_base.{summ.json.gz,muts.json.gz,mutass.zip} "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "
    cmd+="$OUTDIR/$runid.neutree.pickle "

    echo $cmd
  done
}

function main {
  #convert_inputs
  #run_pwgs

  #convert_outputs true
  # Don't run `parallel` with `--halt 1`, as some jobs will fail --
  # `write_results.py` will exit with non-zero when burn-in hasn't finished.
  convert_outputs false | parallel -j80 --eta
}

main
