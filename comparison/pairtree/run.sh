#!/bin/bash
set -euo pipefail
#command -v parallel > /dev/null || module load gnu-parallel

BASEDIR=~/work/pairtree
EXPTSDIR=~/work/pairtree-experiments
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
JOBDIR=~/jobs
PYTHON3=$HOME/.apps/miniconda3/bin/python3

PARALLEL=40
TREE_CHAINS=$PARALLEL
TREES_PER_CHAIN=3000
PHI_ITERATIONS=10000
PHI_FITTER=projection
THINNED_FRAC=1.0
BURNIN=0.333333

BATCH=sims.smallalpha.pairtree
PAIRTREE_INPUTS_DIR=$EXPTSDIR/inputs/$BATCH
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.${PHI_FITTER}

#BATCH=steph.congraph
#PAIRTREE_INPUTS_DIR=$EXPTSDIR/inputs/steph.xeno.pairtree
#PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.${PHI_FITTER}

source $SCRIPTDIR/util.sh

function run_pairtree {
  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    outdir="$PAIRTREE_RESULTS_DIR/$runid"
    resultsfn="$outdir/$runid.results.npz"
    #[[  $runid =~ K100_ ]] || continue
    #is_big_K $runid && continue
    is_run_complete $resultsfn && continue
    #[[ -f $resultsfn ]] || continue

    jobfn=$(mktemp)
    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=$PARALLEL\n"
    cmd+="#SBATCH --mem=20GB\n"
    cmd+="#SBATCH --time=95:59:00\n"
    cmd+="#SBATCH --job-name pairtree_$runid\n"
    cmd+="#SBATCH --output=$JOBDIR/slurm_pairtree_${runid}_%j.txt\n"
    cmd+="#SBATCH --partition=cpu\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="mkdir -p $outdir && cd $outdir && "
    cmd+="TIMEFORMAT='%R %U %S'; time ( "
    cmd+="OMP_NUM_THREADS=1  "
    cmd+="$PYTHON3 $BASEDIR/bin/pairtree "
    cmd+="--seed 1 "
    cmd+="--parallel $PARALLEL "
    cmd+="--tree-chains $TREE_CHAINS "
    cmd+="--trees-per-chain $TREES_PER_CHAIN "
    cmd+="--phi-iterations $PHI_ITERATIONS "
    cmd+="--phi-fitter $PHI_FITTER "
    cmd+="--params $PAIRTREE_INPUTS_DIR/${runid}.params.json "
    cmd+="--thinned-frac $THINNED_FRAC "
    cmd+="--burnin $BURNIN "
    cmd+="--gamma 0.7 "
    cmd+="--zeta 0.7 "
    cmd+="--iota 0.7 "
    #cmd+="--only-build-tensor "
    #cmd+="--verbose "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$resultsfn "
    cmd+=">$runid.stdout "
    cmd+="2>$runid.stderr) 2>$runid.time "

    echo -e $cmd > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function create_neutree {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    basepath="${outdir}/${runid}"
    #is_big_K $runid || continue
    is_run_complete $resultsfn || continue

    (
      cmd="cd $outdir && "
      cmd+="OMP_NUM_THREADS=1 $PYTHON3 $SCRIPTDIR/convert_outputs.py "
      cmd+="$resultsfn "
      cmd+="${PAIRTREE_INPUTS_DIR}/${runid}.hbstruct.params.json "
      cmd+="${basepath}.neutree.pickle "
      echo $cmd
    )
  done
}

function create_mutrel_from_clustrel {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    basepath="${outdir}/${runid}"
    #is_big_K $runid || continue
    is_run_complete $resultsfn || continue

    (
      cmd="cd $outdir && "
      cmd+="OMP_NUM_THREADS=1 $PYTHON3 $SCRIPTDIR/make_mutrel_from_clustrel.py "
      cmd+="$resultsfn "
      cmd+="${basepath}.clustrel_mutrel.npz "
      echo $cmd
    )
  done
}

function main {
  run_pairtree #| grep python3 | parallel -j2 --halt 1 --eta
  #create_neutree | parallel -j80 --halt 1 --eta
  #create_mutrel_from_clustrel | sort --random-sort | parallel -j5 --halt 1 --eta
}

main
