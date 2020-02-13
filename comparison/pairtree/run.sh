#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
JOBDIR=~/jobs

PARALLEL=40
TREE_CHAINS=$PARALLEL
TREES_PER_CHAIN=3000
PHI_ITERATIONS=10000
PHI_FITTER=rprop
THINNED_FRAC=1.0
BURNIN=0.333333

BATCH=sims.smallalpha.pairtree
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/$BATCH
TRUTH_DIR=$BASEDIR/scratch/results/sims.smallalpha.truth
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.rprop

#TREE_TYPE=xeno
#BATCH=steph.${TREE_TYPE}.pairtree.multichain.testlol
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.${TREE_TYPE}.pairtree.nostructs
#PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/$BATCH

source $SCRIPTDIR/util.sh

function run_pairtree {
  logfn=$BASEDIR/scratch/logs/${BATCH}.$(date '+%s').log
  rm -f $logfn

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    outdir="$PAIRTREE_RESULTS_DIR/$runid"
    resultsfn="$outdir/$runid.results.npz"
    #[[  $runid =~ K100_ ]] && continue
    #is_big_K $runid && continue
    is_run_complete $resultsfn && continue
    #[[ -f $resultsfn ]] || continue

    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $runid"
      echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

        #"--only-build-tensor" \
        #"--verbose" \
      echo "mkdir -p $outdir && cd $outdir && " \
        "TIMEFORMAT='%R %U %S'; time (" \
        "OMP_NUM_THREADS=1 " \
        "python3 $BASEDIR/bin/pairtree" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $TREE_CHAINS" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--phi-fitter $PHI_FITTER" \
        "--params $PAIRTREE_INPUTS_DIR/${runid}.params.json" \
        "--thinned-frac $THINNED_FRAC" \
        "--burnin $BURNIN" \
        "--gamma 0.7" \
        "--zeta 0.7" \
        "--iota 0.7" \
        "$PAIRTREE_INPUTS_DIR/$runid.ssm" \
        "$resultsfn" \
        ">$runid.stdout" \
        "2>$runid.stderr) 2>$runid.time"
    ) #> $jobfn
    #sbatch $jobfn
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
      cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/convert_outputs.py "
      cmd+="$resultsfn "
      cmd+="${PAIRTREE_INPUTS_DIR}/${runid}.params.json "
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
      cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutrel_from_clustrel.py "
      cmd+="$resultsfn "
      cmd+="${basepath}.clustrel_mutrel.npz "
      echo $cmd
    )
  done
}

function main {
  run_pairtree #| grep python3 | parallel -j80 --halt 1 --eta
  #create_neutree | parallel -j80 --halt 1 --eta
  #create_mutrel_from_clustrel | sort --random-sort | parallel -j5 --halt 1 --eta
}

main
