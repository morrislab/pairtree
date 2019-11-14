#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
JOBDIR=~/jobs

PARALLEL=40
TREE_CHAINS=$PARALLEL
TREES_PER_CHAIN=3000
PHI_ITERATIONS=10000
PHI_FITTER=projection
THINNED_FRAC=1.0
BURNIN=0.333333

BATCH=sims.smallalpha.pairtree
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/$BATCH
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.lol1
#BATCH=steph.pairtree.multichain
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree
#PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/$BATCH

source $SCRIPTDIR/util.sh

function run_pairtree {
  logfn=$BASEDIR/scratch/logs/${BATCH}.$(date '+%s').log
  rm -f $logfn

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    outdir="$PAIRTREE_RESULTS_DIR/$runid"
    resultsfn="$outdir/$runid.results.npz"
    #is_big_K $runid && continue
    #is_run_complete $resultsfn && continue
    #[[ -f $resultsfn ]] || continue

    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=1:59:00"
      echo "#SBATCH --job-name $runid"
      echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

        #"--only-build-tensor" \
        #"--verbose" \
      echo "mkdir -p $outdir && cd $outdir && " \
        "TIMEFORMAT='%R %U %S'; time (" \
        "OMP_NUM_THREADS=1 " \
        "python3 $PROTDIR/bin/pairtree" \
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
    ) > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function convert_outputs {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    #is_big_K $runid || continue
    is_run_complete $resultsfn || continue

    jobfn=$(mktemp)
    (
        basepath="${outdir}/${runid}"

        cmd="cd $outdir && "
        cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutphis.py "
        cmd+="$resultsfn "
        cmd+="${PAIRTREE_INPUTS_DIR}/${runid}.ssm "
        cmd+="${basepath}.trees_mutphi.npz "
        cmd+=">${basepath}.mutphi_output_conversion.stdout "
        cmd+="2>${basepath}.mutphi_output_conversion.stderr"
        echo $cmd

        cmd="cd $outdir && "
        cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutrels.py "
        cmd+="--trees-mutrel ${basepath}.trees_mutrel.npz "
        cmd+="--clustrel-mutrel ${basepath}.clustrel_mutrel.npz "
        cmd+="$resultsfn "
        cmd+=">${basepath}.mutrel_output_conversion.stdout "
        cmd+="2>${basepath}.mutrel_output_conversion.stderr"
        echo $cmd
    ) #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done
}

function main {
  #run_pairtree #| grep python3 | parallel -j80 --halt 1 --eta
  convert_outputs | grep python3 | grep mutphi | sort --random-sort | parallel -j80 --halt 1 --eta
  convert_outputs | grep python3 | grep mutrel | sort --random-sort | parallel -j5 --halt 1 --eta
}

main
