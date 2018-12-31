#!/bin/bash
set -euo pipefail
module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/sims.pairtree.fixedclusters
PARALLEL=40

function run_pairtree {
  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $runid"
      echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      echo "cd $PAIRTREE_RESULTS_DIR && " \
        "OMP_NUM_THREADS=1 python3 $PROTDIR/basic.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $PARALLEL" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--params $PAIRTREE_INPUTS_DIR/${runid}.params.json" \
        "$PAIRTREE_INPUTS_DIR/$runid.ssm" \
        "$PAIRTREE_RESULTS_DIR/$runid.results.npz" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done
}


function convert_outputs {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    echo "cd $PAIRTREE_RESULTS_DIR &&" \
      "python3 $SCRIPTDIR/convert_outputs.py" \
      "$PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$PAIRTREE_RESULTS_DIR/$runid.results.npz" \
      "$PAIRTREE_RESULTS_DIR/$runid.{pairs,trees}.mutrel.npz"
  done | parallel -j$PARALLEL --halt 1
}

function main {
  run_pairtree
  #convert_outputs
}

main