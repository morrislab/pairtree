#!/bin/bash
set -euo pipefail
#module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
JOBDIR=~/jobs
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/sims.pairtree.fixedclusters

PARALLEL=40
TREE_CHAINS=40
BURNIN_PER_CHAIN=1000
TREES_PER_CHAIN=2000
PHI_ITERATIONS=10000

function is_run_complete {
  if ([[ -f $1 ]] && unzip -l $1 | grep "adjm.npy" > /dev/null); then
    return 0
  else
    return 1
  fi
}

function run_pairtree {
  mkdir -p $PAIRTREE_RESULTS_DIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    resultsfn="$PAIRTREE_RESULTS_DIR/$runid.results.npz"
    is_run_complete $resultsfn && continue

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
        "--tree-chains $TREE_CHAINS" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--burnin-per-chain $BURNIN_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--params $PAIRTREE_INPUTS_DIR/${runid}.params.json" \
        "$PAIRTREE_INPUTS_DIR/$runid.ssm" \
        "$resultsfn" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function convert_outputs {
  for resultsfn in $PAIRTREE_RESULTS_DIR/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    resultsfn="$PAIRTREE_RESULTS_DIR/$runid.results.npz"
    is_run_complete $resultsfn || continue
    #jobfn=$(mktemp)
    (
        #echo "#!/bin/bash"
        #echo "#SBATCH --nodes=1"
        #echo "#SBATCH --ntasks=$PARALLEL"
        #echo "#SBATCH --time=23:59:00"
        #echo "#SBATCH --job-name convert_$runid"
        #echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
        #echo "#SBATCH --mail-type=NONE"

        echo "cd $PAIRTREE_RESULTS_DIR &&" \
          "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/convert_outputs.py" \
          "--clustrel-mutrel $PAIRTREE_RESULTS_DIR/$runid.pairtree_clustrel.mutrel.npz" \
          "$resultsfn" \
          ">$PAIRTREE_RESULTS_DIR/$runid.output_conversion.stdout" \
          "2>$PAIRTREE_RESULTS_DIR/$runid.output_conversion.stderr"
        for tree_weights in uniform llh; do
          echo "cd $PAIRTREE_RESULTS_DIR &&" \
            "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/convert_outputs.py" \
            "--weight-trees-by $tree_weights" \
            "--trees-mutrel $PAIRTREE_RESULTS_DIR/$runid.pairtree_trees_$tree_weights.mutrel.npz" \
            "--phi $PAIRTREE_RESULTS_DIR/$runid.pairtree_trees_$tree_weights.mutphi.npz" \
            "$resultsfn" \
            ">>$PAIRTREE_RESULTS_DIR/$runid.output_conversion.stdout" \
            "2>>$PAIRTREE_RESULTS_DIR/$runid.output_conversion.stderr"
        done
    ) #> $jobfn
    #sbatch $jobfn
    #rm $jobfn
  done
}

function main {
  run_pairtree
  #convert_outputs
}

main
