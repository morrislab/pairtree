#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
JOBDIR=~/jobs

PARALLEL=1
TREE_CHAINS=$PARALLEL
BURNIN_PER_CHAIN=1000
TREES_PER_CHAIN=2000
PHI_ITERATIONS=10000
PHI_FITTER=$1

#BATCH=sims.pairtree.${PHI_FITTER}.singlechain
BATCH=sims.pairtree
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/$BATCH
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree
#PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/steph.xeno.withgarb.pairtree

function is_run_complete {
  if ([[ -f $1 ]] && unzip -l $1 | grep "adjm.npy" > /dev/null); then
    return 0
  else
    return 1
  fi
}

function run_pairtree {
  mkdir -p $PAIRTREE_RESULTS_DIR
  logfn=$BASEDIR/scratch/logs/${BATCH}.$(date '+%s').log
  rm -f $logfn

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    resultsfn="$PAIRTREE_RESULTS_DIR/$runid.results.npz"
    is_run_complete $resultsfn && continue
    #[[ $runid =~ K30 || $runid =~ K100 ]] || continue

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
      echo "cd $PAIRTREE_RESULTS_DIR && " \
        "TIMEFORMAT='%R %U %S'; time (OMP_NUM_THREADS=1 " \
        "PATH=$HOME/tmp/jose/bin:$PATH " \
        "LD_LIBRARY_PATH=$HOME/tmp/jose/bin:$LD_LIBRARY_PATH " \
        "python3 $PROTDIR/basic.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $TREE_CHAINS" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--burnin-per-chain $BURNIN_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--phi-fitter $PHI_FITTER" \
        "--params $PAIRTREE_INPUTS_DIR/${runid}.params.json" \
        "$PAIRTREE_INPUTS_DIR/$runid.ssm" \
        "$resultsfn" \
        ">$runid.stdout" \
        "2>$runid.stderr) 2>$runid.time"
    ) #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done #| grep python3 | parallel -j40 --halt 1 --joblog $logfn
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

        basepath="${PAIRTREE_RESULTS_DIR}/${runid}.pairtree_clustrel"
        echo "cd $PAIRTREE_RESULTS_DIR &&" \
          "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutrels.py" \
          "--clustrel-mutrel ${basepath}.mutrel.npz" \
          "$resultsfn" \
          ">${basepath}.output_conversion.stdout" \
          "2>${basepath}.output_conversion.stderr"

        subset_size=2000
        for use_subset in false; do
          for tree_weights in llh; do
            basepath="${PAIRTREE_RESULTS_DIR}/${runid}.pairtree_trees."
            basepath+=$([[ "$use_subset" == "true" ]] && echo subset || echo all)
            basepath+=".${tree_weights}"

            cmd="cd $PAIRTREE_RESULTS_DIR && "
            cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutphis.py "
            cmd+="--weight-trees-by $tree_weights "
            cmd+="--ssms ${PAIRTREE_INPUTS_DIR}/${runid}.ssm "
            if [[ "$use_subset" == "true" ]]; then
              cmd+="--use-subset $subset_size "
            fi
            cmd+="$resultsfn "
            cmd+="${basepath}.mutphi.npz "
            cmd+=">${basepath}.mutphi_output_conversion.stdout "
            cmd+="2>${basepath}.mutphi_output_conversion.stderr"
            echo $cmd

            cmd="cd $PAIRTREE_RESULTS_DIR && "
            cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutrels.py "
            if [[ "$use_subset" == "true" ]]; then
              cmd+="--use-subset $subset_size "
            fi
            cmd+="--weight-trees-by $tree_weights "
            cmd+="--trees-mutrel ${basepath}.mutrel.npz "
            cmd+="$resultsfn "
            cmd+=">${basepath}.mutrel_output_conversion.stdout "
            cmd+="2>${basepath}.mutrel_output_conversion.stderr"
            echo $cmd
          done
        done
    ) #> $jobfn
    #sbatch $jobfn
    #rm $jobfn
  done
}

function main {
  #run_pairtree
  convert_outputs
}

main
