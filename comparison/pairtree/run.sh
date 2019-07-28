#!/bin/bash
set -euo pipefail
command -v parallel > /dev/null || module load gnu-parallel

PROTDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
JOBDIR=~/jobs

PARALLEL=40
#TREE_CHAINS=$PARALLEL
TREE_CHAINS=1
TREES_PER_CHAIN=3000
PHI_ITERATIONS=10000
PHI_FITTER=projection
THINNED_FRAC=1.0

BATCH=sims.pairtree
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/${BATCH}.lol51
#BATCH=steph.pairtree.multichain
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree
#PAIRTREE_RESULTS_DIR=$BASEDIR/scratch/results/$BATCH

function is_run_complete {
  if ([[ -f $1 ]] && unzip -l $1 | grep "adjm.npy" > /dev/null); then
    return 0
  else
    return 1
  fi
}

function is_big_K {
  runid=$1
  if [[ $runid =~ K30 || $runid =~ K100 ]]; then
    return 0
  else
    return 1
  fi
}

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
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $runid"
      echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

        #"--only-build-tensor" \
      echo "mkdir -p $outdir && cd $outdir && " \
        "TIMEFORMAT='%R %U %S'; time (" \
        "OMP_NUM_THREADS=1 " \
        "PATH=$HOME/tmp/jose/bin:$PATH " \
        "LD_LIBRARY_PATH=$HOME/tmp/jose/bin:$LD_LIBRARY_PATH " \
        "python3 $PROTDIR/pairtree.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $TREE_CHAINS" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--phi-fitter $PHI_FITTER" \
        "--params $PAIRTREE_INPUTS_DIR/${runid}.params.json" \
        "--thinned-frac $THINNED_FRAC" \
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
  for resultsfn in $PAIRTREE_RESULTS_DIR/*/*.results.npz; do
    outdir=$(dirname $resultsfn)
    runid=$(basename $resultsfn | cut -d. -f1)
    is_run_complete $resultsfn || continue

    jobfn=$(mktemp)
    (
        #echo "#!/bin/bash"
        #echo "#SBATCH --nodes=1"
        #echo "#SBATCH --ntasks=$PARALLEL"
        #echo "#SBATCH --time=23:59:00"
        #echo "#SBATCH --job-name convert_$runid"
        #echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
        #echo "#SBATCH --mail-type=NONE"

        #basepath="${outdir}/${runid}.pairtree_clustrel"
        #echo "cd $outdir &&" \
        #  "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutrels.py" \
        #  "--clustrel-mutrel ${basepath}.mutrel.npz" \
        #  "$resultsfn" \
        #  ">${basepath}.output_conversion.stdout" \
        #  "2>${basepath}.output_conversion.stderr"

        subset_size=2000
        for use_subset in false; do
          for tree_weights in llh; do
            basepath="${outdir}/${runid}.pairtree_trees."
            basepath+=$([[ "$use_subset" == "true" ]] && echo subset || echo all)
            basepath+=".${tree_weights}"

            cmd="cd $outdir && "
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

            if ! is_big_K $runid; then
              cmd="cd $outdir && "
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
            fi
          done
        done
    ) #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done
}

function main {
  run_pairtree | grep python3 | parallel -j80 --halt 1 --eta
  convert_outputs | grep python3 | parallel -j10 --halt 1 --eta
}

main
