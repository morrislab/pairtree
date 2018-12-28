#!/bin/bash
set -euo pipefail

INDIR=~/work/pairtree/scratch/inputs/sims
PROTDIR=~/work/pairtree
JOBDIR=/tmp

PARALLEL=4
TREES_PER_CHAIN=1000
PHI_ITERATIONS=10000

function make_simulated_data {
  mkdir -p $INDIR

  M=10
  G=1
  for K in 2 4 8; do
    for S in 1 3 10 30; do
      for T in 50 200 1000; do
        for run in $(seq 3); do
          jobname="sim_K${K}_S${S}_T${T}_M${M}_G${G}_run${run}"
          jobfn=$(mktemp)
          (
            #echo "#!/bin/bash"
            #echo "#SBATCH --nodes=1"
            #echo "#SBATCH --ntasks=$PARALLEL"
            #echo "#SBATCH --time=24:00:00"
            #echo "#SBATCH --job-name $jobname"
            #echo "#SBATCH --output=$JOBDIR/${jobname}.%j.log"
            #echo "#SBATCH --mail-type=FAIL"
            echo "python3 $PROTDIR/make_simulated_data.py" \
              "--write-clusters" \
              "-K $K" \
              "-S $S" \
              "-T $T" \
              "-M $M" \
              "-G $G" \
              "$INDIR/$jobname.data.pickle" \
              "$INDIR/$jobname.params.json" \
              "$INDIR/$jobname.ssm" \
              "> $JOBDIR/$jobname.stdout" \
              "2>$JOBDIR/$jobname.stderr" \
          ) #> $jobfn
          #sbatch $jobfn
          rm $jobfn
        done
      done
    done
  done
}

function run_pairtree {
  OUTDIR=~/work/pairtree/scratch/results/sims

  for ssmfn in $INDIR/*.ssm; do
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

      echo "cd $OUTDIR && " \
        "OMP_NUM_THREADS=1 python3 $PROTDIR/basic.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $PARALLEL" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--params $INDIR/${runid}.params.json" \
        "$INDIR/$runid.ssm" \
        "$OUTDIR/$runid.results.npz" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done
}

function main {
  #make_simulated_data
  run_pairtree
}

main
