#!/bin/bash
set -euo pipefail

PROTDIR=~/work/pairtree
INDIR=$PROTDIR/scratch/inputs/sims.pairtree
JOBDIR=/tmp

PARALLEL=40
TREES_PER_CHAIN=1000
PHI_ITERATIONS=10000
FRACGARB=20

function make_simulated_data {
  mkdir -p $INDIR
  module load gnu-parallel

  for K in 2 4 8; do
  for S in 1 3 10 30; do
  for T in 50 200 1000; do
  for M in 20 50 100; do
  for run in $(seq 1); do
    G=$(expr $M / $FRACGARB) # Takes floor of dividend
    jobname="sim_K${K}_S${S}_T${T}_M${M}_G${G}_run${run}"
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
        "> $INDIR/$jobname.stdout" \
        "2>$INDIR/$jobname.stderr" \
    )
  done
  done
  done
  done
  done | parallel -j$PARALLEL --halt 1
}

function main {
  make_simulated_data
}

main
