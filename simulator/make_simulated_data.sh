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

  for K in 3 10; do
  for S in 1 3 10 30 100; do
  for T in 50 200 1000; do
  for M_per_cluster in 10 20 100; do
  for G_frac in 0 0.01 0.1 1; do
  for run in $(seq 1); do
    M=$(echo "$K * $M_per_cluster" | bc)
    G=$(echo "(($G_frac * $M) + 0.5) / 1" | bc)

    jobname="sim_K${K}_S${S}_T${T}_M${M}_G${G}_run${run}"
    (
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
  done
  done | parallel -j$PARALLEL --halt 1
}

function main {
  make_simulated_data
}

main
