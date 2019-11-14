#!/bin/bash
set -euo pipefail

PROTDIR=~/work/pairtree
SIMDIR=~/work/pearsim
INDIR=$PROTDIR/scratch/inputs/sims.smallalpha.pairtree
JOBDIR=/tmp

PARALLEL=40

function make_simulated_data {
  mkdir -p $INDIR
  module load gnu-parallel

  for K in 3 10 30 100; do
  for S in 1 3 10 30 100; do
    (( $K >= 30 && $S < 10 )) && continue
  for T in 50 200 1000; do
  for M_per_cluster in 10 20 100; do
  for G_frac in 0; do
  for run in $(seq 1 4); do
    M=$(echo "$K * $M_per_cluster" | bc)
    G=$(echo "(($G_frac * $M) + 0.5) / 1" | bc)

    jobname="sim_K${K}_S${S}_T${T}_M${M}_G${G}_run${run}"
    (
      echo "python3 $SIMDIR/make_simulated_data.py" \
        "--write-clusters" \
        "--write-numpy $INDIR/$jobname.truth.npz" \
        "-K $K" \
        "-S $S" \
        "-T $T" \
        "-M $M" \
        "-G $G" \
        "--alpha 0.1" \
        "$INDIR/$jobname.truth.pickle" \
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
