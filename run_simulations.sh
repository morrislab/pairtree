#!/bin/bash
set -euo pipefail

RESULTSDIR=~/work/pairtree/scratch/sims
PROTDIR=~/work/pairtree
JOBDIR=$SCRATCH/jobs

PARALLEL=40

function main {
  mkdir -p $RESULTSDIR

  M=50
  G=10
  for K in 2 4 8; do
    for S in 1 3 10 30; do
      for T in 50 200 1000; do
        for run in $(seq 3); do
          jobname="sim_K${K}_S${S}_T${T}_M${M}_G${G}_run${run}"
          (
            echo "#!/bin/bash"
            echo "#SBATCH --nodes=1"
            echo "#SBATCH --ntasks=$PARALLEL"
            echo "#SBATCH --time=24:00:00"
            echo "#SBATCH --job-name $jobname"
            echo "#SBATCH --output=$JOBDIR/${jobname}.%j.log"
            echo "#SBATCH --mail-type=FAIL"
            echo "python3 $PROTDIR/run_simulation.py" \
              "--parallel $PARALLEL" \
              "-K $K" \
              "-S $S" \
              "-T $T" \
              "-M $M" \
              "-G $G" \
              "$RESULTSDIR/$jobname.data.pickle" \
              "$RESULTSDIR/$jobname.pairs.npz" \
              "$RESULTSDIR/$jobname.results.html" \
              "> $JOBDIR/$jobname.stdout" \
              "2>$JOBDIR/$jobname.stderr" \
          ) > $JOBDIR/job.sh
          sbatch $JOBDIR/job.sh
          rm $JOBDIR/job.sh
        done
      done
    done
  done
}

main
