#!/bin/bash
set -euo pipefail

INPUTSDIR=~/work/pairtree/scratch/inputs/sims
PROTDIR=~/work/pairtree
JOBDIR=$SCRATCH/jobs

function make_simulated_data {
  mkdir -p $INPUTSDIR

  M=100
  G=10
  for K in 2 4 8; do
    for S in 1 3 10 30; do
      for T in 50 200 1000; do
        for run in $(seq 3); do
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
              "-K $K" \
              "-S $S" \
              "-T $T" \
              "-M $M" \
              "-G $G" \
              "$INPUTSDIR/$jobname.data.pickle" \
              "$INPUTSDIR/$jobname.ssm" \
              "> $JOBDIR/$jobname.stdout" \
              "2>$JOBDIR/$jobname.stderr" \
          ) #> $JOBDIR/job.sh
          #sbatch $JOBDIR/job.sh
          #rm $JOBDIR/job.sh
        done
      done
    done
  done
}

function main {
  make_simulated_data
}

main
