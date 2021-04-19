#!/bin/sh
set -euo pipefail
shopt -s nullglob

BASEDIR=~/work/pairtree
PEARSIMDIR=~/work/pearsim
PAIRTREEDIR=$BASEDIR/bin
INDIR=$BASEDIR/scratch/inputs/garbdetect
RESULTSDIR=$BASEDIR/scratch/results/garbdetect
JOBDIR=~/jobs
PARA=40
PYTHON=python3

function commafy {
  echo $(echo $1 | sed 's/\./,/g')
}


function make_inputs {
  mkdir -p $INDIR && cd $INDIR

  for K in 30; do
  for S in 30; do
  for T in 1000; do
  for M_per_cluster in 20; do
  for G_frac in 0.02; do
  for run in $(seq 20); do
  for alpha in 1 0.1; do
  for garbtype in uniform acquired_twice wildtype_backmut missed_cna; do
    M=$(echo "$K * $M_per_cluster" | bc)
    G=$(echo "(($G_frac * $M) + 0.5) / 1" | bc)

    jobname="sim_$(echo $garbtype | sed 's/_//')_alpha$(commafy $alpha)_run${run}"
    (
      echo "$PYTHON $PEARSIMDIR/make_simulated_data.py" \
        "--seed $run" \
        "--write-clusters" \
        "-K $K" \
        "-S $S" \
        "-T $T" \
        "-M $M" \
        "-G $G" \
        "--alpha $alpha" \
        "--garbage-type $garbtype" \
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
  done
  done
  done | parallel -j$PARA --halt 2 --eta
}

function init_pairwise {
  # Initially, produce a file storing the pairwise evidence. Delete the
  # resulting params.json output since we won't use it. I could modify the
  # `removegarbage` program to have a mode for only creating the pairwise
  # evidence tensor without producing params.json, but this option will
  # presumably be used only rarely, and it complicates the workflow. Over the
  # longer term, I should have a separate program for creating the evidence
  # tensor, then reuse that for removing garbage, clustering variants, and
  # building trees. This would require changes throughout Pairtree, however --
  # namely, we'd return to work with a tensor build around individual variants,
  # not supervariants. This would also require a more efficient means of
  # computing the tensor -- e.g., we could pre-compute the tensor for all
  # likely pairwise values, then determine when the binomial recision becomes
  # high enough that we don't need to cache any more values.
  for ssmfn in $INDIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    results=$RESULTSDIR/$runid
    mkdir -p $results && cd $results
    jobfn=$(mktemp)

    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=$PARA\n"
    cmd+="#SBATCH --time=23:59:00\n"
    cmd+="#SBATCH --job-name removegarb_$runid\n"
    cmd+="#SBATCH --output=$JOBDIR/slurm_pairtree_${runid}_%j.txt\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="$PYTHON $PAIRTREEDIR/removegarbage"
    cmd+=" --pairwise-results $results/$runid.pairwise.npz"
    cmd+=" $INDIR/$runid.ssm"
    cmd+=" $INDIR/$runid.params.json"
    cmd+=" $results/$runid.params.json"
    cmd+=" > $results/$runid.stdout"
    cmd+=" 2>$results/$runid.stderr"
    # Remove unused params.json file.
    cmd+=" && rm $results/$runid.params.json"

    echo -e $cmd > $jobfn
    sbatch $jobfn
    rm $jobfn
  done #| parallel -j$PARA --halt 2 --eta
}

function detect_garb {
  for pairwisefn in $RESULTSDIR/sim*/*.pairwise.npz; do
    results=$(dirname $pairwisefn)
    runid_base=$(basename $results)

    for garb_prior in 0.001 0.01 0.1 0.2 0.5; do
      for max_garb in 0.001 0.01 0.1 0.2 0.5; do
        runid="${runid_base}_prior$(commafy $garb_prior)_maxgarb$(commafy $max_garb)"
        cmd=""
        cmd+="$PYTHON $PAIRTREEDIR/removegarbage"
        cmd+=" --prior $garb_prior"
        cmd+=" --max-garb-prob $max_garb"
        cmd+=" --pairwise-results $pairwisefn"
        cmd+=" $INDIR/$runid_base.ssm"
        cmd+=" $INDIR/$runid_base.params.json"
        cmd+=" $results/$runid.params.json"
        cmd+=" > $results/$runid.stdout"
        cmd+=" 2>$results/$runid.stderr"
        echo $cmd
      done
    done
  done
}

function eval_garbdetect {
  for results in $RESULTSDIR/sim*; do
    runid_base=$(basename $results)
    for jsonfn in $results/*.params.json; do
      runid=$(basename $jsonfn | cut -d. -f1)
      cmd="$PYTHON $BASEDIR/comparison/pairtree/eval_garbage.py"
      cmd+=" $INDIR/$runid_base.params.json"
      cmd+=" $results/$runid.params.json"
      cmd+=" >$results/$runid.eval.json"
      echo $cmd
    done
  done #| parallel -j$PARA --halt 2 --eta
}

function main {
  #make_inputs
  init_pairwise
  #detect_garb
  #eval_garbdetect
}

main
