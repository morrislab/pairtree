#!/bin/sh
set -euo pipefail
shopt -s nullglob

BASEDIR=~/work/pairtree
PEARSIMDIR=~/work/pearsim
PAIRTREEDIR=$BASEDIR/bin
JOBDIR=~/jobs
PYTHON=python3

S=30
T=1000
M_per_cluster=20
G_per_cluster=2
ALPHA=0.1
PARA=80
INDIR=$BASEDIR/scratch/inputs/garbdetect
RESULTSDIR=$BASEDIR/scratch/results/garbdetect

function commafy {
  echo $(echo $1 | sed 's/\./,/g')
}

function make_inputs {
  mkdir -p $INDIR && cd $INDIR

  for K in 10 30; do
  M=$(echo "$K * $M_per_cluster" | bc)
  G=$(echo "$K * $G_per_cluster" | bc)
  for MIN_GARB_PHI_DELTA in 0.0005 0.001 0.01 0.05 0.1; do
  for garbtype in wildtype_backmut missed_cna acquired_twice uniform;  do
  for run in $(seq 20); do
    jobname="sim_$(echo $garbtype | sed 's/_//')_K${K}_mindelta$(commafy $MIN_GARB_PHI_DELTA)_run${run}"
    echo "$PYTHON $PEARSIMDIR/make_simulated_data.py" \
      "--seed $run" \
      "--write-clusters" \
      "-K $K" \
      "-S $S" \
      "-T $T" \
      "-M $M" \
      "-G $G" \
      "--alpha $ALPHA" \
      "--garbage-type $garbtype" \
      "--min-garb-phi-delta $MIN_GARB_PHI_DELTA" \
      "--min-garb-pairs 3" \
      "--min-garb-samps 3" \
      "$INDIR/$jobname.truth.pickle" \
      "$INDIR/$jobname.params.json" \
      "$INDIR/$jobname.ssm" \
      "> $INDIR/$jobname.stdout" \
      "2>$INDIR/$jobname.stderr"
  done
  done
  done
  done
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
    pairwisefn=$results/$runid.pairwise.npz

    #[[ $runid == sim_wildtypebackmut_alpha0,1_run1 ]] || continue
    [[ -f $pairwisefn ]] && continue
    #[[ $(basename $RESULTSDIR) =~ "phidelta0,05" ]] || continue

    mkdir -p $results && cd $results
    jobfn=$(mktemp)
    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --qos=nopreemption\n"
    cmd+="#SBATCH --partition=cpu\n"
    cmd+="#SBATCH --mem=4GB\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=$PARA\n"
    cmd+="#SBATCH --time=23:59:00\n"
    cmd+="#SBATCH --job-name removegarb_$runid\n"
    cmd+="#SBATCH --output=$JOBDIR/slurm_pairtree_${runid}_%j.txt\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="$PYTHON $PAIRTREEDIR/removegarbage"
    cmd+=" --parallel $PARA"
    cmd+=" --pairwise-results $pairwisefn"
    cmd+=" $INDIR/$runid.ssm"
    cmd+=" $INDIR/$runid.params.json"
    cmd+=" $results/$runid.params.json"
    cmd+=" > $results/$runid.stdout"
    cmd+=" 2>$results/$runid.stderr"
    # Remove unused params.json file.
    cmd+=" && rm $results/$runid.params.json"

    echo -e $cmd #> $jobfn
    #sbatch $jobfn
    rm $jobfn
  done
}

function detect_garb {
  for pairwisefn in $RESULTSDIR/sim*/*.pairwise.npz; do
    results=$(dirname $pairwisefn)
    runid_base=$(basename $results)

    for garb_prior in 0.01 0.1 0.2 0.5; do
      for max_garb in 0.001 0.01 0.1 0.5 0.75; do
        runid="${runid_base}_prior$(commafy $garb_prior)_maxgarb$(commafy $max_garb)"
        cmd=""
        cmd+="$PYTHON $PAIRTREEDIR/removegarbage"
        cmd+=" --prior $garb_prior"
        cmd+=" --max-garb-prob $max_garb"
        cmd+=" --pairwise-results $pairwisefn"
        cmd+=" --verbose"
        cmd+=" $INDIR/$runid_base.ssm"
        cmd+=" $INDIR/$runid_base.params.json"
        cmd+=" $results/$runid.params.json"
        cmd+=" > $results/$runid.stdout"
        cmd+=" 2>$results/$runid.stderr"
        echo $cmd
      done
    done
  done | parallel -j$PARA --halt 2 --eta
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
  done | parallel -j$PARA --halt 2 --eta
}

function plot_garbdetect {
  cd $RESULTSDIR
  echo "$PYTHON $BASEDIR/comparison/pairtree/plot_garbage.py" \
    "--plot-fn $RESULTSDIR/garb.html" \
    "--plot-bars" \
    "--plot-2dhist" \
    "--filter-maxgarb 0.01" \
    "--filter-prior 0.2" \
    "--filter-mindelta 0.01" \
    "sim*/*.eval.json" \
    "> $RESULTSDIR/garb.json" | parallel -j$PARA --halt 2 --eta
}

function main {
  #make_inputs
  #init_pairwise
  #detect_garb
  #eval_garbdetect
  plot_garbdetect
}

main
