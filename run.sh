#!/bin/bash
#set -euo pipefail

PROTDIR=~/work/pairtree
RESULTSDIR=$PROTDIR/scratch/results
HANDBUILTDIR=~/work/steph/data/handbuilt_trees
JOBDIR=$SCRATCH/jobs

STEPH_INDIR=$PROTDIR/scratch/inputs/steph.xeno.withgarb.pairtree
STEPH_OUTDIR=$RESULTSDIR/steph.xeno.withgarb.pairtree

PARALLEL=40
TREES_PER_CHAIN=1000
PHI_ITERATIONS=10000

function run_sims {
  INDIR=$PROTDIR/scratch/inputs/sims
  OUTDIR=$RESULTSDIR/sims

  for ssmfn in $INDIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $runid"
      echo "#SBATCH --output=$JOBDIR/slurm_${runid}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      echo "cd $OUTDIR && " \
        "python3 $PROTDIR/basic.py" \
        "--parallel $PARALLEL" \
          "$runid" \
          "$INDIR/$runid.ssm" \
          ">$runid.stdout" \
          "2>$runid.stderr"
    ) #> $JOBDIR/job.sh
    #sbatch $JOBDIR/job.sh
    #rm $JOBDIR/job.sh
  done
}

function run_steph {
  mkdir -p $STEPH_OUTDIR
  #rm -f $STEPH_OUTDIR/SJ*.{html,json,csv,stdout}

  for ssmfn in $STEPH_INDIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    jobname="steph_pairtree_${runid}"
    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $jobname"
      echo "#SBATCH --output=$JOBDIR/slurm_${jobname}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      echo "cd $STEPH_OUTDIR && " \
        "python3 $PROTDIR/basic.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $PARALLEL" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--params $STEPH_INDIR/${runid}.params.json" \
        "$STEPH_INDIR/$runid.ssm" \
        "$STEPH_OUTDIR/$runid.results.npz" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function plot_steph {
  module load gnu-parallel

  for resultsfn in $STEPH_OUTDIR/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    (
      echo "cd $STEPH_OUTDIR && " \
        "python3 $PROTDIR/plotter.py" \
        "$runid" \
        "$STEPH_INDIR/$runid.ssm" \
        "$STEPH_INDIR/${runid}.params.json" \
        "$STEPH_OUTDIR/$runid.results.npz" \
        "$STEPH_OUTDIR/$runid.results.html" \
        ">>$runid.stdout" \
        "2>>$runid.stderr"
    )
  done | parallel -j$PARALLEL
}

function write_indices {
  cd $STEPH_OUTDIR

  #for stdoutfn in *.stdout; do
  #  runid=$(basename $stdoutfn | cut -d. -f1)
  #  (
  #    echo "cluster1,cluster2,distance,norm_distance"
  #    cat $stdoutfn | grep ^cluster_dist | cut -d, -f2- | sort -nrk4 -t,
  #  ) > $runid.cluster_distances.csv
  #  (
  #    echo "phi,cluster,distance,norm_distance"
  #    cat $stdoutfn | grep ^phi_dist | cut -d, -f2- | sort -nrk4 -t,
  #  ) > $runid.phi_distances.csv
  #done

  for resultsfn in *.results.html; do
    runid=$(echo $resultsfn | cut -d. -f1)
    echo "<h3>$runid</h3>"
    echo "<ul>"
    echo "<li><a href=$runid.results.html>$runid results</a></li>"
    #echo "<li><a href=$runid.cluster_distances.csv>$runid cluster distances</a></li>"
    #echo "<li><a href=$runid.phi_distances.csv>$runid phi distances</a></li>"
    echo "</ul>"
  done > index.html
}

function main {
  #run_steph
  plot_steph
  write_indices
  #run_sims
}

main
