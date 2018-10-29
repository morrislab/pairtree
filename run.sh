#!/bin/bash
set -euo pipefail

PROTDIR=~/work/pairtree
RESULTSDIR=~/work/pairtree/scratch/results
INPUTSDIR=~/work/pairtree/inputs
HANDBUILTDIR=~/work/steph/data/handbuilt_trees
JOBDIR=$SCRATCH/jobs

PARALLEL=40
TREE_SAMPLES=1000
PHI_ITERATIONS=10000

function make_trees {
  #for treetype in patient xeno; do
  for treetype in xeno; do
    OUTDIR=$RESULTSDIR/$treetype
    mkdir -p $OUTDIR

    #rm -f $OUTDIR/SJ*.{html,json,csv,stdout}
    mkdir -p $OUTDIR

    if [[ $treetype == "xeno" ]]; then
      INDIR=$INPUTSDIR/steph.xenos.nocns
    else
      INDIR=$INPUTSDIR/steph.vanilla
    fi

    for ssmfn in $INDIR/*.ssm; do
      runid=$(basename $ssmfn | cut -d. -f1)
      jobname="${treetype}_${runid}"
      (
        echo "#!/bin/bash"
        echo "#SBATCH --nodes=1"
        echo "#SBATCH --ntasks=$PARALLEL"
        echo "#SBATCH --time=24:00:00"
        echo "#SBATCH --job-name $jobname"
        echo "#SBATCH --output=$JOBDIR/slurm_${jobname}_%j.txt"
        echo "#SBATCH --mail-type=FAIL"

        echo "cd $OUTDIR && " \
          "python3 $PROTDIR/pairtree.py" \
          "--tree-samples $TREE_SAMPLES" \
          "--phi-iterations $PHI_ITERATIONS" \
          "--tree-chains $PARALLEL" \
          "--parallel $PARALLEL" \
            "$runid" \
            "$INDIR/$runid.{ssm,params.json}" \
            "$HANDBUILTDIR/$runid.json handbuilt.$treetype" \
            ">$runid.stdout" \
            "2>$runid.stderr"
      ) #> $JOBDIR/job.sh
      #sbatch $JOBDIR/job.sh
      #rm $JOBDIR/job.sh
    done
  done
}

function write_indices {
  for treetype in patient xeno; do
    OUTDIR=$RESULTSDIR/$treetype
    cd $OUTDIR

    for stdoutfn in *.stdout; do
      runid=$(basename $stdoutfn | cut -d. -f1)
      (
	echo "cluster1,cluster2,distance,norm_distance"
	cat $stdoutfn | grep ^cluster_dist | cut -d, -f2- | sort -nrk4 -t,
      ) > $runid.cluster_distances.csv
      (
	echo "phi,cluster,distance,norm_distance"
	cat $stdoutfn | grep ^phi_dist | cut -d, -f2- | sort -nrk4 -t,
      ) > $runid.phi_distances.csv
    done

    for resultsfn in *.results.html; do
      runid=$(echo $resultsfn | cut -d. -f1)
      echo "<h3>$runid</h3>"
      echo "<ul>"
      echo "<li><a href=$runid.results.html>$runid results</a></li>"
      echo "<li><a href=$runid.cluster_distances.csv>$runid cluster distances</a></li>"
      echo "<li><a href=$runid.phi_distances.csv>$runid phi distances</a></li>"
      echo "</ul>"
    done > index.html
  done
}

function main {
  make_trees
  #write_indices
}

main
