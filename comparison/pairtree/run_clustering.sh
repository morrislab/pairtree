#!/bin/sh
set -euo pipefail

CLUSTERS=$SCRATCH/tmp/clusters
BASEDIR=$HOME/work/pairtree
INPUTSDIR=$HOME/work/pairtree-experiments/inputs/sims.smallalpha.pairtree

function run_clustering {
  for model in socks; do
    clusterdir=$CLUSTERS/$model
    mkdir -p $clusterdir

    for ssmfn in $INPUTSDIR/sim_K{30,100}_*.ssm; do
      runid=$(basename $ssmfn | cut -d. -f1)
      echo "echo $runid \$(NUMBA_DISABLE_JIT=0 python3 $BASEDIR/bin/clustervars --model $model $INPUTSDIR/$runid.{ssm,params.json} $clusterdir/$runid.clustered.params.json)"
    done
  done | parallel -j40 --halt 1 --eta > /dev/null
}

function eval_clustering {
  for model in pants socks; do
    clusterdir=$CLUSTERS/$model
    (
      echo "runid,true_clusters,found_clusters,true_llh,found_llh,true_nlglh,found_nlglh,homogeneity,completeness,vmeasure,ami"
      for paramsfn in $clusterdir/*.params.json; do
        runid=$(basename $paramsfn | cut -d. -f1)
        # Disable JIT to speed this up.
        echo "echo $runid,\$(NUMBA_DISABLE_JIT=1 python3 $BASEDIR/comparison/pairtree/compare_clusterings.py $INPUTSDIR/$runid.ssm $INPUTSDIR/$runid.params.json $paramsfn)"
      done | parallel -j40 --halt 1 --eta 
    ) > $CLUSTERS/$model.txt

    python3 $BASEDIR/comparison/pairtree/plot_clustering.py $CLUSTERS/$model.{txt,html}
  done
}

function main {
  #run_clustering
  eval_clustering
}

main
