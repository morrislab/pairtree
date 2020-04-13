#!/bin/sh
set -euo pipefail

CLUSTERS=$SCRATCH/tmp/clusters
BASEDIR=$HOME/work/pairtree
INPUTSDIR=$HOME/work/pairtree-experiments/inputs/sims.smallalpha.pairtree

function run_clustering {
  for model in pants socks; do
    clusterdir=$CLUSTERS/$model
    mkdir -p $clusterdir

    for ssmfn in $INPUTSDIR/sim_K{3,10}_*.ssm; do
      sampid=$(basename $ssmfn | cut -d. -f1)
      echo "echo $sampid \$(NUMBA_DISABLE_JIT=0 python3 $BASEDIR/bin/clustervars --model $model $INPUTSDIR/$sampid.{ssm,params.json} $clusterdir/$sampid.clustered.params.json)"
    done | parallel -j40 --halt 1 --eta > /dev/null
  done
}

function eval_clustering {
  for model in pants socks; do
    clusterdir=$CLUSTERS/$model
    (
      echo "sampid,true_clusters,found_clusters,true_llh,found_llh,true_nlglh,found_nlglh,homogeneity,completeness,vmeasure,ami"
      for paramsfn in $clusterdir/*.params.json; do
        sampid=$(basename $paramsfn | cut -d. -f1)
        # Disable JIT to speed this up.
        echo "echo $sampid,\$(NUMBA_DISABLE_JIT=1 python3 $BASEDIR/comparison/pairtree/compare_clusterings.py $INPUTSDIR/$sampid.ssm $INPUTSDIR/$sampid.params.json $paramsfn)"
      done | parallel -j40 --halt 1 --eta 
    ) > $CLUSTERS/$model.txt

    python3 $BASEDIR/comparison/pairtree/plot_clustering.py $CLUSTERS/$model.{txt,html}
  done
}

function main {
  run_clustering
  eval_clustering
}

main
