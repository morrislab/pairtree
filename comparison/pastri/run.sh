#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
PARALLEL=40
NUM_ITERS=10000
PASTRI_DIR=$HOME/.apps/pastri

#BATCH=sims.pastri.informative
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
BATCH=steph.xeno.pastri.hbclusters
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
      #"--uniform-proposal" \
    echo "python3 $SCRIPTDIR/convert_inputs.py " \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.counts" \
      "$INDIR/$sampid.proposal"
  done #| parallel -j$PARALLEL --halt 1
}

function run_pastri {
  mkdir -p $OUTDIR

  for countsfn in $INDIR/*.counts; do
    runid=$(basename $countsfn | cut -d. -f1)
    jobname="steph_pastri_${runid}"

    [[ -f $OUTDIR/$runid.trees ]] && continue

    (
      # Must source ~/.bash_host_specific to get PATH set properly for
      # Miniconda.
      echo "source $HOME/.bash_host_specific && " \
        "cd $OUTDIR && " \
        "python2 $PASTRI_DIR/src/RunPASTRI.py" \
        "--output_prefix $runid" \
        "--num_iters $NUM_ITERS" \
        "$INDIR/${runid}.counts" \
        "$INDIR/${runid}.proposal" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) 
  done #| parallel -j$PARALLEL --joblog $SCRATCH/tmp/$BATCH.log
}

function get_F_and_C {
  for treesfn in $OUTDIR/*.trees; do
    runid=$(basename $treesfn | cut -d. -f1)
    # Trees in $treesfn are ranked by likelihood. `get_F_and_C.py` takes tree
    # rank as a parameter. Thus, if we count N trees with LH > 0, we know their
    # ranks are [0, 1, ..., N-1].
    valid_count=$(cat $treesfn | grep '^>' | grep -v -- '-inf$' | wc -l)
    for idx in $(seq $valid_count); do
      echo "source $HOME/.bash_host_specific && " \
        "cd $OUTDIR && " \
        "python2 $PASTRI_DIR/src/get_F_and_C.py" \
        "-i $idx" \
        "-o $OUTDIR/$runid" \
        "$INDIR/${runid}.counts" \
        "$treesfn" \
        "$OUTDIR/${runid}.fsamples" \
        "> $OUTDIR/${runid}.${idx}.output_conversion.stdout" \
        "2>$OUTDIR/${runid}.${idx}.output_conversion.stderr"
    done
  done | parallel -j$PARALLEL
}

function convert_outputs {
  for tree_weights in llh uniform; do
    for treesfn in $OUTDIR/*.trees; do
      runid=$(basename $treesfn | cut -d. -f1)
      echo "source $HOME/.bash_host_specific && " \
        "cd $OUTDIR && " \
        "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutrels.py" \
        "--weight-trees-by $tree_weights" \
        "$runid" \
        "$PAIRTREE_INPUTS_DIR/$runid.params.json" \
        "$treesfn" \
        "${OUTDIR}/${runid}.pastri_trees_${tree_weights}.mutrel.npz" #\
        #"> $OUTDIR/${runid}.mutrel_output_conversion.stdout" \
        #"2>$OUTDIR/${runid}.mutrel_output_conversion.stderr"
      echo "source $HOME/.bash_host_specific && " \
        "cd $OUTDIR && " \
        "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_mutphis.py" \
        "--weight-trees-by $tree_weights" \
        "$runid" \
        "${PAIRTREE_INPUTS_DIR}/${runid}.ssm" \
        "${PAIRTREE_INPUTS_DIR}/${runid}.params.json" \
        "$treesfn" \
        "${OUTDIR}/${runid}.pastri_trees_${tree_weights}.mutphi.npz" #\
        #"> $OUTDIR/${runid}.mutphi_output_conversion.stdout" \
        #"2>$OUTDIR/${runid}.mutphi_output_conversion.stderr"
    done
  done #| parallel -j$PARALLEL
}

function main {
  #convert_inputs
  #run_pastri
  get_F_and_C
  convert_outputs
}

main
