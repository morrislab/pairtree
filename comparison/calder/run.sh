#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel
command -v java > /dev/null     || module load java

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
PARA=4
CALDER_DIR=$HOME/.apps/calder
GUROBI_DIR=$HOME/.apps/gurobi801
JOB_DIR=$HOME/jobs
PYTHON=$HOME/.apps/miniconda3/bin/python3

BATCH=sims.smallalpha.calder
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.calder

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTDIR=$BASEDIR/scratch/results/$BATCH

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "$PYTHON3 $SCRIPTDIR/convert_inputs.py " \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$PAIRTREE_INPUTS_DIR/$sampid.calder"
  done | parallel -j$PARA --halt 2
}

function run_calder {
  module load java
  for infn in $INDIR/*.calder; do
    runid=$(basename $infn | cut -d. -f1)
    outpath=$OUTDIR/$runid
    mkdir -p $outpath

    cmd="GRB_LICENSE_FILE=$GUROBI_DIR/gurobi.lic LD_LIBRARY_PATH=$GUROBI_DIR/linux64/lib"
    cmd+=" java -jar $CALDER_DIR/calder.jar"
    cmd+=" -i $infn"
    cmd+=" -o $outpath"
    echo $cmd
  done | parallel -j$PARA --halt 2
}

function main {
  convert_inputs
  run_calder
}

main
