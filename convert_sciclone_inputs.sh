#!/bin/sh
set -euo pipefail

PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=~/work/steph
SSMDIR=$BASEDIR/data/inputs/steph.xenos.nocns
SCINPUTSDIR=$BASEDIR/data/sciclone/inputs

function convert_inputs {
  mkdir -p $SCINPUTSDIR

  cd $SSMDIR
  for ssmfn in *.sampled.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    paramsfn=$SSMDIR/$sampid.params.json
    out=$SCINPUTSDIR/$sampid
    mkdir -p $out
    rm -f $out/*.dat
    python3 "$PROTDIR/convert_to_sciclone.py" \
      "$sampid" \
      "$ssmfn" \
      "$paramsfn" \
      "$out"
  done
}

convert_inputs
