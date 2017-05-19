#!/bin/sh
set -euo pipefail

function main {
  PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  BASEDIR=~/work/steph
  SSMDIR=$BASEDIR/data/inputs/steph.xenos
  OUTDIR=$BASEDIR/data/pairwise

  for ssmfn in $SSMDIR/*.sampled.ssm; do
    sampid=$(basename $ssmfn | cut -d . -f1)
    echo "python3 $PROTDIR/pairwise.py "\
      "$ssmfn" \
      "$OUTDIR/$sampid.pairwise.json" \
      "> $OUTDIR/$sampid.stdout" \
      "2> $OUTDIR/$sampid.stderr"
  done | parallel -j40 --halt 1 --joblog /tmp/pairwise.jobs.txt

  #for jsonfn in $OUTDIR/*.pairwise.json; do
  #  sampid=$(basename $jsonfn | cut -d . -f1)
  #  cmd="python3 $PROTDIR/plot.py"
  #  echo "$cmd --cluster $sampid $jsonfn $OUTDIR/$sampid.clustered.pairwise.html"
  #  echo "$cmd           $sampid $jsonfn $OUTDIR/$sampid.unclustered.pairwise.html"
  #done | parallel -j40 --halt 1
}

main
