#!/bin/sh
set -euo pipefail

function main {
  PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  BASEDIR=~/work/steph
  SSMDIR=$BASEDIR/data/inputs/steph.xenos
  OUTDIR=$BASEDIR/data/pairwise

  mkdir -p $OUTDIR
  rm -f $OUTDIR/*.{pairwise.json,stdout,stderr,js}
  rm -f $OUTDIR/*.pairwise.html

  for ssmfn in $SSMDIR/*.sampled.ssm; do
    sampid=$(basename $ssmfn | cut -d . -f1)
    echo "python3 $PROTDIR/pairwise.py "\
      "$ssmfn" \
      "$OUTDIR/$sampid.pairwise.json" \
      "> $OUTDIR/$sampid.stdout" \
      "2> $OUTDIR/$sampid.stderr"
  done | parallel -j40 --halt 1

  cp -a $PROTDIR/highlight_table_labels.js $OUTDIR/
  for jsonfn in $OUTDIR/*.pairwise.json; do
    sampid=$(basename $jsonfn | cut -d . -f1)
    cmd="python3 $PROTDIR/plot.py"
    ssmfn=$SSMDIR/$sampid.sampled.ssm
    paramsfn=$SSMDIR/$sampid.params.json
    spreadsheetfn=$BASEDIR/data/ssms/$sampid.csv

    echo "$cmd --cluster $sampid $jsonfn $ssmfn $paramsfn $spreadsheetfn $OUTDIR/$sampid.clustered.pairwise.html"
    echo "$cmd           $sampid $jsonfn $ssmfn $paramsfn $spreadsheetfn $OUTDIR/$sampid.unclustered.pairwise.html"
  done | parallel -j40 --halt 1

  cd $OUTDIR
  for status in clustered unclustered; do
    for htmlfn in S*.$status.pairwise.html; do
      sampid=$(basename $htmlfn | cut -d. -f1)
      echo "<a href=$htmlfn>$sampid ($status)</a><br>"
    done
    echo "<br>"
  done > index.html
}

main
