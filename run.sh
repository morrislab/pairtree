#!/bin/sh
set -euo pipefail

PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=~/work/steph
SSMDIR=$BASEDIR/data/inputs/steph.xenos
OUTDIR=$BASEDIR/data/pairwise
RENAMEDSAMPS=$BASEDIR/misc/renamed.txt
HIDDENSAMPS=$BASEDIR/misc/hidden.txt

function rename_samples {
  for paramsfn in $SSMDIR/*.params.json; do
    sampid=$(basename $paramsfn | cut -d . -f1)
    echo "python3 $PROTDIR/rename_samples.py" \
      "$sampid" \
      "$HIDDENSAMPS" \
      "$RENAMEDSAMPS" \
      "$paramsfn"
  done | parallel -j40 --halt 1
}

function calc_pairwise {
  rm -f $OUTDIR/*.{pairwise.json,stdout,stderr}

  for ssmfn in $SSMDIR/*.sampled.ssm; do
    sampid=$(basename $ssmfn | cut -d . -f1)
    echo "python3 $PROTDIR/pairwise.py "\
      "$ssmfn" \
      "$OUTDIR/$sampid.pairwise.json" \
      "> $OUTDIR/$sampid.stdout" \
      "2> $OUTDIR/$sampid.stderr"
  done | parallel -j40 --halt 1
}

function plot {
  rm -f $OUTDIR/*.{pairwise.html,js}

  cp -a $PROTDIR/highlight_table_labels.js $OUTDIR/
  for jsonfn in $OUTDIR/*.pairwise.json; do
    sampid=$(basename $jsonfn | cut -d . -f1)
    ssmfn=$SSMDIR/$sampid.sampled.ssm
    paramsfn=$SSMDIR/$sampid.params.json
    spreadsheetfn=$BASEDIR/data/ssms/$sampid.csv

    for output_type in clustered unclustered condensed; do
      echo "python3 $PROTDIR/plot.py " \
	"--output-type $output_type " \
	"$sampid" \
	"$jsonfn" \
	"$ssmfn" \
	"$paramsfn" \
	"$spreadsheetfn" \
	"$OUTDIR/$sampid.$output_type.pairwise.html"
    done
  done | parallel -j40 --halt 1
}

function write_index {
  cd $OUTDIR
  for status in clustered unclustered condensed; do
    echo "<h3>$status</h3>"
    for htmlfn in S*.$status.pairwise.html; do
      sampid=$(basename $htmlfn | cut -d. -f1)
      echo "<a href=$htmlfn>$sampid</a><br>"
    done
  done > index.html
}

function main {
  mkdir -p $OUTDIR

  #rename_samples
  calc_pairwise
  plot
  write_index
}

main
