#!/bin/sh
set -euo pipefail

PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=~/work/steph
SSMDIR=$BASEDIR/data/inputs/steph.xenos.nocns
RUN=xeno
SCINPUTSDIR=$BASEDIR/data/sciclone/inputs/$RUN
SCRESULTSDIR=$BASEDIR/data/sciclone/results/$RUN
HANDBUILTDIR=$BASEDIR/data/handbuilt_trees

function makedirs {
  mkdir -p $SCINPUTSDIR $SCRESULTSDIR
}

function convert_inputs {
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

function run_sciclone {
  cd $SCINPUTSDIR
  for sampid in *; do
    rm -f "$SCRESULTSDIR/$sampid.results.txt"
    echo "Rscript --vanilla $PROTDIR/run_sciclone.R" \
      "$PWD/$sampid/*.dat" \
      "$SCRESULTSDIR/$sampid.results.txt" \
      "$SCRESULTSDIR/$sampid.sampnames.json" \
      ">  $SCRESULTSDIR/$sampid.stdout" \
      "2> $SCRESULTSDIR/$sampid.stderr"
  done | parallel -j40 --halt 1
}

function convert_outputs {
  cd $SCRESULTSDIR
  for resultsfn in *.results.txt; do
    sampid=$(basename $resultsfn | cut -d. -f1)
    python3 "$PROTDIR/convert_sciclone_results.py" \
      "$sampid" \
      "$SSMDIR/$sampid.sampled.ssm" \
      "$resultsfn" \
      "$sampid.handbuilt.json"
  done
}

function plot {
  rm -f $SCRESULTSDIR/*.{pairwise.html,css,js}

  cp -a $PROTDIR/*.{css,js} $SCRESULTSDIR/
  cd $SCRESULTSDIR
  for handbuiltfn in $SCRESULTSDIR/*.handbuilt.json; do
    sampid=$(basename $handbuiltfn | cut -d . -f1)
    jsonfn=
    ssmfn=$SSMDIR/$sampid.sampled.ssm
    paramsfn=$SSMDIR/$sampid.params.json
    spreadsheetfn=$BASEDIR/data/ssms/$sampid.csv
    pairwisefn=$BASEDIR/data/pairwise.$RUN.nocns/$sampid.pairwise.json
    orig_handbuiltfn=$HANDBUILTDIR/$sampid.json

    cmd="python3 $PROTDIR/plot.py
      --tree-type handbuilt.$RUN
      $sampid
      $pairwisefn
      $ssmfn
      $paramsfn
      $spreadsheetfn
      $handbuiltfn
      $SCRESULTSDIR/$sampid.pairwise.html
      $SCRESULTSDIR/$sampid.summ.json
      $SCRESULTSDIR/$sampid.muts.json
      $SCRESULTSDIR/$sampid.phi.json
      $SCRESULTSDIR/$sampid.clustermat.json
      >  $SCRESULTSDIR/$sampid.plot.stdout
      2> $SCRESULTSDIR/$sampid.plot.stderr"
    if [[ -f $orig_handbuiltfn ]]; then
      cmd+=" && python3 $PROTDIR/clustermat.py
	--sampid $sampid
	--ssms $ssmfn
	--params $paramsfn
	--handbuilt1 $orig_handbuiltfn
	--handbuilt2 $handbuiltfn
	--treetype1 handbuilt.xeno
	--treetype2 handbuilt.xeno
	$SCRESULTSDIR/$sampid.sciclone_vs_handbuilt.clustermat.json
	>> $SCRESULTSDIR/$sampid.pairwise.html"
    fi
    echo $cmd
  done | parallel -j40 --halt 1
}

function main {
  makedirs
  convert_inputs
  run_sciclone
  convert_outputs
  plot
}

main
