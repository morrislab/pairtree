#!/bin/sh
set -euo pipefail

PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=~/work/steph
SSMDIR=$BASEDIR/data/inputs/steph.xenos.nocns
SCINPUTSDIR=$BASEDIR/data/sciclone/inputs

RUN=$1
HANDBUILTDIR=$BASEDIR/data/handbuilt_trees
SCRESULTSDIR=$BASEDIR/data/sciclone/results/$RUN

function run_sciclone {
  mkdir -p $SCRESULTSDIR

  cd $SCINPUTSDIR
  for sampid in *; do
    rm -f "$SCRESULTSDIR/$sampid.results.txt"
    if [[ $RUN == "xeno" ]]; then
      inputs="$PWD/$sampid/*.dat"
    else
      inputs=$(for foo in "$PWD/$sampid/"*.dat; do
	if [[ $foo =~ "Xeno" ]]; then
	  continue
	fi
	echo \"$foo\"
      done | paste -sd' ')
    fi

    echo "Rscript --vanilla $PROTDIR/run_sciclone.R" \
      "$inputs" \
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
      "$RUN" \
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
	--treetype1 handbuilt.$RUN
	--treetype2 handbuilt.$RUN
	$SCRESULTSDIR/$sampid.sciclone_vs_handbuilt.clustermat.json
	>> $SCRESULTSDIR/$sampid.pairwise.html"
    fi
    echo $cmd
  done | parallel -j40 --halt 1
}

function main {
  #run_sciclone
  convert_outputs
  #plot
}

main
