#!/bin/sh
set -euo pipefail

#RUNNAME=patient
RUNNAME=$1

PROTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=~/work/steph
SSMDIR=$BASEDIR/data/inputs/steph.xenos.nocns
OUTDIR=$BASEDIR/data/pairwise.$RUNNAME.nocns
HANDBUILTDIR=$BASEDIR/data/handbuilt_trees
RENAMEDSAMPS=$BASEDIR/misc/renamed.txt
HIDDENSAMPS=$BASEDIR/misc/hidden.txt
PWGSDIR=~/.apps/phylowgs

OUTPUT_TYPES="clustered"

function remove_samples {
  for paramsfn in $SSMDIR/*.params.json; do
    sampid=$(basename $paramsfn | cut -d . -f1)
    echo "python3 $PROTDIR/remove_samples.py" \
      "$sampid" \
      "$SSMDIR/$sampid.sampled.ssm" \
      "$paramsfn"
  done | parallel -j40 --halt 1
}

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
      "$sampid" \
      "$ssmfn" \
      "$OUTDIR/$sampid.pairwise.json" \
      "> $OUTDIR/$sampid.stdout" \
      "2> $OUTDIR/$sampid.stderr"
  done | parallel -j40 --halt 1
}

function plot {
  rm -f $OUTDIR/*.{pairwise.html,css,js}

  cp -a $PROTDIR/*.{css,js} $OUTDIR/
  for jsonfn in $OUTDIR/*.pairwise.json; do
    sampid=$(basename $jsonfn | cut -d . -f1)
    ssmfn=$SSMDIR/$sampid.sampled.ssm
    paramsfn=$SSMDIR/$sampid.params.json
    spreadsheetfn=$BASEDIR/data/ssms/$sampid.csv
    handbuiltfn="$HANDBUILTDIR/$sampid.json"

    [ -f "$handbuiltfn" ] || continue
    cp -a "$handbuiltfn" "$OUTDIR/$sampid.handbuilt.json"

    for output_type in $OUTPUT_TYPES; do
      echo "python3 $PROTDIR/plot.py " \
	"--output-type $output_type " \
	"--tree-type handbuilt.$RUNNAME " \
	"$sampid" \
	"$jsonfn" \
	"$ssmfn" \
	"$paramsfn" \
	"$spreadsheetfn" \
	"$handbuiltfn" \
	"$OUTDIR/$sampid.$output_type.pairwise.html" \
	"$OUTDIR/$sampid.summ.json" \
	"$OUTDIR/$sampid.muts.json" \
	"$OUTDIR/$sampid.phi.json" \
	">  $OUTDIR/$sampid.plot.stdout" \
	"2> $OUTDIR/$sampid.plot.stderr"
    done
  done | parallel -j40 --halt 1
}

function write_plot_index {
  cd $OUTDIR
  for status in $OUTPUT_TYPES; do
    echo "<h3>$status</h3>"
    for htmlfn in S*.$status.pairwise.html; do
      sampid=$(basename $htmlfn | cut -d. -f1)
      echo "<a href=$htmlfn>$sampid</a><br>"
    done
  done > index.html
}

function add_tree_indices {
  for jsonfn in $OUTDIR/*.summ.json; do
    sampid=$(basename $jsonfn | cut -d . -f1)
    gzip --force "$OUTDIR/$sampid.summ.json" "$OUTDIR/$sampid.muts.json"
    echo "PYTHONPATH=$PWGSDIR python2 $PROTDIR/add_tree_indices.py" \
      "$OUTDIR/$sampid.summ.json.gz" \
      "$OUTDIR/$sampid.muts.json.gz"
  done | parallel -j40 --halt 1
  gunzip $OUTDIR/*.{summ,muts}.json.gz

}

function calc_concordance {
  cd $OUTDIR
  for jsonfn in S*summ.json; do
    sampid=$(echo $jsonfn | cut -d. -f1)
    for foo in $sampid.*.json; do
      echo -n "gzip < $foo > $foo.gz && "
    done
    echo "PYTHONPATH=~/.apps/phylowgs python2 ~/work/steph/protocols/pairwise/calc_concordance.py" \
      "$sampid.{summ,muts}.json.gz" \
      "$sampid.concord.tsv" \
      "$sampid.concord.html"
  done | parallel -j40 --halt 1
}

function write_concord_index {
  cd $OUTDIR
  for htmlfn in S*.concord.html; do
    sampid=$(basename $htmlfn | cut -d. -f1)
    echo "<a href=$htmlfn>$sampid</a><br>"
  done > index.html
}

function add_to_witness {
  witnessdir=$PWGSDIR/witness/data/steph.$RUNNAME.$(date '+%Y%m%d')
  mkdir -p $witnessdir
  cp -a $OUTDIR/*.{summ,muts}.json $witnessdir
  cd $PWGSDIR/witness
  python2 index_data.py
}

function match_clusters {
  cd $BASEDIR/data/pairwise.xeno.nocns
  for ssmfn in $SSMDIR/*.sampled.ssm; do
    sampid=$(basename $ssmfn | cut -d . -f1)
    if [[ $sampid == SJBALL022610 || $sampid == SJETV010 ]]; then
      continue
    fi
    echo "python3 $PROTDIR/match_clusters.py "\
      "$sampid" \
      "$HANDBUILTDIR/$sampid.json" \
      "$ssmfn" \
      "$SSMDIR/$sampid.params.json" \
      "handbuilt.xeno" \
      "handbuilt.patient" \
      "> $sampid.matched.txt"
  done | parallel -j40 --halt 1
}

function main {
  mkdir -p $OUTDIR

  #rename_samples
  #remove_samples

  #calc_pairwise
  plot
  add_tree_indices
  write_plot_index
  add_to_witness

  #calc_concordance
  #write_concord_index

  match_clusters
}

main
