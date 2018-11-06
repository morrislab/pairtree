#!/bin/bash
set -euo pipefail
module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
RESULTSDIR=$BASEDIR/scratch/results

function main {
  cd $RESULTSDIR
  prefix=steph.xeno.nogarb

  for mutrelfn in ${prefix}.truth/*.mutrel.npz; do
    runid=$(basename $mutrelfn | cut -d. -f1)
    mutrels="truth=${prefix}.truth/$runid.mutrel.npz "
    mutrels+="pwgs=${prefix}.pwgs/$runid/$runid.mutrel.npz "
    mutrels+="pairtree_trees=${prefix}.pairtree/$runid.trees.mutrel.npz"

    for M in $(echo $mutrels | tr ' ' '\n' | cut -d= -f2); do
      [[ ! -f $M ]] && continue 2
    done

    echo "cd $RESULTSDIR && python3 $SCRIPTDIR/eval.py $mutrels > $runid.score.txt"
  done | parallel -j40 --halt 0

  (
    echo 'runid,'$(head -n1 $(ls *.score.txt | head -n1))
    for foo in *.score.txt; do
      echo $(echo $foo | cut -d. -f1),$(tail -n+2 $foo) 
    done
  ) | grep -v ',$' | curl -F c=@- https://ptpb.pw
}

main
