#!/bin/sh
set -eu -o pipefail

function make_results {
  cd ~/work/pairtree/scratch/results/sims.smallalpha.pairtree.truestructs.projection
  for dir in . ../sims.smallalpha.pairtree.projection; do
    for foo in sim*; do
      echo python3 ~/work/pairtree/bin/plottree --plot tree ~/work/pairtree-experiments/inputs/sims.smallalpha.pairtree/$foo.{ssm,params.json} $dir/$foo/$foo.{results.npz,toptree.html}
      echo python3 ~/work/pairtree/bin/summposterior ~/work/pairtree-experiments/inputs/sims.smallalpha.pairtree/$foo.{ssm,params.json} $dir/$foo/$foo.{results.npz,summ.html}
    done
  done | parallel -j40 --halt 1 --eta
}

function extract_json {
  pattern=$1
  egrep --only-matching "tree_json = '[^']+" | cut -c14- | jq $pattern | awk '{printf "%.3f", $0}'
}

function make_index {
  cd ~/work/pairtree/scratch/results
  (
    echo '<!doctype html><html lang="en"><head><meta charset="utf-8">'
    echo '<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css"></head>'
    echo '<body><table class="table table-striped"><thead><tr><th>Run ID</th><th>True tree</th><th>Pairtree tree</th><th>True summary</th><th>Paitree summary</th><th>True nlglh</th><th>Pairtree nlglh</th><th>True prob</th><th>Pairtree prob</th><th>JSD</th></tr></thead><tbody>'
    for foo in sims.smallalpha.pairtree.truestructs.projection/sim*; do
      runid=$(basename $foo)
      tt_truth="sims.smallalpha.pairtree.truestructs.projection/$runid/$runid.toptree.html"
      tt_pairtree="sims.smallalpha.pairtree.projection/$runid/$runid.toptree.html"

      echo "<tr>"
      echo "<td>$runid</td>"
      echo "<td><a href=$tt_truth>true top tree</a></td>"
      echo "<td><a href=$tt_pairtree>pairtree top tree</a></td>"
      echo "<td><a href=sims.smallalpha.pairtree.truestructs.projection/$runid/$runid.summ.html>true summ</a></td>"
      echo "<td><a href=sims.smallalpha.pairtree.projection/$runid/$runid.summ.html>pairtree summ</a></td>"
      for pattern in .nlglh .prob; do
        for fn in $tt_truth $tt_pairtree; do
          echo "<td>$(cat $fn | extract_json $pattern)</td>"
        done
      done
      echo "<td>$(cat ~/work/pairtree/scratch/scores/entropy.txt | grep "^$runid" | awk -F, '{print $16}')</td>"
      echo "</tr>"
    done
    echo "</tbody></table></body></html>"
  ) > treecmp.html
}

function main {
  #make_results
  make_index
}

main
