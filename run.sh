RESULTSDIR=~/tmp/pairtree
PROTDIR=~/work/pairtree
DATADIR=~/work/steph/data

function make_trees {
  #for treetype in patient xeno; do
  for treetype in xeno; do
    OUTDIR=$RESULTSDIR/$treetype
    rm -f $OUTDIR/SJ*.{html,json,csv,stdout}
    mkdir -p $OUTDIR
    cd $OUTDIR

    if [[ $treetype == "xeno" ]]; then
      INPUTSDIR=$DATADIR/inputs/steph.xenos.nocns
    else
      INPUTSDIR=$DATADIR/inputs/steph.vanilla
    fi

    for ssmfn in $INPUTSDIR/SJ*.sampled.ssm; do
      runid=$(basename $ssmfn | cut -d. -f1)
      echo "python3 $PROTDIR/pairtree.py" \
        "--tree-samples 1000" \
        "--phi-iterations 10000" \
        "--tree-chains 20" \
        "--parallel 20" \
          "$runid" \
          "$INPUTSDIR/$runid.{sampled.ssm,params.json}" \
          "$DATADIR/handbuilt_trees/$runid.json handbuilt.$treetype" \
          ">$runid.stdout"
    done | grep -e SJETV010steph -e SJBALL022610steph | parallel -j6 --halt 1

    for stdoutfn in *.stdout; do
      runid=$(basename $stdoutfn | cut -d. -f1)
      (
	echo "cluster1,cluster2,distance,norm_distance"
	cat $stdoutfn | grep ^cluster_dist | cut -d, -f2- | sort -nrk4 -t,
      ) > $runid.cluster_distances.csv
      (
	echo "phi,cluster,distance,norm_distance"
	cat $stdoutfn | grep ^phi_dist | cut -d, -f2- | sort -nrk4 -t,
      ) > $runid.phi_distances.csv
    done

    for resultsfn in *steph*.results.html; do
      runid=$(echo $resultsfn | cut -d. -f1)
      echo "<h3>$runid</h3>"
      echo "<ul>"
      echo "<li><a href=$runid.results.html>$runid results</a></li>"
      echo "<li><a href=$runid.cluster_distances.csv>$runid cluster distances</a></li>"
      echo "<li><a href=$runid.phi_distances.csv>$runid phi distances</a></li>"
      echo "</ul>"
    done > index.html
  done
}

function main {
  make_trees
}

main
