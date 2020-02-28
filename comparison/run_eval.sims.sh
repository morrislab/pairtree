#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
RESULTSDIR=$BASEDIR/scratch/results
SCORESDIR=$BASEDIR/scratch/scores
PARALLEL=80

function para {
  parallel -j$PARALLEL --halt 1 --eta
}

function make_sims_truth {
  mkdir -p $TRUTH_DIR

  for truthfn in $PAIRTREE_INPUTS_DIR/*.truth.pickle; do
    runid=$(basename $truthfn | cut -d. -f1)
    truthdir=$TRUTH_DIR/$runid
    mkdir -p $truthdir/$runid

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/enum_true_trees.py "
    cmd+="--check-delta "
    cmd+="$truthfn "
    cmd+="$truthdir/$runid.results.npz "
    echo $cmd

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutphi.py "
    cmd+="$truthfn "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$truthdir/$runid.mutphi.npz"
    echo $cmd

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_baseline_mutdist_from_sims.py "
    cmd+="$truthfn "
    cmd+="$truthdir/$runid.phi.npz"
    echo $cmd
  done | parallel -j80 --halt 1 --eta

  for resultsfn in $TRUTH_DIR/*/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    truthdir=$TRUTH_DIR/$runid

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/pairtree/convert_outputs.py "
    cmd+="$truthdir/$runid.results.npz "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "
    cmd+="$truthdir/$runid.neutree.npz "
    cmd+="&& OMP_NUM_THREADS=1 python3 $SCRIPTDIR/neutree/make_mutrels.py "
    cmd+="$truthdir/$runid.neutree.npz "
    cmd+="$truthdir/$runid.mutrel.npz "
    echo $cmd
  done | parallel -j10 --halt 1 --eta
}

function make_mle_mutphis {
  mkdir -p $MLE_MUTPHIS_DIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)

    echo "python3 $SCRIPTDIR/make_mle_mutphis.py" \
      "$ssmfn" \
      "$MLE_MUTPHIS_DIR/$runid.mutphi.npz"
  done | parallel -j80 --halt 1 --eta
}

function make_results_paths {
  runid=$1
  result_type=$2
  paths=""

  if [[ $result_type == mutphi ]]; then
    paths+="mle_unconstrained=${BATCH}.mle_unconstrained/$runid.$result_type.npz "
  fi

  if [[ $result_type == mutphi || $result_type == mutrel ]]; then
    paths+="truth=${TRUTH_DIR}/$runid/${runid}.${result_type}.npz "
  fi
  paths+="pwgs_supervars=${BATCH}.pwgs.supervars/$runid/$runid.$result_type.npz "
  paths+="pastri=${BATCH}.pastri/$runid/$runid.$result_type.npz "
  if [[ $result_type == mutrel ]]; then
    paths+="pairtree_clustrel=${BATCH}.pairtree/${runid}/${runid}.clustrel_mutrel.npz "
  fi
  paths+="pairtree=${BATCH}.pairtree/${runid}/${runid}.${result_type}.npz "
  paths+="lichee=${BATCH}.lichee/$runid/$runid.$result_type.npz "
  paths+="citup=${BATCH}.citup.rawvars.qip/$runid/$runid.$result_type.npz "

  echo $paths
}

function eval_mutrels {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for truthfn in $(ls $TRUTH_DIR/*/*.mutrel.npz | sort --random-sort); do
    runid=$(basename $truthfn | cut -d. -f1)
    mutrels=$(make_results_paths $runid mutrel)

    cmd="cd $RESULTSDIR && "
    cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutrels.py "
    cmd+="--discard-garbage "
    cmd+="--params $PAIRTREE_INPUTS_DIR/$runid.params.json "
    cmd+="$mutrels "
    cmd+="> $SCORESDIR/$BATCH/$runid.mutrel_score.txt"
    echo $cmd
  done
}

function eval_mutphis {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for mutphifn in $(ls $TRUTH_DIR/*/*.mutphi.npz | sort --random-sort); do
    runid=$(basename $mutphifn | cut -d. -f1)
    mutphis=$(make_results_paths $runid mutphi)

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutphis.py " \
      "--params $PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$mutphis " \
      "> $SCORESDIR/$BATCH/$runid.mutphi_score.txt"
  done
}

function eval_mutdists {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for phifn in $(ls $TRUTH_DIR/*/*.phi.npz | sort --random-sort); do
    # We don't actually use `phifn`, but that's okay -- we just want a
    # reference for what runs exist.
    runid=$(basename $phifn | cut -d. -f1)
    mutdists=$(make_results_paths $runid mutdist)

    for P in 1 2; do
      cmd="cd $RESULTSDIR && "
      cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutdists.py "
      cmd+="-p $P "
      cmd+="--params $PAIRTREE_INPUTS_DIR/$runid.params.json "
      cmd+="$mutdists "
      cmd+="> $SCORESDIR/$BATCH/$runid.mutdistl${P}_score.txt "
      echo $cmd
    done
  done
}

function compile_scores {
  score_type=$1
  suffix=${score_type}_score.txt
  outfn=$SCORESDIR/$BATCH.$score_type.txt
  (
    cd $SCORESDIR/$BATCH
    methods=$(head -n1 $(ls *.$suffix | head -n1))
    echo 'runid,'$methods
    for foo in *.$suffix; do
      if [[ $(head -n1 $foo) != $methods ]]; then
        echo "Methods in $foo ($(head -n1 $foo)) don't match expected $methods" >&2
        exit 1
      fi
      S=$(echo $foo | cut -d. -f1)

      echo $S,$(tail -n+2 $foo)
    done
  ) > $outfn
  #cat $outfn | curl -F c=@- https://ptpb.pw >&2
}

function eval_runtime {
  for datafn in $PAIRTREE_INPUTS_DIR/*.truth.pickle; do
    runid=$(basename $datafn | cut -d. -f1)

    for timetype in cpu wall; do
      cmd="python3 $SCRIPTDIR/eval_runtime.py "
      cmd+="--time-type $timetype "
      cmd+="citup=$RESULTSDIR/sims.citup.rawvars.qip/$runid/$runid.time "
      cmd+="lichee=$RESULTSDIR/sims.lichee/$runid/$runid.time "
      cmd+="pairtree=$RESULTSDIR/sims.pairtree/$runid/$runid.time "
      cmd+="pastri=$RESULTSDIR/sims.pastri.informative/$runid/$runid.time "
      cmd+="pwgs_supervars=$RESULTSDIR/sims.pwgs.supervars/$runid/$runid.time "
      cmd+="> $SCORESDIR/$BATCH/$runid.${timetype}time.txt "
      echo $cmd
    done
  done
}

function compile_runtime {
  for timetype in cpu wall; do
    outfn=$SCORESDIR/$BATCH.${timetype}time.txt
    methods=$(head -n1 $(ls $SCORESDIR/$BATCH/*.${timetype}time.txt | head -n1))
    (
      echo 'runid,'$methods
      for foo in $SCORESDIR/$BATCH/*.${timetype}time.txt; do
        if [[ $(head -n1 $foo) != $methods ]]; then
          echo "Methods in $foo don't match expected $methods" >&2
          exit 1
        fi
        runid=$(basename $foo | cut -d. -f1)
        echo $runid,$(tail -n+2 $foo)
      done
    ) > $outfn
  done
}

function partition_by_K {
  ptype=$1
  threshold=30

  srcfn=$SCORESDIR/$BATCH.$ptype.txt
  dstfn_bigK=$SCORESDIR/$BATCH.$ptype.bigK.txt
  dstfn_smallK=$SCORESDIR/$BATCH.$ptype.smallK.txt

  header=$(head -n1 $srcfn)
  for dstfn in $dstfn_bigK $dstfn_smallK; do
    echo $header > $dstfn
  done

  tail -n+2 $srcfn | while read line; do
    runid=$(echo $line | cut -d, -f1)
    K=$(echo $runid | egrep --only-matching 'K[[:digit:]]+' | cut -c2-)
    if [[ $K -le $threshold ]]; then
      dstfn=$dstfn_smallK
    else
      dstfn=$dstfn_bigK
    fi
    echo $line >> $dstfn
  done
}

function plot_single_method {
  method=$1
  score_type=$2
  infn=$3
  outfn=$4

  cmd="python3 $SCRIPTDIR/plotter/plot_single_method.py "
  cmd+="--template plotly_white "
  cmd+="--score-type $score_type "
  if [[ $score_type == mutphi ]]; then
    cmd+="--baseline truth "
  fi
  cmd+="$infn $method $outfn"
  echo $cmd
}

function make_index {
  cd $SCORESDIR
  (
    for anal in mutdistl1 mutdistl2 mutrel mutphi; do
      echo "<h3>$anal</h3><ul>"
      for foo in $BATCH.*$anal*.html; do
        echo "<li><a href=$foo>$foo</a></li>"
      done
      echo "</ul>"
    done
  ) > index.${BATCH}.html
}

function plot_sims {
  score_type=$1
  infn=$2
  outfn=$3

  cmd="python3 $SCRIPTDIR/plotter/plot_sims.py "
  cmd+="--template seaborn "
  cmd+="--score-type $score_type "
  if [[ $score_type == mutphi ]]; then
    cmd+="--baseline truth "
  fi
  cmd+="$infn $outfn"
  echo $cmd
}


function plot_results_sims {
  #(
  #  eval_mutphis
  #  eval_mutdists
  #) | para
  #eval_mutrels | parallel -j2 --halt 1 --eta

  #for type in mutrel mutphi mutdistl1 mutdistl2; do
  #  compile_scores $type
  #done

  (
    basefn="$SCORESDIR/$BATCH"

    for score_type in mutphi mutrel mutdistl1 mutdistl2; do
      plot_sims  $score_type ${basefn}.${score_type}.txt ${basefn}.lol.${score_type}.html

      for method in pairtree pairtree_clustrel lichee pwgs_supervars citup pastri; do
        [[ $method == pairtree_clustrel && $score_type != mutrel ]] && continue
        plot_single_method $method $score_type ${basefn}.${score_type}.txt ${basefn}.${method}.${score_type}.html
      done
    done
  ) | para

  make_index
}

function remove_missing {
  basefn=$SCORESDIR/$BATCH
  for ttype in cputime walltime; do
    cmd="python3 $SCRIPTDIR/intersect_scores.py "
    cmd+="$basefn.$ttype.txt "
    cmd+="$basefn.mutphi.txt "
    echo $cmd
  done | bash
}

function plot_runtime {
  eval_runtime | para
  compile_runtime
  remove_missing
  partition_by_K cputime
  partition_by_K walltime
  for ksize in bigK smallK; do
    for ttype in cputime walltime; do
      basefn="$SCORESDIR/$BATCH.$ttype.$ksize"
      # I removed the `plot_all_methods` function because I hate this code.
      # I'll need to change this.
      plot_all_methods $ttype violin $basefn.{txt,html}
    done
  done | para
}


function main {
  export BATCH=sims.smallalpha
  export PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/${BATCH}.pairtree
  export TRUTH_DIR=$RESULTSDIR/${BATCH}.truth
  export MLE_MUTPHIS_DIR=$RESULTSDIR/${BATCH}.mle_unconstrained
  #make_sims_truth
  #make_mle_mutphis
  plot_results_sims
  #plot_runtime
}

main
