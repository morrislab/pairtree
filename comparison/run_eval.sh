#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
RESULTSDIR=$BASEDIR/scratch/results
SCORESDIR=$BASEDIR/scratch/scores
PARALLEL=80

function para {
  parallel -j$PARALLEL --halt 1
}

function make_sims_truth {
  mkdir -p $TRUTH_DIR

  for datafn in $PAIRTREE_INPUTS_DIR/*.truth.pickle; do
    runid=$(basename $datafn | cut -d. -f1)
    truthdir=$TRUTH_DIR/$runid

    cmd="mkdir -p $truthdir/$runid"

    cmd+="&& OMP_NUM_THREADS=1 python3 $SCRIPTDIR/enum_true_trees.py "
    cmd+="--check-delta "
    cmd+="$datafn "
    cmd+="$truthdir/$runid.results.npz "

    cmd+="&& OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutphi.py "
    cmd+="$datafn "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$truthdir/$runid.trees_mutphi.npz"
    echo $cmd
  done | parallel -j80 --halt 1 --eta

  for resultsfn in $TRUTH_DIR/*/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    truthdir=$TRUTH_DIR/$runid

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/pairtree/make_mutrels.py "
    cmd+="--trees-mutrel $truthdir/$runid.trees_mutrel.npz "
    cmd+="$truthdir/$runid.results.npz "
    echo $cmd
  done | parallel -j10 --halt 1 --eta
}

function make_mle_mutphis {
  mkdir -p $MLE_MUTPHIS_DIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    #[[ $runid =~ K30 || $runid =~ K100 ]] || continue

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

  if [[ $BATCH == steph ]]; then
    paths+="truth=${TRUTH_DIR}/$runid/${runid}.trees_${result_type}.npz "
    paths+="pairtree_multi=${BATCH}.pairtree.multichain/${runid}/${runid}.trees_${result_type}.npz "
    #paths+="pairtree_single=${BATCH}.pairtree.singlechain/${runid}/${runid}.pairtree_trees.all.llh.$result_type.npz "
    #paths+="pwgs_allvars=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    #paths+="pplus_allvars=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
    #paths+="pwgs_supervars=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    #paths+="pplus_supervars=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
    #paths+="lichee=${BATCH}.xeno.lichee/$runid/$runid.$result_type.npz "
    #if [[ $result_type == mutrel ]]; then
    #  paths+="pairtree_tensor=${BATCH}.pairtree.onlytensor/${runid}/${runid}.pairtree_clustrel.mutrel.npz "
    #fi
  elif [[ $BATCH == sims.smallalpha ]]; then
    paths+="truth=${TRUTH_DIR}/$runid/${runid}.trees_${result_type}.npz "
    #paths+="pwgs_supervars=sims.pwgs.supervars/$runid/$runid.pwgs_trees_single.$result_type.npz "
    #paths+="pastri=${BATCH}.pastri.informative/$runid/$runid.pastri_trees_llh.$result_type.npz "
    #for foo in $(seq 69 77); do
    #  paths+="pairtree_lol${foo}=sims.pairtree.lol${foo}/${runid}/${runid}.trees_${result_type}.npz "
    #done
    #paths+="pairtree_single=sims.pairtree.singlechain/${runid}/${runid}.pairtree_trees.all.llh.$result_type.npz "
    #if [[ $result_type == mutrel ]]; then
    #  paths+="pairtree_clustrel=sims.pairtree.multichain/${runid}/${runid}.clustrel_${result_type}.npz "
    #fi
    paths+="pairtree_multi=${BATCH}.pairtree.multichain/${runid}/${runid}.trees_${result_type}.npz "
    #paths+="pairtree_quad=sims.pairtree.quadchain/${runid}/${runid}.pairtree_trees.all.llh.$result_type.npz "
    #paths+="pairtree_single_old=sims.pairtree.projection.singlechain.old_proposals/${runid}/${runid}.pairtree_trees.all.llh.$result_type.npz "
    #paths+="pairtree_multi_old=sims.pairtree.projection.multichain.old_proposals/${runid}/${runid}.pairtree_trees.all.llh.$result_type.npz "
    #paths+="lichee=sims.lichee/$runid/$runid.$result_type.npz "
    #paths+="citup=sims.citup.rawvars.qip/$runid/$runid.$result_type.npz "
  fi

  echo $paths
}

function eval_mutrels {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for truthfn in $(ls $TRUTH_DIR/*/*.trees_mutrel.npz | sort --random-sort); do
    runid=$(basename $truthfn | cut -d. -f1)
    mutrels=$(make_results_paths $runid mutrel)

    cmd="cd $RESULTSDIR && "
    cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutrels.py "
    cmd+="--discard-garbage "
    if [[ $BATCH == steph ]]; then
      for method in {pwgs,pplus}_allvars; do
        cmd+="--ignore-garbage-for $method "
      done
    fi
    cmd+="--params $PAIRTREE_INPUTS_DIR/$runid.params.json "
    cmd+="$mutrels "
    cmd+="> $SCORESDIR/$BATCH/$runid.mutrel_score.txt"
    echo $cmd
  done
}

function eval_mutphis {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for mutphifn in $(ls $TRUTH_DIR/*/*.trees_mutphi.npz | sort --random-sort); do
    runid=$(basename $mutphifn | cut -d. -f1)
    mutphis=$(make_results_paths $runid mutphi)

    echo "cd $RESULTSDIR && " \
      "OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutphis.py " \
      "--params $PAIRTREE_INPUTS_DIR/$runid.params.json" \
      "$mutphis " \
      "> $SCORESDIR/$BATCH/$runid.mutphi_score.txt"
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

      # These are the runs we did not use in the paper.
      for bad_sampid in SJETV010{,nohypermut,stephR1,stephR2} SJBALL022610; do
        if [[ $S == $bad_sampid ]]; then
          continue 2
        fi
      done
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
      #cmd+="pairtree_tensor=$RESULTSDIR/sims.pairtree.onlytensor/$runid/$runid.time "
      #cmd+="pairtree_single=$RESULTSDIR/sims.pairtree.singlechain/$runid/$runid.time "
      cmd+="pairtree_multi=$RESULTSDIR/sims.pairtree.multichain/$runid/$runid.time "
      #cmd+="pairtree_quad=$RESULTSDIR/sims.pairtree.quadchain/$runid/$runid.time "
      #cmd+="pairtree_single_old=$RESULTSDIR/sims.pairtree.projection.singlechain.old_proposals/$runid/$runid.time "
      #cmd+="pairtree_multi_old=$RESULTSDIR/sims.pairtree.projection.multichain.old_proposals/$runid/$runid.time "
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
  threshold=10

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

function plot_individual {
  ptype=$1
  plot_type=$2
  infn=$3
  outfn=$4

  cd $SCORESDIR

  if [[ $(cat $infn | wc -l) -le 1 ]]; then
    # Since we don't generate mutrel for bigK, it will be an empty file
    # with just the header.
    continue
  fi

  cmd="python3 $SCRIPTDIR/plot_individual.py "
  cmd+="--template plotly_white "
  cmd+="--score-type $ptype "
  cmd+="--plot-type $plot_type "
  if [[ $BATCH == "sims" ]]; then
    cmd+="--S-threshold 10 "
  fi
  if [[ $ptype == "mutphi" || $ptype == "mutrel" ]]; then
    cmd+="--hide-method mle_unconstrained "
  fi
  if [[ $ptype == "mutrel" ]]; then
    cmd+="--bandwidth 1 "
  elif [[ $ptype == "mutphi" ]]; then
    cmd+="--bandwidth 0.13 "
  elif [[ $ptype =~ time$ ]]; then
    cmd+="--bandwidth 0.18 "
  fi
  if [[ $ptype == "mutphi" && $BATCH == "sims" ]]; then
    cmd+="--baseline truth "
  fi
  if [[ $ptype == "mutphi" && $BATCH == "steph" ]]; then
    cmd+="--baseline truth "
  fi
  if [[ $ptype =~ time$ ]]; then
    cmd+="--log-y-axis "
  fi
  cmd+="$infn $outfn "

  echo $cmd
}

function plot_results_sims {
  eval_mutrels | parallel -j2 --halt 1 --eta
  eval_mutphis | para

  for type in mutrel mutphi; do
    compile_scores $type
    partition_by_K $type
  done

  (
    basefn="$SCORESDIR/$BATCH"
    for ksize in smallK bigK; do
      #plot_individual mutrel box $basefn.mutrel.$ksize.{txt,html}
      plot_individual mutphi box $basefn.mutphi.$ksize.{txt,html}
    done
  ) | para
}

function plot_results_steph {
  #make_mle_mutphis

  (eval_mutrels; eval_mutphis) | para
  compile_scores mutrel
  compile_scores mutphi
  (
    basefn="$SCORESDIR/$BATCH"
    plot_individual mutrel box $basefn.mutrel.{txt,html}
    plot_individual mutphi box $basefn.mutphi.{txt,html}
  ) | para
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
  #eval_runtime | para
  #compile_runtime
  #remove_missing
  #partition_by_K cputime
  #partition_by_K walltime
  for ksize in bigK smallK; do
    for ttype in cputime walltime; do
      basefn="$SCORESDIR/$BATCH.$ttype.$ksize"
      plot_individual $ttype violin $basefn.{txt,html}
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

  #export BATCH=steph
  #export PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree
  #export TRUTH_DIR=$RESULTSDIR/steph.pairtree.hbstruct
  #export MLE_MUTPHIS_DIR=$RESULTSDIR/${BATCH}.mle_unconstrained
  #plot_results_steph
}

main
