#!/bin/bash
command -v parallel > /dev/null || module load gnu-parallel

SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
RESULTSDIR=$BASEDIR/scratch/results
SCORESDIR=$BASEDIR/scratch/scores
PARALLEL=40

function para {
  parallel -j$PARALLEL --halt 1
}

function make_sims_truth {
  mkdir -p $TRUTH_DIR

  for datafn in $PAIRTREE_INPUTS_DIR/*.data.pickle; do
    runid=$(basename $datafn | cut -d. -f1)

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutrel.py "
    if [[ $BATCH == "sims" && !($runid =~ K30 || $runid =~ K100) ]]; then
      cmd+="--enumerate-trees "
    fi
    cmd+="$datafn "
    cmd+="$TRUTH_DIR/$runid.mutrel.npz"
    echo $cmd

    cmd="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/make_truth_mutphi.py "
    cmd+="$datafn "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$TRUTH_DIR/$runid.mutphi.npz"
    echo $cmd
  done
}

function make_mle_mutphis {
  mkdir -p $MLE_MUTPHIS_DIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    [[ $runid =~ K30 || $runid =~ K100 ]] || continue

    echo "python3 $SCRIPTDIR/make_mle_mutphis.py" \
      "$ssmfn" \
      "$MLE_MUTPHIS_DIR/$runid.mutphi.npz"
  done
}

function make_results_paths {
  runid=$1
  result_type=$2
  paths=""

  if [[ $result_type == mutphi ]]; then
    paths+="mle_unconstrained=${BATCH}.mle_unconstrained/$runid.$result_type.npz "
  fi

  if [[ $BATCH == steph ]]; then
    truth_path="${BATCH}.pairtree.hbstruct/${runid}.pairtree_trees.all.llh.${result_type}.npz"
    paths+="truth=$truth_path "
    paths+="pairtree_trees_llh=${BATCH}.xeno.withgarb.pairtree/${runid}.pairtree_trees.all.llh.$result_type.npz "
    paths+="pairtree_trees_uniform=${BATCH}.xeno.withgarb.pairtree/${runid}.pairtree_trees.all.uniform.$result_type.npz "
    paths+="singlepairtree_trees_llh=${BATCH}.xeno.withgarb.pairtree/${runid}.pairtree_trees.subset.llh.$result_type.npz "
    paths+="singlepairtree_trees_uniform=${BATCH}.xeno.withgarb.pairtree/${runid}.pairtree_trees.subset.uniform.$result_type.npz "
    paths+="pairtree_handbuilt=$truth_path "
    paths+="pwgs_allvars_single_uniform=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_single_uniform.$result_type.npz "
    paths+="pwgs_allvars_single_llh=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    paths+="pwgs_allvars_multi_uniform=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_multi_uniform.$result_type.npz "
    paths+="pwgs_allvars_multi_llh=${BATCH}.pwgs.allvars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
    paths+="pwgs_supervars_single_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_uniform.$result_type.npz "
    paths+="pwgs_supervars_single_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    paths+="pwgs_supervars_multi_uniform=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_uniform.$result_type.npz "
    paths+="pwgs_supervars_multi_llh=${BATCH}.pwgs.supervars/$runid/$runid.pwgs_trees_multi_llh.$result_type.npz "
  elif [[ $BATCH == sims ]]; then
    paths+="truth=${TRUTH_DIR}/$runid.$result_type.npz "
    paths+="pwgs_supervars_single_llh=sims.pwgs.supervars/$runid/$runid.pwgs_trees_single_llh.$result_type.npz "
    paths+="pastri_llh=${BATCH}.pastri.informative/$runid/$runid.pastri_trees_llh.$result_type.npz "
    paths+="pairtree_proj_single_llh=sims.pairtree.projection.singlechain/${runid}.pairtree_trees.all.llh.$result_type.npz "
    paths+="pairtree_proj_multi_llh=sims.pairtree.projection.multichain/${runid}.pairtree_trees.all.llh.$result_type.npz "
    paths+="pairtree_clustrel=sims.pairtree.onlytensor/$runid.pairtree_clustrel.$result_type.npz "
    paths+="lichee_llh=sims.lichee/$runid.llh.$result_type.npz "
  fi

  echo $paths
}

function eval_mutrels {
  cd $RESULTSDIR
  mkdir -p $SCORESDIR/$BATCH

  for truthfn in $(ls $TRUTH_DIR/*.mutrel.npz | sort --random-sort); do
    runid=$(basename $truthfn | cut -d. -f1)
    mutrels=$(make_results_paths $runid mutrel)

    cmd="cd $RESULTSDIR && "
    cmd+="OMP_NUM_THREADS=1 python3 $SCRIPTDIR/eval_mutrels.py "
    cmd+="--discard-garbage "
    if [[ $BATCH == steph ]]; then
      for method in pwgs_allvars_{single,multi}_{uniform,llh}; do
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

  for mutphifn in $(ls $TRUTH_DIR/*.mutphi.npz | sort --random-sort); do
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
        echo "Methods in $foo don't match expected $methods" >&2
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
  for datafn in $PAIRTREE_INPUTS_DIR/*.data.pickle; do
    runid=$(basename $datafn | cut -d. -f1)

    for timetype in cpu wall; do
      cmd="python3 $SCRIPTDIR/eval_runtime.py "
      cmd+="--time-type $timetype "
      cmd+="lichee=$RESULTSDIR/sims.lichee/$runid.time "
      cmd+="pairtree_tensor=$RESULTSDIR/sims.pairtree.onlytensor/$runid.time "
      cmd+="pairtree_single=$RESULTSDIR/sims.pairtree.projection.singlechain/$runid.time "
      cmd+="pairtree_multi=$RESULTSDIR/sims.pairtree.projection.multichain/$runid.time "
      cmd+="pastri=$RESULTSDIR/sims.pastri.informative.complete/$runid/$runid.time "
      cmd+="pwgs=$RESULTSDIR/sims.pwgs.supervars/$runid/$runid.time "
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

  cd $SCORESDIR
  for ktype in bigK smallK; do
    basefn="$SCORESDIR/$BATCH.$ptype.$ktype"
    infn=$basefn.txt
    outfn=$basefn.html

    if [[ $(cat $infn | wc -l) -le 1 ]]; then
      # Since we don't generate mutrel for bigK, it will be an empty file
      # with just the header.
      continue
    fi

    cmd="python3 $SCRIPTDIR/plot_individual.py "
    cmd+="--template plotly_white "
    cmd+="--score-type $ptype "
    cmd+="--plot-type $plot_type "
    cmd+="--S-threshold 10 "
    cmd+="--bandwidth 0.07 "
    if [[ $ptype == "mutphi" && $BATCH == "sims" ]]; then
      cmd+="--baseline truth "
    fi
    if [[ $ptype == "mutphi" && $BATCH == "steph" ]]; then
      cmd+="--baseline pairtree_handbuilt "
    fi
    if [[ $ptype =~ time$ ]]; then
      cmd+="--log-y-axis "
    fi
    cmd+="$( [[ $BATCH == sims ]] && echo --partition-by-samples) "
    cmd+="$infn $outfn "
    echo $cmd
  done
}

function run_batch {
  export MLE_MUTPHIS_DIR=$RESULTSDIR/${BATCH}.mle_unconstrained
  #[[ $BATCH == "sims" ]] && make_sims_truth
  #make_mle_mutphis

  #(eval_mutrels; eval_mutphis) | para 2>$SCRATCH/tmp/eval.log
  #compile_scores mutrel
  #compile_scores mutphi
  #partition_by_K mutrel
  #partition_by_K mutphi
  #(plot_individual mutrel violin; plot_individual mutphi box) | para

  #eval_runtime | para
  #compile_runtime
  #partition_by_K cputime
  #partition_by_K walltime
  (plot_individual cputime violin; plot_individual walltime violin) | para
}

function main {
  export BATCH=sims
  export PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.pairtree
  export TRUTH_DIR=$RESULTSDIR/sims.truth
  run_batch

  #export BATCH=steph
  #export PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree
  #export TRUTH_DIR=$RESULTSDIR/steph.pairtree.hbstruct
  #run_batch
}

main
