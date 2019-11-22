#!/bin/bash
command -v java > /dev/null || module load java

BASEDIR=~/work/pairtree
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
JAVA=$HOME/.apps/jre8/bin/java
JOB_DIR=$HOME/jobs
LICHEE_DIR=$HOME/.apps/lichee/LICHeE/release

NUM_TREES=3000

BATCH=sims.smallalpha.lichee
PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/sims.smallalpha.pairtree
#BATCH=steph.xeno.lichee
#PAIRTREE_INPUTS_DIR=$BASEDIR/scratch/inputs/steph.xeno.withgarb.pairtree

INDIR=$BASEDIR/scratch/inputs/$BATCH
OUTBASE=$BASEDIR/scratch/results/$BATCH

export OMP_NUM_THREADS=1

function convert_inputs {
  mkdir -p $INDIR

  for ssmfn in $PAIRTREE_INPUTS_DIR/*.ssm; do
    sampid=$(basename $ssmfn | cut -d. -f1)
    echo "python3 $SCRIPTDIR/convert_inputs.py " \
      "$PAIRTREE_INPUTS_DIR/$sampid.ssm" \
      "$PAIRTREE_INPUTS_DIR/$sampid.params.json" \
      "$INDIR/$sampid.snv" \
      "$INDIR/$sampid.cluster"
  done
}

function run_lichee {
  mkdir -p $OUTBASE
  cd $OUTBASE

  for snvfn in $INDIR/*.snv; do
    runid=$(basename $snvfn | cut -d. -f1)
    outd="$OUTBASE/$runid"
    outfn="$outd/$runid.trees"
    [[ -f $outfn ]] && continue

    cmd=""
    cmd+="#!/bin/bash\n"
    cmd+="#SBATCH --nodes=1\n"
    cmd+="#SBATCH --ntasks=1\n"
    cmd+="#SBATCH --mem=20GB\n"
    cmd+="#SBATCH --time=23:59:00\n"
    cmd+="#SBATCH --job-name lichee_$runid\n"
    cmd+="#SBATCH --output=$JOB_DIR/slurm_lichee_${runid}_%j.txt\n"
    cmd+="#SBATCH --partition=cpu\n"
    cmd+="#SBATCH --mail-type=NONE\n"

    cmd+="mkdir -p $outd && cd $outd && "
    cmd+="TIMEFORMAT='%R %U %S'; time ($JAVA -jar $LICHEE_DIR/lichee.jar "
    cmd+="-build "
    cmd+="-i $INDIR/$runid.snv "
    cmd+="-o $outfn "
    cmd+="-clustersFile $INDIR/$runid.cluster "
    cmd+="-maxVAFAbsent 0.005 "
    cmd+="-minVAFPresent 0.005 "
    cmd+="-n 0 "
    cmd+="-s $NUM_TREES "
    cmd+=">  $outd/$runid.stdout "
    cmd+="2> $outd/$runid.stderr) 2>$outd/$runid.time"

    jobfn=$(mktemp)
    echo -e $cmd > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function convert_outputs {
  cd $OUTBASE
  for treesfn in $OUTBASE/*/*.trees; do
    runid=$(basename $treesfn | cut -d. -f1)
    outd=$(dirname $treesfn)

    cmd="python3 $SCRIPTDIR/convert_outputs.py "
    cmd+="--mutrel $outd/$runid.mutrel.npz "
    cmd+="--structures $outd/$runid.params.json "
    cmd+="$INDIR/$runid.snv "
    cmd+="$outd/$runid.trees "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.params.json "

    echo $cmd
  done
}

function compute_phis {
  cd $OUTBASE
  for structfn in $OUTBASE/*/*.params.json; do
    runid=$(basename $structfn | cut -d. -f1)
    outd=$(dirname $structfn)
    resultsfn=$outd/$runid.results.npz

    # `replace_llh.py` assumes that Pairtree will write trees (with computed
    # phis) in the same order as the structures given in the `$structfn` file.
    # We must disable posterior sorting to preserve this property.
    cmd="python3 $BASEDIR/bin/pairtree "
    cmd+="--seed 1 "
    cmd+="--phi-fitter projection "
    cmd+="--parallel 0 "
    cmd+="--params $structfn "
    cmd+="--disable-posterior-sorting "
    cmd+="$PAIRTREE_INPUTS_DIR/$runid.ssm "
    cmd+="$resultsfn "

    cmd+="&& python3 $SCRIPTDIR/replace_llh.py "
    cmd+="$structfn "
    cmd+="$resultsfn "

    echo $cmd
  done
}

function compute_mutphis {
  cd $OUTBASE

  for resultsfn in $OUTBASE/*/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    outd=$(dirname $resultsfn)
    ssmfn=${PAIRTREE_INPUTS_DIR}/${runid}.ssm 
    mutphifn=$outd/$runid.mutphi.npz

    cmd="python3 $BASEDIR/comparison/pairtree/make_mutphis.py "
    cmd+="$resultsfn "
    cmd+="$ssmfn "
    cmd+="$mutphifn "

    cmd+="&& python3 $BASEDIR/comparison/impute_missing_mutphis.py "
    cmd+="$ssmfn "
    cmd+="$PAIRTREE_INPUTS_DIR/${runid}.params.json "
    cmd+="$mutphifn "

    echo $cmd
  done
}

function main {
  #convert_inputs
  run_lichee
  #convert_outputs | sort --random-sort | parallel -j8 --halt 1 --eta
  #compute_phis | parallel -j40 --halt 1 --eta
  #compute_mutphis | parallel -j80 --halt 1 --eta
}

main
