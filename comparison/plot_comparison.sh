#!/bin/bash
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
BASEDIR=~/work/pairtree
SCORESDIR=$BASEDIR/scratch/scores
BATCH=sims.smallalpha

function plot_comparison {
  # To redirect port 80 to 8000:
  #   iptables -A PREROUTING -t nat -i eth0 -p tcp --dport 80 -j REDIRECT --to-port 8000
  cd $SCRIPTDIR

  export MUTREL_RESULTS=$SCORESDIR/$BATCH.mutrel.txt
  export MUTPHI_RESULTS=$SCORESDIR/$BATCH.mutphi.txt
  export SIM_PARAMS="GKMST"
  production=false

  if [[ $production == true ]]; then
    gunicorn -w 4 -b 0.0.0.0:8000 plot_comparison:server
  else
    python3 plot_comparison.py
  fi
}

plot_comparison
