function is_run_complete {
  if ([[ -f $1 ]] && unzip -l $1 | grep "adjm.npy" > /dev/null); then
    return 0
  else
    return 1
  fi
}

function is_big_K {
  runid=$1
  if [[ $runid =~ K30 || $runid =~ K100 ]]; then
    return 0
  else
    return 1
  fi
}

