#!/bin/bash

# Date: 03/21/2022
# Author: Ethan Kulman
# This script runs the pipeline of Pairtree scripts described 
# in the Pairtree STAR Methods protocol paper.

# This script expects itself to be in the /pairtree/star_methods/ directory
# along with the following files: methods.ssm, methods.params.json 
 
# At each step, we create a folder to place temporary files in. If the 
# argument 'clean' is passed in (see usage), then the script will automatically 
# remove the temporary folders and files and leave only the results folder.

# Usage: $ ./run_pipeline.sh /path/to/pairtree
# Usage: $ ./run_pipeline.sh /path/to/pairtree clean

#######################################################################
# -------------------- Constants: directories ----------------------- #
#######################################################################

# folder name this script lives in 
FOLDER_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# They are named after the scripts which produced the files they contain
DIR_FBVRP=$FOLDER_PATH/"fix_bad_var_read_prob"
DIR_RMG=$FOLDER_PATH/"removegarbage"
DIR_CV=$FOLDER_PATH/"clustervars"
DIR_PAIRTREE=$FOLDER_PATH/"pairtree"
DIR_PLOTTREE=$FOLDER_PATH/"plottree"
DIR_SUMMPOSTERIOR=$FOLDER_PATH/"summposterior"

# directories that will be created on startup if they don't already exist
DIRS_TO_CREATE=($DIR_FBVRP $DIR_RMG $DIR_CV \
                $DIR_PAIRTREE $DIR_PLOTTREE $DIR_SUMMPOSTERIOR)

# directories that will be removed if the 'clean' flag is passed
DIRS_TO_CLEAN=($DIR_FBVRP $DIR_RMG $DIR_CV)


#######################################################################
# --------------------- Constants: file names ----------------------- #
#######################################################################

# starting files for method
FN_SSM="methods.ssm"
FN_PARAMS="methods.params.json"

# temp files for fix_bad_var_read_prob
FN_SSM_FBVRP="methods.fbvrp.ssm"
FN_PARAMS_FBVRP="methods.fbvrp.params.json"

# temp files for removegarbage
FN_SSM_RMG="methods.rmg.ssm"
FN_PARAMS_RMG="methods.rmg.params.json"

# temp files for clustervars 
FN_SSM_CV="methods.cv.ssm"
FN_PARAMS_CV="methods.cv.params.json"

# files for pairtree 
FN_SSM_PAIRTREE="methods.results.ssm"
FN_PARAMS_PAIRTREE="methods.results.params.json"
FN_NPZ_PAIRTREE="methods.results.npz"

# files for plottree 
FN_HTML_PLOTTREE="methods.plottree.html"

# files for summposterior 
FN_HTML_SUMMPOSTERIOR="methods.summposterior.html"
                      

#######################################################################
# ------------- Functions: create/remove directories ---------------- #
#######################################################################

# function for creating directories
create_dirs () {
      
  for dir in "${DIRS_TO_CREATE[@]}"; do
    [ -d $dir ] && rm -r $dir
    mkdir $dir
  done
  
} 

# function for cleaning directories 
clean_dirs () {
  
  for dir in "${DIRS_TO_CLEAN[@]}"; do
    [ -d $dir ] && rm -r $dir
  done
  
} 

#######################################################################
# ----------- Start of pipeline for STAR Methods protocol ----------- #
#######################################################################

# ------------- Print out welcome message =) -------------- #
echo ""
echo "---------------------------------------------------------------"
echo "Running pairtree pipeline script for the STAR Methods protocol!"
echo "---------------------------------------------------------------"

# ------------- Check that pairtree directory is valid -------------- #
PTDIR=$PWD
([ ! -d  $PTDIR ] || [ ! -d $PTDIR/bin ] || [ ! -d $PTDIR/util ]) && \
  echo "Script must be run from /path/to/pairtree directory. Usage: run_pipeline.sh [optional clean]" && exit 1


# ----------- Check that files for running protocol exist ----------- #
([ ! -f $FOLDER_PATH/$FN_SSM ] || [ ! -f $FOLDER_PATH/$FN_PARAMS ]) && \
  echo "methods.ssm or methods.params.json is not in the same directory as run_pipeline.sh" && exit 1


# ----------------------- Create directories ------------------------ #
create_dirs && echo "Creating directories..."
 

# -------------------- (1) Run fix_bad_var_read_prob -------------------- #
echo && echo "Running utils/fix_bad_var_read_prob.py"
python3 $PTDIR/util/fix_bad_var_read_prob.py $FOLDER_PATH/$FN_SSM $FOLDER_PATH/$FN_PARAMS \
                                             $DIR_FBVRP/$FN_SSM_FBVRP $DIR_FBVRP/$FN_PARAMS_FBVRP \
                                             --action modify_var_read_prob


# -------------------- (2) Run removegarbage -------------------- #
echo && echo "Running bin/removegarbage"
python3 $PTDIR/bin/removegarbage $DIR_FBVRP/$FN_SSM_FBVRP $DIR_FBVRP/$FN_PARAMS_FBVRP \
                                 $DIR_RMG/$FN_PARAMS_RMG
                         
# since we don't modify the .ssm, just copy it to the removegarbage directory
cp $DIR_FBVRP/$FN_SSM_FBVRP $DIR_RMG/$FN_SSM_RMG


# -------------------- (3) Run clustervars ---------------------- #
echo && echo "Running bin/clustervars"
python3 $PTDIR/bin/clustervars $DIR_RMG/$FN_SSM_RMG $DIR_RMG/$FN_PARAMS_RMG \
                               $DIR_CV/$FN_PARAMS_CV
                         
# since we don't modify the .ssm, just copy it to the clustervars directory
cp $DIR_RMG/$FN_SSM_RMG $DIR_CV/$FN_SSM_CV


# -------------------- (4) Run pairtree -------------------- #
echo && echo "Running bin/pairtree"
python3 $PTDIR/bin/pairtree --params $DIR_CV/$FN_PARAMS_CV --seed 5555 \
                            $DIR_CV/$FN_SSM_CV $DIR_PAIRTREE/$FN_NPZ_PAIRTREE 
                         
# since we don't modify the .ssm, just copy it to the pairtree directory
cp $DIR_CV/$FN_SSM_CV $DIR_PAIRTREE/$FN_SSM_PAIRTREE

# since we don't modify the .params.json, just copy it to the pairtree directory
cp $DIR_CV/$FN_PARAMS_CV $DIR_PAIRTREE/$FN_PARAMS_PAIRTREE


# -------------------- (5) Run plottree -------------------- #
echo && echo "Running bin/plottree"
python3 $PTDIR/bin/plottree --runid star_methods \
                            $DIR_PAIRTREE/$FN_SSM_PAIRTREE $DIR_PAIRTREE/$FN_PARAMS_PAIRTREE \
                            $DIR_PAIRTREE/$FN_NPZ_PAIRTREE $DIR_PLOTTREE/$FN_HTML_PLOTTREE
                      
                      
# -------------------- (6) Run summposterior --------------- #
echo && echo "Running bin/summposterior"
python3 $PTDIR/bin/summposterior --runid star_methods \
                                $DIR_PAIRTREE/$FN_SSM_PAIRTREE $DIR_PAIRTREE/$FN_PARAMS_PAIRTREE \
                                $DIR_PAIRTREE/$FN_NPZ_PAIRTREE $DIR_SUMMPOSTERIOR/$FN_HTML_SUMMPOSTERIOR
                      

# ------------------------ Clean directories ------------------------ #
[[ $1 == "clean" ]] && clean_dirs && echo "Cleaning directories..."


# ----------------------- Completion message ------------------------ #
echo && echo "Completed."
