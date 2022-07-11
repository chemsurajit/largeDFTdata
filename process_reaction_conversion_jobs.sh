#!/bin/bash

# This script will prepare and submit the jobs related to reaction data processing.

# first parameter processing
###############################################################################
nnode=1 #Number of nodes; default 1
nprocs=1 # number of processors; default 1
mol_csv="./csvs/qm9_bonds_energies.csv" #default name of the csv file
rean_csv="./csvs/reactions.csv" #default name of the reaction indices csv file
# The pyscript_path is set to the directory where this script is kept by default:
pyscript_path="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";
outdir="./Data"
submit_script="$pyscript_path/submit.sh"
verbose="info"
parallel_main_python_script="makeReactionDB.py"

help () {
  echo "Usage: $0
  Description:
    This script will first run the split_reactionscsv.py script to create json
    files containing index of the reactions.csv file. The number of this file
    is equal to number of nodes. It then runs the make_reactions_parallel.py
    script to create the Reactions_n.csv files under the directory Nodes.

  Optional parameter:
    [-nnodes/--nnodes nnodes] [-nprocs/--nprocs nprocs]
    [-mol_csv/--mol_csv mol_csv] [-rean_csv/--rean_csv rean_csv]
    [-h/--help]

  Parameters:
    -nnode/--nnodes: Number of nodes. Default: ${nnode}
    -nprocs/--nprocs: Number of processors for each nodes. Default: ${nprocs}
    -mol_csv/--mol_csv: CSV file containing QM9 mol data.
            Default: ${mol_csv}
    -rean_csv/--rean_csv: CSV file containing only indices for the reactions.
            Default: ${rean_csv}
    -pyscript_path/--pyscript_path: The PATH where the python scripts are kept.
            Default: ${pyscript_path}
    -verbose/--verbose: The verbosity of parallel python script.
            Default: ${verbose}
    -h/--help: To print this help
  "
}

while [[ "$#" -gt 0 ]]; do
  case $1 in
    -nnode|--nnode) nnode="$2"; shift ;;
    -nprocs|--nprocs) nprocs="$2" ; shift ;;
    -mol_csv|--mol_csv) mol_csv="$2"; shift ;;
    -rean_csv|--rean_csv) rean_csv="$2"; shift ;;
    -pyscript_path|--pyscript_path) pyscript_path="$2"; shift ;;
    -verbose|--verbose) verbose="$2"; shift ;;
    -h|--help) help; exit 0 ;;
    *) echo "Unknown parameter: $1"; help; exit 1 ;;
  esac
  shift
done
# optional parameter processing completed.
###########################################################################

# Sanity test:
# Print out the parameters.
echo "The parameters are set to:"
echo "nnode: $nnode"
echo "nprocs: $nprocs"
echo "mol_csv: $mol_csv"
echo "rean_csv: $rean_csv"
echo "python script path: $pyscript_path"
echo "parallel python verbose: $verbose"

# Step 1. Creation of the Nodes_<n>.json files.
run_json="y"
if ls Node_*.json 1> /dev/null 2>&1; then
  echo "Json file exists."
  read -p "Want to create new json files? (Y/N): " run_json
  if [[ "$run_json" == [yY] ]]; then
    echo "Json file will be created."
    rm Node_*.json
    python3 "$pyscript_path"/split_reactionscsv.py --reaction_csv "$rean_csv" --nnode "$nnode"
  else
    echo "json files will be read."
    njson=$(ls *.json| wc -l)
    if [ $njson -ne $nnode ]; then
      echo "Number of nodes are not equal to number of json files."
      echo "Number of nodes will be set to number of json file"
      nnode=$njson
      echo "Number of nodes set to: $nnode"
    fi
  fi
else
  echo "Running split file."
  python3 "$pyscript_path"/split_reactionscsv.py --reaction_csv "$rean_csv" --nnode "$nnode"
fi


# Step 2. Create the directories according to the number of nodes.
if [ ! -d "$outdir" ]; then
  mkdir -p $outdir
fi


# step 3: submitting jobs.
# Decide whether to remove the csv files (if exists) or to run the script
# a. check if one file exists in the outputs file.
run_calc="y"
if ls $outdir/Node_*_core_*.csv 1> /dev/null 2>&1; then
  echo "Csv file with reaction data exists inside $outdir "
  echo "If you chose Y/y, all the files will be deleted and new calculations will be performed."
  read -p "Want to run new calculations? (Y/N): " run_calc && [[ $run_calc == [yY] ]] || exit 1
else
  echo "Running calculations with $nnode nodes and $nprocs processors."
fi
if [ "${run_calc^^}" == "Y" ]; then
  echo "Running calculations"
  for json in *.json; do
    echo "submitting jobs with indices from: $json"
    echo ${json%.json}
    sbatch -J ${json%.json} $submit_script \
              $pyscript_path/$parallel_main_python_script \
              $rean_csv $mol_csv \
              $outdir $nprocs \
              $json $verbose
  done
fi
