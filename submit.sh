#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --output=log_%x.log
#SBATCH --job-name=test
#SBATCH --time=2-2:00:00
# for node fixing
#SBATCH --exclusive # for exclusively one node
#SBATCH --partition=xeon40
#SBATCH -N 1 # minimum of 1 nodes
#SBATCH -n 40 # 40 mpi processor total for one node

ulimit -s unlimited

# setting up of anaconda environment. Optional
source /home/energy/surna/anaconda3/etc/profile.d/conda.sh
conda activate sure_svol

pyscript=$1
verbose="info" # change the verbose to other level for parallel python

#The following parameters are set to default.
# The values will only change if number of input arguments are > 2.
reactions_csv="./test_small_data/reactions.csv"
mol_csv="./Data/molecules_qm9.csv"
output_dir="./Data/"
nprocs=1
indices="Node_0.json"

if [[ $# -gt 2 ]]; then
	reactions_csv=$2
	mol_csv=$3
	output_dir=$4
	nprocs=$5
	indices=$6
    verbose=$7
fi # To understand this, refer to the automated running bash script.

echo "In submit script, parallel python reaction: $pyscript"
echo "In submit script, reactions_csv: $reactions_csv"
echo "In submit script, molecules_csv: $mol_csv"
echo "In submit script, output_dir : $output_dir"
echo "In submit script, verbose : $verbose"
echo "In submit script, indices file: $indices"
echo "In submit script, nprocs: $nprocs"



python3 $pyscript --rid_csv $reactions_csv \
    --mol_data $mol_csv \
    --out_dir $output_dir --logging $verbose \
    --nprocs $nprocs --json_id_file $indices 

echo "finished"
