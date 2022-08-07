# DFTbondDependency
Repository for bond dependency paper

## Dependency:
3) pandas
4) csv
5) Numpy
6) requests


All the packages are part of standard Python library and can be installed with either Conda, PIP, etc.

If this repository is used, please cite us. 
1. Citation to the code can be downloaded by clicking: "Cite this repository" in the right side panel.
2. Citation to the preprint in: https://doi.org/10.26434/chemrxiv-2022-9prf3

## Description of the Data.
The data is publicly available now. Link: https://doi.org/10.26434/chemrxiv-2022-9prf3

## Description of the codes.
The repository contains python scripts for the calculations described in the paper: https://doi.org/10.26434/chemrxiv-2022-9prf3
1) get_data.py: This script will download data from the DTU Data website. The public link to the website
will be available upon acceptance of the paper. This file will download either all the files from the 
database (if -all/--all option is given) Or it will only download the xyzfiles and the log files of the
energy calculations

2) makeMolDB.py: This file makes a csv file containing the atomization energy values from the molecule and atoms log files.

3) make_reaction_ids.py: This script will create a csv file with only two columns: 'reactantindex','pdtindex'.
The indices are the index of the molecules in the csv file made by the make_molecule_bond_en_csv.py script.

4) process_reaction_conversion_jobs.sh: This is a bash script to create the final csv file containing
all the information related to the reactions. It takes the csv file containing molecular data (created by the
script makeMolDB.py), the indices of the reactants and products in form of a csv file
(created by using the script make_reaction_ids.py), The G4MP2 energies of the molecules as csv file (with
index and energy), the path of the python script make_reactions_parallel.py, number of Nodes to be used,
and number of processors per each nodes.
It first split the csv file containing indices of the reactants and products according to the number of Nodes
and saves those in a json file with names Node_n.json with n from {1,2,...n} if n number of nodes are used.

5) make_reactions_parallel.py: This file takes csv file containing indices for the reactions, csv file containing
all the data of the molecules, csv file containing G4MP2 energy, number of processors, json file
containing the indices of the csv file with "reactantindex","pdtindex".

6) submit.sh: An example submit script to run make_reactions_parallel.py in a single node with multiple processors.
It is called from the script process_reaction_conversion_jobs.sh. It is written for the slurm scheduler.


10) CITATION.cff: This file is to provide citation data for this repository in bibtex or APA format.

## How to run:
The help message for each of the files (except submit.sh) can be obtained by running the corresponding
script with -h.


The steps described in the paper can be followed by running the below scripts in the following sequence: 
1. get_data.py
2. makeMolDB.py
3. make_reaction_ids.py
4. process_reaction_conversion_jobs.sh, make_reactions_parallel.py, submit.sh

## License:
All the scripts in this repository are covered under the MIT license terms (LICENSE.txt).
