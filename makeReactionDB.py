import argparse
import pandas as pd

"""
This is a program to generate database for reaction from the QM9 dataset. 
Or dataset similar to that. The database will contain reactant index,
product index, reactant coordinates, product coordinates, reaction SMILES,
reactant smile, product smile, 76 different reaction energies, GFNXTB reaction
energies. 
"""


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Make sqlite3 format db file using ASE for the reactions."
    )
    parser.add_argument(
        "-rid_csv", "--rid_csv",
        type=str,
        required=True,
        help="CSV file containing reactant and product indices."
    )
    parser.add_argument(
        "-mol_data", "--mol_data",
        type=str,
        required=True,
        help="DB file containing DFT & XTB energies of the QM9 molecules."
    )
    parser.add_argument(
        "-nprocs", "--nprocs",
        type=int,
        required=False,
        default=1,
        help="Number of processors to be used."
    )
    parser.add_argument(
        "-out_dir", "--out_dir",
        type=str,
        required=False,
        default="outputs",
        help="Directory where the Reactions_n.csv files will be saved."
    )
    parser.add_argument(
        "-logging", "--logging",
        type=str,
        required=False,
        default="info",
        choices=["debug", "info", "warning", "error", "critical"],
        help="Provide logging level. Default is warning."
    )
    return

def main():
    pass


if __name__ == "__main__":
    main()
