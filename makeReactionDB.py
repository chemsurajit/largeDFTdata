import argparse
from email import header
import json
import logging
from operator import index
import os.path
import sys
import time
import concurrent.futures as confut
import csv
import numpy as np
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
        help="CSV file containing DFT & XTB energies of the QM9 molecules."
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
    parser.add_argument(
        "-json_id_file", "--json_id_file",
        type=str,
        required=True,
        help="Node JSON file. Containing indices to be read in a particular Node."
    )
    return parser.parse_args()


def load_csv_data(rid_csv, json_id_file, mol_data_csv):
    """
    This function load the csv data for molecular atomization energies, and
    reactant, pdt indices from the corresponding csv files: mol_data, rid_csv.
    It will return the whole mol_data and part of the reaction indices which is
    specific to a node.
    return type: pandas dataframe
    """
    # first load the molecular data.
    # The csv file has a column name "index"
    logging.info("Reading molecular data from: %s ..." % mol_data_csv)
    mol_data_pd = pd.read_csv(
        mol_data_csv, keep_default_na=False, na_values=np.nan
    ).dropna()
    # Now load the reaction data from the indices of the json file and rid_csv
    logging.info("Done reading molecular data.")
    logging.info("Reading the Node specific reaction indices...")
    with open(json_id_file) as fp:
        json_dict = json.load(fp)
    index_list = json_dict["indices"]
    del json_dict
    total_rean_pd = pd.read_csv(
        rid_csv, keep_default_na=False, na_values=np.nan
    ).dropna()
    selected_df = total_rean_pd.iloc[index_list]
    return mol_data_pd, selected_df


def check_create_outdir(out_dir):
    """
    This function will check if the out_dir exists already.
    If it exists, the program will exit.
    """
    if os.path.isdir(out_dir):
        logging.warning("The output data directory %s exists." % out_dir)
        logging.warning("If you want to run this program, please delete it.")
    else:
        try:
            os.mkdir(out_dir)
        except OSError as error:
            print(error)
    return


def process_reaction_data(rids_pd, 
                          coreno, 
                          nodeno, 
                          molecule_data_pd, 
                          outdir):
    """
    The main function for the parallel run where all the reactions will be computed.
    return: Integer 0 upon completion.
    """
    pid = os.getpid()
    ppid = os.getppid()
    logging.info("pid, ppid info: %s %s" % (pid, ppid))
    logging.debug("Inside the process_reaction_data function")
    # print the basic info with log
    start_rean_index = list(rids_pd.index.values)[0]
    end_rean_index = list(rids_pd.index.values)[-1]
    logging.info("Start index: %d, pid: %d" % (start_rean_index, pid))
    logging.info("End index: %d, pid: %d" % (end_rean_index, pid))
    logging.info("pid: %d, nreaction to be processed: %d" % (pid, rids_pd.shape[0]))
    #
    output_csv_file = os.path.join(outdir, "Reactions_" + str(nodeno) + "_core_" + str(coreno) + ".csv")
    output_db_file = os.path.join(outdir, "Reactions_" + str(nodeno) + "_core_" + str(coreno) + ".db")
    # This information is important information as it might save time later
    # when we want to see a particular reaction by index.
    # here, startindex, endindex corresponds to index of the reaction
    logging.info("startindex, endindex, csvfile, dbfile: %s %s %s %s" % (start_rean_index, end_rean_index, output_csv_file, output_db_file))
    # First get all the columns corresponds to energies
    energy_columns = molecule_data_pd.columns[molecule_data_pd.columns.str.endswith("_SZ")].to_list() + \
                    molecule_data_pd.columns[molecule_data_pd.columns.str.endswith("_DZP")].to_list() + \
                    molecule_data_pd.columns[molecule_data_pd.columns.str.endswith("_TZP")].to_list() + \
                    ["GFNXTB"]
    logging.info("The energy columns are: %s" % energy_columns)
    # set the start time for this core:
    start_core_time = time.time()
    counter = 0
    chunk_tocsv = [] # chunk data for the csv file writing.
    # for db file, has to be written for every time in loop.
    for rowid, row in rids_pd.iterrows():
        logging.debug("loopstart pid, rowid, row: %d, %d, %s" % (pid, rowid, row.to_string()))
        reactant_index = row.reactindex
        pdt_index = row.pdtindex
        logging.debug("pid, reactant_index: %d %d" %(pid, reactant_index))
        logging.debug("pid, pdt index: %d %d" % (pid, pdt_index))
        try:
            reactant_row = molecule_data_pd.loc[molecule_data_pd['index'] == reactant_index]
        except Exception as ex:
            print("The following exception occurs: %s " % str(ex))
            print("No row in molecule_data_pd for pid, reactant index: %d %d" % (pid, reactant_index))
            continue
        logging.debug("pid, react row: %d %s" % (pid, reactant_row.to_string()))
        try:
            pdt_row = molecule_data_pd.loc[molecule_data_pd['index'] == pdt_index]
        except Exception as ex:
            print("The following exception occurs: %s " % str(ex))
            print("No row in molecule_data_pd for pid, pdt index: %d %d" % (pid, pdt_index))
            continue
        logging.debug("pid, pdt row: %d %s" % (pid, pdt_row.to_string()))
        react_smi = reactant_row["smiles"].values[0]
        logging.debug("pid, react_smi %d %s" % (pid, react_smi))
        pdt_smi = pdt_row["smiles"].values[0]
        logging.debug("pid, pdt_smi: %d %s" % (pid, pdt_smi))
        reaction_properties = pdt_row[energy_columns] - reactant_row[energy_columns].values
        logging.debug("pid, reaction_prop_diff1: %d %s" % (pid, reaction_properties))
        reaction_properties["react_smi"], reaction_properties["pdt_smi"] = [react_smi, pdt_smi]
        logging.debug("pid, reaction_prop_diff2: %d, %s" % (pid, reaction_properties.to_string()))
        #reaction_properties["chemformula"] = reactant_row["chemformula"].values[0]
        logging.debug("pid, reaction_prop_diff4: %d, %s" % (pid, reaction_properties.to_string()))
        reaction_properties["reactindex"], reaction_properties["pdtindex"] = [reactant_index, pdt_index]
        logging.debug("pid, reaction_prop_diff5: %d, %s" % (pid, reaction_properties.to_string()))
        logging.debug("loopend pid, rowid: %d, %d" % (pid, rowid))
        # Now do the hard part.
        # Check how to include two coordinates in ASE db format.
        # 1) Get the coordinates, atoms, in the same order from the xyz files
        #reactant_atoms, reactant_coordinates = get_coordinate_from_xyz(reactant_index, xyzdir)
        #pdt_atoms, pdt_coordinates = get_coordinates_from_xyz(pdt_index, xyzdir)
        # Then include the reactant coordinates and pdt coordinates.
        #
        # reactant_coords, pdt_coords, atoms = (reactant_xyz, pdt_xyz)
        #
        #
        # For now, copy from the bond dependency paper. Add the db format later
        # No need to go for if else in case of db file.
        #The if else will be somewhat cheaper for csv file.
        if counter == 0:
            logging.info("New csv file will created: %s, pid: %d" % (output_csv_file, pid))
            logging.debug("csv, pid, rowid in if: %s, %d, %d" %(output_csv_file, pid, rowid))
            reaction_properties.to_csv(output_csv_file,
                                       mode="w", index=False,
                                       quoting=csv.QUOTE_MINIMAL,
                                       sep=",")
            counter += 1
            continue
        chunk_tocsv.append(reaction_properties)
        if (counter+1) % 10000 == 0:
            logging.info("Converted: %d reactions to %s with pid %d" % (counter, output_csv_file, pid))
            logging.info("Will update the datachunk to csv file: %s" % output_csv_file)
            pd.concat(chunk_tocsv).to_csv(output_csv_file,
                                       mode="a", index=False,
                                       quoting=csv.QUOTE_MINIMAL,
                                       header=False, sep=",")
            chunk_tocsv = []
        counter += 1
    # Now write the remainder data at the end of forloop to the csv file:
    logging.info("Will update the remainder of the datachunk outside for loop to %s, pid: %d" % (output_csv_file, pid))
    pd.concat(chunk_tocsv).to_csv(output_csv_file,
                                mode="a", index=False,
                                quoting=csv.QUOTE_MINIMAL,
                                header=False, sep=",")
    stop_core_time = time.time()
    completed_in = round((stop_core_time-start_core_time)/3600.0, 2)
    logging.info("Loop with pid %d completed in %s hr" % (pid, completed_in))
    return

def main():
    # multi threading credit:
    # https://betterprogramming.pub/pandas-how-to-process-a-dataframe-in-parallel-make-pandas-lightning-fast-669978cf5356
    args = get_arguments()
    log_level = args.logging.upper()
    logging.basicConfig(
        format="[%(asctime)s] %(levelname)s: %(message)s",
        level=log_level,
        datefmt="%H:%M:%S",
    )
    # The directory will be created just to avoid flooding with multiple file
    # inside the user provided location.
    out_dir = os.path.join(os.path.abspath(args.out_dir), "reaction_data")
    # The check create outdir function will check if the directory already
    # exists. If it exists, the program will stop.
    check_create_outdir(out_dir)
    logging.info("Output databases will be saved in %s." % out_dir)
    nprocs = args.nprocs
    logging.info("Number of processors chose: %s" % nprocs)
    # The node number is derived from the JSON file name.
    # The file names are named as: Node_<n>.json, where <n> = 1..N
    # for a N node parallel job submission.
    node_no = args.json_id_file.split(".")[0]
    # load the reactions.csv file as a dataframe
    logging.info("loading indices for %s from %s" % (node_no, args.json_id_file))
    logging.info("loading data...")
    molecule_data_pd, bigchunk_rid_pd = load_csv_data(rid_csv= args.rid_csv,
                                                      json_id_file = args.json_id_file,
                                                      mol_data_csv=args.mol_data)
    # Now split the node specific indices to the number of processors per node:
    splitted_rid_pd = np.array_split(bigchunk_rid_pd, nprocs)
    # Since the csv files are large, try to free as much memory as possible.
    del bigchunk_rid_pd
    # Now come the multi-threading part
    main_func_result = []
    start = time.time()
    logging.info("Starting parallel run in Node: %s" % node_no)
    with confut.ProcessPoolExecutor(max_workers=nprocs) as executor:
        results = [executor.submit(process_reaction_data, rid_pd, coreno, node_no, molecule_data_pd, out_dir)
                   for coreno, rid_pd in enumerate(splitted_rid_pd)]
        for result in confut.as_completed(results):
            try:
                main_func_result.append(result.result())
            except Exception as ex:
                # taken from: https://www.codegrepper.com/code-examples/python/python+exception+with+line+number
                logging.error("Exception occurs: %s" % str(ex))
                ex_type, ex_obj, ex_trace = sys.exc_info()
                line_number = ex_trace.tb_lineno
                logging.error("ERROR: %s; LINE NO: %s" % (str(ex), line_number))
                pass
    end = time.time()
    logging.info("JOB COMPLETED.")
    logging.info("PPID %s completed in %se hr" % (os.getpid(), round((end-start)/3600.0, 2)))
    pass


if __name__ == "__main__":
    main()
