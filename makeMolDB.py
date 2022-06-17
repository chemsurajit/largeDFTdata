import argparse
import os
import re
import sys
import fnmatch
from collections import Counter
import ase.db
import xyz2mol
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import pandas as pd
import csv
import logging


def check_create_outdir(outdir=None):
    """
    This function will create a output directory if it does not exists.
    Will not create if the outdirectory exists.
    """
    if outdir is not None:
        if os.path.isdir(outdir):
            logging.info("Directory %s exists. Will not create new." % outdir)
        else:
            os.mkdir(outdir)
    else:
        logging.error("variable outdir can't be None")
    return


def get_files(directory, match=""):
    """This function will return a list of files as list. The list elements
    will be in the form of path object. The filenames will be searched with the
    match string.
    args:
        directory - Path from which files need to be found.
        match - a matching string for the files. Else, all the files will be returned.
    returns:
        files - List of files with names as match string after being sorted."""
    allfiles = []
    for files in os.listdir(directory):
        # Take the xyz files which has 'dsgdb9nsd_' in their name
        if fnmatch.fnmatch(files, match):
            allfiles.append(os.path.abspath(os.path.join(directory, files)))
    # Check if the list xyz_files is empty. If it is empty, then program exits.
    if not allfiles:
        print("No files with match %s found in directory: %s " % (match, directory))
        print("Program exit now.")
        sys.exit()
    return sorted(allfiles)


def get_dft_energies(logfile):
    """This function takes a ADF logfile and return the dft energies as dictionary.
    args:
        logfile - ADF output file
    returns:
        dft_energies - dictionary of DFT energies. Keys as functional, value as energy in eV.
    """
    dft_energies = {}
    with open(logfile, 'r') as flog:
        for line in flog:
            if "FR:" in line:
                FR_cont = True
                func = line[4:19].strip().upper()
                en_ev = float(line.split("=")[1].split()[1])
                dft_energies[func] = en_ev
    if not FR_cont:
        logging.warning("No data found for DFT functional in: %s" % logfile)
        logging.info("The reason could be that the ADF log file prints the energy with different format.")
    return dft_energies


def mol_to_bonds_list(molobj, filename):
    """This function takes a RDKit mol object as input with xyzfile and makes
    a list of bonds. The list of bonds will always be confined in the variable
    defined as bonds_list. If some bonds not present, it will be assigned to zero.
    args:
        molobj - RDKit format mol object
        filename - The xyz file name
    returns:
        bonds_list - list of bonds after counting each of the bond types.
    """
    bonds_list = {'C_s_C': 0, 'C_d_C': 0, 'C_t_C': 0, 'C_C_A': 0,
                  'C_s_H': 0, 'C_s_O': 0, 'C_d_O': 0, 'C_t_O': 0,
                  'C_O_A': 0, 'C_s_N': 0, 'C_d_N': 0, 'C_t_N': 0,
                  'C_N_A': 0, 'C_s_F': 0, 'O_s_O': 0, 'O_d_O': 0,
                  'O_O_A': 0, 'O_s_H': 0, 'O_s_N': 0, 'O_d_N': 0,
                  'O_t_N': 0, 'O_N_A': 0, 'O_s_F': 0, 'N_s_N': 0,
                  'N_d_N': 0, 'N_t_N': 0, 'N_N_A': 0, 'N_s_H': 0,
                  'N_s_F': 0, 'F_s_H': 0
                  }

    xyzfile = os.path.basename(filename)
    if molobj is not None:
        for bond in molobj.GetBonds():
            btype = bond.GetBondType()
            a1 = molobj.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()
            a2 = molobj.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()
            #
            # C-C bonds
            #
            if a1+a2 == "CC":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_C'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['C_d_C'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['C_t_C'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['C_C_A'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
                    raise ValueError
            #
            # C-H bonds
            #
            elif a1+a2 == "CH" or a2+a1 == "CH":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_H'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
                    raise ValueError
            #
            # C-O bonds
            #
            elif a1+a2 == "CO" or a2+a1 == "CO":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_O'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['C_d_O'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['C_t_O'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['C_O_A'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # C-N bonds
            #
            elif a1+a2 == "CN" or a2+a1 == "CN":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_N'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['C_d_N'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['C_t_N'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['C_N_A'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # C-F bonds
            #
            elif a1+a2 == "CF" or a2+a1 == "CF":
                if str(btype) == "SINGLE":
                    bonds_list['C_s_F'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # O-O bonds
            #
            elif a1+a2 == "OO":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_O'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['O_d_O'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['O_O_A'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # O-H bond
            #
            elif a1+a2 == "OH" or a2+a1 == "OH":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_H'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # O-N bonds
            #
            elif a1+a2 == "ON" or a2+a1 == "ON":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_N'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['O_d_N'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['O_t_N'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['O_N_A'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # O-F bond
            #
            elif a1+a2 == "OF" or a2+a1 == "OF":
                if str(btype) == "SINGLE":
                    bonds_list['O_s_F'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # N-N bonds
            #
            elif a1+a2 == "NN":
                if str(btype) == "SINGLE":
                    bonds_list['N_s_N'] += 1
                elif str(btype) == "DOUBLE":
                    bonds_list['N_d_N'] += 1
                elif str(btype) == "TRIPLE":
                    bonds_list['N_t_N'] += 1
                elif str(btype) == "AROMATIC":
                    bonds_list['N_N_A'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # N-H bond
            #
            elif a1+a2 == "NH" or a2+a1 == "NH":
                if str(btype) == "SINGLE":
                    bonds_list['N_s_H'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # N-F bonds
            #
            elif a1+a2 == "NF" or a2+a1 == "NF":
                if str(btype) == "SINGLE":
                    bonds_list['N_s_F'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
            #
            # F-H bond
            #
            elif a1+a2 == "FH" or a2+a1 == "FH":
                if str(btype) == "SINGLE":
                    bonds_list['F_s_H'] += 1
                else:
                    logging.error("Don't know bond type: %s %s" % (str(btype), (a1+"-"+a2)))
                    logging.error("Exiting for: %s" % xyzfile)
    return bonds_list


def get_properties_combined(index, smiles, chem_formula, bonds_list):
    """This function takes an index, smiles string and a list of bonds and return
    a dictionary with the proper key names. This key names will be used homogeneously
    in all the programs that follows.
    args:
        index - index of the molecule in the QM9_GMP2 dataset
        smiles - smiles string as calculated using the xyz2mol program
        chem_formula - Formula as type string.
        bonds_list - A dictionary with all the bond names as keys and their numbers
                     in the molecule as values.
    returns:
        dictionary with keys as below.
    """
    all_properties = {"index": index, "smiles": smiles, "chemformula":chem_formula}
    return {**all_properties, **bonds_list}


def get_smiles_from_xyz(ixyzfile):
    """Function that takes xyzfile path and return the smiles"""
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(ixyzfile)
    return_code = True
    mols = None
    try:
        mols = xyz2mol.xyz2mol(atoms, xyz_coordinates, charge=0,
                               use_graph=True, allow_charged_fragments=False,
                               embed_chiral=True, use_huckel=False)
    except:
        mols = None
        logging.warning("No mol object from the xyz file: %s" % ixyzfile)
    return mols


def update_failed_indices(outputfile, indices):
    if not indices:
        logging.info("Number of failed indices = %d" % len(indices))
        logging.info("No file will be created.")
    else:
        with open(outputfile, 'w') as fp:
            for index in indices:
                fp.write("%s\n" % str(index))
    return


def update_pd_df(inpdf, index=None, smiles=None, chemformula=None, bonds=None, g4mp2_energy=None, energies=None):
    """Function to update to dataframe
    """
    row_dict = {"index":index, "smiles":smiles, "chemformula":chemformula, "G4MP2":g4mp2_energy}
    row_dict.update(energies)
    row_dict.update(bonds)
    df = pd.DataFrame(row_dict, index=[0])
    newdf = pd.concat([inpdf, df])
    return newdf


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Make db file of the QM9 dataset with 76 DFT and 3 basis sets."
    )
    parser.add_argument(
        "-xyz_dir", "--xyz_dir",
        type=str,
        required=True,
        help="Location of the xyz directory. Default is current directory."
                        )
    parser.add_argument(
        "-dft_log_dir", "--dft_log_dir",
        type=str,
        required=True,
        help="Location of the directory from where the logfiles will be read."
    )
    parser.add_argument(
        "-output_dir", "--output_dir",
        type=str,
        required=False,
        default="output",
        help="Name of the directory where all the outputs from this file will be saved."
    )
    parser.add_argument(
        "-logging", "--logging",
        type=str,
        required=False,
        default="warning",
        choices=["debug", "info", "warning", "error", "critical"],
        help="Provide logging level. Default is warning."
    )
    parser.add_argument(
        "-output_format", "--output_format",
        type=str,
        required=False,
        default="csv",
        choices=["csv", "db"],
        help="Output format for the data. db stands for sqlite3 db format using ASE."
    )
    return parser.parse_args()


def get_energies(adf_out_file):
    """
    This function will return the energies of each of the 76 DFT functionals as dictionary.
    dictionary format: {FUNCTIONAL_name : energy_value}
    """
    dft_energies = {}
    FR_cont = False # A flag to print if the outfile does not contain energies with expected format
    with open(adf_out_file) as fp:
        for line in fp.readlines():
            if "FR:" in line:
                FR_cont = True
                func = line[4:19].strip().upper()
                en_ev = float(line.split("=")[1].split()[1])
                dft_energies[func] = en_ev
    if not FR_cont:
        logging.error("ADF outfile %s doesn't contain energies in required form." % adf_out_file)
    return dft_energies

def create_atoms_csv(dft_log_dir):
    """
    This function will create atoms.csv files inside the following directories:
    PATH/SZ/atoms/, PATH/DZP/atoms/, and PATH/TZP/atoms/
    It is also assumed that the extension of the adf log files are in .out format, ie,
    C.out, F.out, H.out, N.out, O.out
    """
    atoms = ["C", "F", "H", "N", "O"]
    atom_dirs = [
        os.path.join(dft_log_dir, "SZ", "atoms"),
        os.path.join(dft_log_dir, "DZP", "atoms"),
        os.path.join(dft_log_dir, "TZP", "atoms")
    ]
    for atom_dir in atom_dirs:
        atoms_data = []
        out_csv = os.path.join(atom_dir, "atoms.csv")
        if os.path.isfile(out_csv):
            logging.warning("atoms csv file %s exists. Skipping." % out_csv)
        else:
            for atom in atoms:
                adf_atom_out = os.path.join(atom_dir, atom+".out")
                energies = get_energies(adf_atom_out)
                energies["atom"] = atom
                atoms_data.append(energies)
            with open(out_csv, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=atoms_data[0].keys())
                writer.writeheader()
                for row in atoms_data:
                    writer.writerow(row)
            logging.info("atoms.csv file created inside: %s" % atom_dir)
    return

def create_molecules_csv(dft_log_dir):
    """
    This function will create molecules.csv files inside the following directories:
    <PATH>/SZ/molecules/, <PATH>/DSP/molecules/, <PATH>/TZP/molecules/
    It is also assumed that the extension of the adf log files are in .out format, ie,
    mol1.out, mol2.out, etc.
    """
    molecule_dirs = [
        os.path.join(dft_log_dir, "SZ", "molecules"),
        os.path.join(dft_log_dir, "DZP", "molecules"),
        os.path.join(dft_log_dir, "TZP", "molecules")
    ]
    for molecule_dir in molecule_dirs:
        out_csv = os.path.join(molecule_dir, "molecules.csv")
        logging.info("Creating molecules.csv file.... ")
        logging.debug("creating molecules.csv file at: %s" % molecule_dir)
        if os.path.isfile(out_csv):
            logging.info("molecules csv file %s exist. Skipping." % out_csv)
        else:
            # The adf output file names should to same as match.
            logfiles = get_files(molecule_dir, match="*_xyz.out")
            loop_counter = 0
            with open(out_csv, 'w') as csvfile:
                for logfile in logfiles:
                    energies = get_energies(logfile)
                    logindex = int(os.path.basename(logfile).split("_")[0])
                    if not energies:
                        logging.warning("No dft energies found for index: %s" % logindex)
                        continue
                    energies["index"] = logindex
                    if loop_counter == 0:
                        writer = csv.DictWriter(csvfile, fieldnames=energies.keys())
                        writer.writeheader()
                    writer.writerow(energies)
                    loop_counter += 1
        logging.info("Done creating molecule csv file: %s" % out_csv)
    return


def get_xyz_info(xyzfile):
    """
    This function will return information from xyz formated coordinate file.

    It assume the file name to be dsgdb9nsd_<i>.xyz where <i> is the index in QM9_G4MP2 set.
    """
    index = int(os.path.basename(xyzfile).split(".")[0].split("_")[1])
    atoms = []
    coords = []
    with open(xyzfile) as fp:
        natoms = int(fp.readline())
        comment = fp.readline()
        for i in range(natoms):
            line = fp.readline()
            atoms.append(line.split()[0])
            coords.append([float(i) for i in line.split()[1:]])
    return index, atoms, coords


def get_atomization_energy(mol_e_row, atoms_energies, atoms_list):
    """
    This function will compute the atomization energies of the molecules from the molecules
    and atoms energies from the corresponding molecules.csv and atoms.csv files.

    mol_e_row is a pd row object.
    atoms_energies is the whole pd dataframe for all the atoms
    """
    atoms_count = Counter(atoms_list)
    mol_e_row = mol_e_row.set_index('index')
    mol_e_dict = mol_e_row.to_dict('records')[0]
    for atoms, natoms in atoms_count.items():
        atoms_energies_dict = atoms_energies.loc[atoms_energies["atom"] == atoms].drop(["atom"],
                                                        axis=1).to_dict('records')[0]
        mol_e_dict = {key: mol_e_dict[key] - natoms*atoms_energies_dict.get(key, 0)
                      for key in mol_e_dict.keys()}
    return mol_e_dict


def modify_dict_keys(atomization_en_dict):
    """
    This function will modify the functional names by replacing the numbers,
    etc and add the corresponding "SZ" or "DZP" or "TZP" suffix to the names.
    Also, this function changes all non-alphanumeric characters of the functionals
    to underscore which is supported in the ASE db format.
    """
    modified_dict = {}
    for func, en in atomization_en_dict.items():
        new_key = re.sub("[^0-9a-zA-Z]+", "_", func)
        modified_dict[new_key] = en
    return modified_dict


def create_molecules_db(outdb_file, xyzfiles, failed_indices_outfile, dft_log_dir):
    """
    This function creates the sqlite3 db file using ASE in the output directory.

    The molecules.csv and atoms.csv files will be read from the dft_log_dir directory.

    The failed_indices is a file where the indices of the molecules will be saved for
    which the xyz2mol program fails to convert xyz to RDKit format mol object.

    energies are saved as external table in ASE database format.
    The structure is: {"SZ": {"func1":E1, "func2":E2...}, "DZP":{"func1":E1, "func2":E2...}..}
    """
    failed_mols_indices = []
    logging.info("creating db file...")
    logging.info("loading atomic csv files...")
    sz_atoms_energies = pd.read_csv(
        os.path.join(dft_log_dir, "SZ", "atoms", "atoms.csv"), index_col=False
    )
    dzp_atoms_energies = pd.read_csv(
        os.path.join(dft_log_dir, "DZP", "atoms", "atoms.csv"), index_col=False
    )
    tzp_atoms_energies = pd.read_csv(
        os.path.join(dft_log_dir, "TZP", "atoms", "atoms.csv"), index_col=False
    )
    logging.info("loading molecules csv files...")
    sz_mol_pd = pd.read_csv(
        os.path.join(dft_log_dir, "SZ", "molecules", "molecules.csv")
    )
    dzp_mol_pd = pd.read_csv(
        os.path.join(dft_log_dir, "DZP", "molecules", "molecules.csv")
    )
    tzp_mol_pd = pd.read_csv(
        os.path.join(dft_log_dir, "TZP", "molecules", "molecules.csv")
    )
    row_count = 0
    with ase.db.connect(outdb_file, append=False) as asedb:
        for xyzfile in xyzfiles:
            external_table = {}
            key_val_pairs = {}
            index, atoms, coords = get_xyz_info(xyzfile)
            mols = get_smiles_from_xyz(xyzfile)
            if mols is None:
                logging.warning("Failed to convert xyz to RDkit mols object: %s" % xyzfile)
                failed_mols_indices.append(index)
                continue
            #bonds_list = mol_to_bonds_list(mols[0], xyzfile)
            smiles = Chem.MolToSmiles(mols[0], isomericSmiles=True)
            #chemformula = CalcMolFormula(mols[0])
            sz_e_row = sz_mol_pd.loc[sz_mol_pd["index"] == index]
            dzp_e_row = dzp_mol_pd.loc[dzp_mol_pd["index"] == index]
            tzp_e_row = tzp_mol_pd.loc[tzp_mol_pd["index"] == index]
            sz_atomization_e_dict = get_atomization_energy(sz_e_row, sz_atoms_energies, atoms)
            dzp_atomization_e_dict = get_atomization_energy(dzp_e_row, dzp_atoms_energies, atoms)
            tzp_atomization_e_dict = get_atomization_energy(tzp_e_row, tzp_atoms_energies, atoms)
            # modify the keys of the functional name so that it can be updated to ase db
            external_table["SZ"] = modify_dict_keys(sz_atomization_e_dict)
            external_table["DZP"] = modify_dict_keys(dzp_atomization_e_dict)
            external_table["TZP"] = modify_dict_keys(tzp_atomization_e_dict)
            key_val_pairs["index"] = index
            key_val_pairs["smiles"] = smiles
            ase_atoms_obj = ase.Atoms(symbols=atoms, positions=coords, pbc=False)
            asedb.write(ase_atoms_obj,
                        key_value_pairs=key_val_pairs,
                        external_tables=external_table)
            #asedb.write(ase_atoms_obj, external_tables=external_table)
            row_count += 1
            if (row_count % 10000) == 0:
                logging.info("Finished %d rows." % row_count)
    logging.info("Finished creating db file: %s", outdb_file)
    update_failed_indices(failed_indices_outfile, failed_mols_indices)
    return


def make_db_mols_ens(
        output_dir,
        xyzfiles,
        failed_indices_outfile,
        dft_log_dirs,
        output_format
):
    if output_format == "csv":
        print("call for csv file making")
    elif output_format == "db":
        output_file = os.path.join(output_dir, "molecules_qm9.db")
        create_molecules_db(output_file, xyzfiles, failed_indices_outfile, dft_log_dirs)
    pass


def main():
    args = get_arguments()

    # logging setup
    log_level = args.log.upper()
    logging.basicConfig(
        format="[%(levelname)s: %(message)s",
        level=log_level,
        datefmt="%H:%M:%S",
    )

    xyz_dir = os.path.abspath(args.xyz_dir)
    dft_log_dir = os.path.abspath(args.dft_log_dir)
    output_dir = os.path.abspath(args.output_dir)
    # create the output directory if doesn't exists
    check_create_outdir(output_dir)
    failed_indices_outfile = os.path.join(output_dir, "failed_xyz2mol.dat")
    #outdb_file = os.path.join(output_dir, "molecules_qm9.db")
    xyzfiles = get_files(xyz_dir, match="dsgdb9nsd_*.xyz")
    logging.info("Number of xyz files: %d" % len(xyzfiles))
    # First, create atoms.csv files inside each of the SZ, DZP, and TZP directories.
    create_atoms_csv(dft_log_dir)
    create_molecules_csv(dft_log_dir)
    make_db_mols_ens(
        output_dir,
        xyzfiles,
        failed_indices_outfile,
        dft_log_dir,
        args.output_format
                     )
    #create_molecules_db(output_file, xyzfiles, failed_indices_outfile, dft_log_dir)
    return


if __name__ == "__main__":
    main()
