import json
import argparse
import numpy as np
import pandas as pd


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Script to generate indices of the large reactions.csv file according to number of nodes."
    )
    parser.add_argument(
        "--reaction_csv",
        required=True,
        type=str,
        help="csv file name of the file where reactant indices, and pdt indices are saved."
    )
    parser.add_argument(
        "--nnode",
        required=False,
        default=1,
        type=int,
        help="Number of nodes to be used. Default 1"
    )
    return parser.parse_args()


def split_list(indices, nnode):
    """
    split the list of indices according to the number of noes.
    :param indices: list of indices of the reactions.csv file
    :param nnode: Number of nodes
    :return: list of lists containing the indices
    """
    #link: https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
    k, m = divmod(len(indices), nnode)
    return (indices[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(nnode))


def main():
    args = get_arguments()
    nnode = args.nnode
    reactions_pd = pd.read_csv(args.reaction_csv,
                               keep_default_na=False, na_values=np.nan
                               ).dropna()
    indices = reactions_pd.index.to_list()
    splitted_indices = list(split_list(indices, nnode))
    for node in range(nnode):
        out_json = "Node_" + str(node+1) + ".json"
        indices = splitted_indices[node]
        with open(out_json, 'w') as oj:
            jdata = {
                "indices" : indices
            }
            json.dump(jdata, oj, indent=1)
    return


if __name__ == "__main__":
    main()
