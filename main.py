#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 13-11-2025

"""

import argparse as arg


def arg_reader():
    """ Reads arguments from command line

    :return class containing the arguments
    """
    arg_parser = arg.ArgumentParser(
        description="Runs GEMMA and downstream processing"
    )
    arg_parser.add_argument(
        "phenotype_file",
        help="The phenotype file containing sample name as first column and " \
        "phenotypes on subsequent columns"
    )
    arg_parser.add_argument(
        "bed_files",
        help="Prefix for the bed files, fam will be " \
        "modified to include phenotypes"
    )
    arg_parser.add_argument(
        "kinship_matrix",
        help="Kinship matrix for population structure correction"
    )

    gemma_options = arg_parser.add_argument_group(
        title="GEMMA optional parameters",
        description="These parameters are fed directly into GEMMA, see GEMMA " \
        "manual for further details."
    )

    gemma_options.add_argument(
        "--maf",
        help="Minor allele frequency threshold, default = 0.05",
        type = float, 
        default=0.05
    )
    gemma_options.add_argument(
        "--miss",
        help="Max fraction of missing calls per variant, default = " \
        "0.01",
        type=float,
        default=0.01
    )
    gemma_options.add_argument(
        "--lmm",
        help="Linear mixed model, default = 1",
        type=int,
        default=1
    )
    return arg_parser


def main():
    args = arg_reader().parse_args()

    # Merge phenotypes into fam file
    

    # Run GEMMA

    # Make manhattan plots 


if __name__ == "__main__":
    main()
