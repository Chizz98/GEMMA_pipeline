#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 13-11-2025

"""

import argparse as arg
import multiprocessing
import scripts.merge_fam as merge_fam
import scripts.run_gemma as run_gemma
import sys
import os
import shutil


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

    arg_parser.add_argument(
        "-p", "--parallel",
        help="Number of parallel processes for GEMMA, default = 1",
        type=int,
        default=1
    )
    arg_parser.add_argument(
        "-o", "--out_dir",
        help="Output directory, default = gemma_pipe_out",
        type=str,
        default="gemma_pipe_out"
    )
    return arg_parser


def confirm_parralel(num_parallel: int) -> bool:
    return_val = True
    cpu_count = multiprocessing.cpu_count()
    if num_parallel > multiprocessing.cpu_count():
        raise Exception(
            f"parralel argument exceeds number of available CPUs ({cpu_count})"
        )
    elif num_parallel > 0.5 * cpu_count:
        print("Warning: parallel processes exceeds half of available CPUs "
        "({cpu_count})")
        response = input("Continue? ([Y/N])").strip().lower()
        if response not in ("y", "yes"):
            return_val = False
    return return_val


def out_dir_handler(out_dir: str) -> bool:
    return_val = True
    if os.path.exists(out_dir):
        print("Warning: out dir already exists")
        response = input("Continue? ([Y/N])").strip().lower()
        if response not in ("y", "yes"):
            return_val = False
    else:
        os.mkdir(out_dir)
    return return_val


def gemma_worker(gemma_args):
    run_gemma.main(
        (
            gemma_args["script"],
            gemma_args["bfiles"],
            gemma_args["kinship"],
            gemma_args["maf"],
            gemma_args["missing"],
            gemma_args["lmm"],
            gemma_args["phenotype"],
            gemma_args["pheno_col"]
        )
    )

    phenotype = gemma_args["phenotype"]
    default_out = "output"
    out_dir = gemma_args["out_dir"]

    # Move output files from default out dir to specified subfolder
    output_files = [
        file_name for file_name in os.listdir(default_out) if 
        file_name.startswith(phenotype)
        ]

    source = [
        os.path.join(default_out, file_name) for file_name in output_files
        ]
    dest = [
        os.path.join(out_dir, file_name) for file_name in output_files
        ]
    for source_fn, dest_fn in zip(source, dest):
        shutil.move(source_fn, dest_fn)


def main():
    args = arg_reader().parse_args()

    # Check parallel processes against available CPU
    num_parallel=args.parallel
    parralel_confirmed = confirm_parralel(num_parallel)
    if not parralel_confirmed:
        sys.exit(1)


    # Check GEMMA default out_dir
    if os.path.exists("output"):
        raise Exception(
            "./output dir exists, GEMMA uses this as default output. Run " \
            "script in dir without ./output present."
        )

    # Make out dir
    out_dir = args.out_dir
    out_dir_handled = out_dir_handler(out_dir)
    if not out_dir_handled:
        sys.exit(1)
        

    # Merge phenotypes into fam file
    pheno_file = args.phenotype_file
    fam_file = args.bed_files + ".fam"
    phenotypes = merge_fam.main(
        argv=[
            "merge_fam.py", 
            pheno_file, 
            fam_file
            ]
        )

    # Run GEMMA
    gemma_out_dir = os.path.join(out_dir, "gemma_output")
    if not os.path.exists(gemma_out_dir):
        os.mkdir(gemma_out_dir)

    bed_files = args.bed_files
    kinship_mat = args.kinship_matrix
    maf_th = args.maf
    missing_th = args.miss
    lmm = args.lmm
    gemma_arg_list = []
    for i, phenotype in enumerate(phenotypes):
            gemma_args = {
                "script": "run_gemma.py",
                "bfiles": bed_files,
                "kinship": kinship_mat,
                "maf": str(maf_th),
                "missing": str(missing_th),
                "lmm": str(lmm),
                "phenotype": phenotype,
                "pheno_col": i + 1,
                "out_dir": gemma_out_dir
            }
            gemma_arg_list.append(gemma_args)

    num_children = min(num_parallel, len(phenotypes))
    with multiprocessing.Pool(num_children) as pool:
        pool.map(gemma_worker, gemma_arg_list)

    # Remove default GEMMA output folder
    os.rmdir("output")

    # Make manhattan plots 


if __name__ == "__main__":
    main()
