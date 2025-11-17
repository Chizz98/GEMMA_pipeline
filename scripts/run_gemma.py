#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 13-11-2025

Runs GEMMA with specified arguments, supports custom output folder
"""
import subprocess
import os
import sys
import re


def run_gemma(param_dict: dict) -> None:
    # Check if biles exist
    for bfile_ext in [".fam", ".bim", ".bed"]:
        bfile_path = param_dict["bfile"] + bfile_ext
        bfile_exists = os.path.exists(bfile_path)
        if not bfile_exists:
            raise Exception(
                f"{bfile_path} does not exist"
            )

    gemma_cmd_blocks = [
        f"-{param} {val}" for param, val in param_dict.items()
        ]
    gemma_cmd = f"gemma {' '.join(gemma_cmd_blocks)}"
    subprocess.run(
        gemma_cmd, 
        shell=True, 
        text=False, 
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
        )


def parse_log(log_fn: str) -> list[str]:
    with open(log_fn) as infile:
        log_text = infile.read()
    out_dict = {
        "phenotype": os.path.basename(log_fn).split(".")[0],
        "total_ind": re.search(
             r"total individuals = ([0-9]+)", log_text).groups()[0],
        "analyzed_ind": re.search(
             r"analyzed individuals = ([0-9]+)", log_text).groups()[0],
        "total_var": re.search(
             r"total SNPs\/var = ([0-9]+)", log_text).groups()[0],
        "analyzed_var": re.search(
             r"analyzed SNPs\/var = ([0-9]+)", log_text).groups()[0]
    }
    return out_dict
    
    

def main(argv = None):
    argv = argv or sys.argv
    if len(argv) != 8:
            print("Usage: python run_gemma.py <bfiles> <kinship_mat> <maf> "
                  "<missing_threshold> <lmm> <output filename> <phenotype_num>")
            sys.exit(1)
    
    # Make parameter_dictionary
    gemma_params = {
        "bfile": argv[1],
        "k": argv[2],
        "maf": argv[3],
        "miss": argv[4],
        "lmm": argv[5],
        "o": argv[6],
        "n": argv[7]
    }

    # Run gemma
    run_gemma(gemma_params)

    # Parse log
    log_out = parse_log(os.path.join("output", gemma_params["o"]) + ".log.txt")
    print(
         f'Info: Finished {log_out["phenotype"]}\n'\
         f'Total individuals:\t{log_out["total_ind"]}\n'\
         f'Analyzed individuals:\t{log_out["analyzed_ind"]}\n'\
         f'Variants loaded:\t{log_out["total_var"]}\n'\
         f'Variants analyzed:\t{log_out["analyzed_var"]}\n'\
         )
    

if __name__ == "__main__":
    main()
