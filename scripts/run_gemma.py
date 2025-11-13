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


def main():
    if len(sys.argv) != 9:
            print("Usage: python run_gemma.py <bfiles> <kinship_mat> <maf> "
                  "<missing_threshold> <lmm> <output filename> <phenotype_num>")
            sys.exit(1)
    # Make parameter_dictionary
    gemma_params = {
        "bfile": sys.argv[1],
        "k": sys.argv[2],
        "maf": sys.argv[3],
        "miss": sys.argv[4],
        "lmm": sys.argv[5],
        "o": sys.argv[6],
        "n": sys.argv[7]
    }

    # Run gemma
    run_gemma(gemma_params)


if __name__ == "__main__":
    main()
