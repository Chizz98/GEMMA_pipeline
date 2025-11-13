#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 12-11-2025

Takes a .fam file and a phenotype file, merges phenotype into fam file. 
"""

import shutil
import sys
import pathlib


def parse_phenotype_file(
    infile_fn: str, 
    delim: str="\t"
    ) -> tuple[list[str], dict[str, list[str]]]:
    """Parses phenonotype file into dict and list of phenotypes

    :param infile_fn: Input phenotype file, with sample id as first column, 
        phenotypes as all other cols. Header should contain phenotype names.
    :type infile_fn: str
    :param delim: Delimiter of phenotype file, defaults to "\t"
    :type delim: str, optional
    :return: list containting phenotypes, 
        dictionary with sample ids as keys and phenotypes as values
    :rtype: tuple[list[str], dict[str, list[str]]]
    """
    out_dict = {}

    with open(infile_fn, "r") as infile:
        header = infile.readline().strip().split(delim)
        pheno_header = header[1:]
        
        for line in infile:
            line = line.strip().split(delim)
            sample_id = line[0]
            phenotypes = line[1:]

            for i, phenotype in enumerate(phenotypes):
                try:
                    float(phenotype)
                except ValueError:
                    if "," in phenotype:
                        raise Exception(
                            f"Comma in: {phenotype}, for " 
                            f"sample {sample_id} and phenotype "
                            f"{pheno_header[i]}, use periods as decimal symbol"
                            f" instead"
                        )
                    else:
                        raise Exception(
                            f"Non numeric phenotype detected: {phenotype}, for " 
                            f"sample {sample_id} and phenotype "
                            f"{pheno_header[i]}"
                            )

            if sample_id in out_dict:
                raise Exception(
                    f"Sample ID {sample_id} exists multiple times in phenotype "
                    f"file"
                )

            if len(phenotypes) != len(pheno_header):
                raise Exception(
                    f"Mismatch: {len(phenotypes)} phenotypes counted "
                    f"for {sample_id}, {len(pheno_header)} phenotypes in header"
                )
            
            out_dict[sample_id] = phenotypes
    return pheno_header, out_dict


def modify_fam(
    fam_fn: str,
    pheno_dict: dict,
    phenotypes: list,
    out_fn: str
    ) -> None:
    """Merges phenotypes and an input fam file into out fam

    :param fam_fn: Filename of input fam
    :type fam_fn: str
    :param pheno_dict: Pheno dict output by parse_phenotype_file
    :type pheno_dict: dict
    :param phenotypes: List of phenotypes
    :type phenotypes: list
    :param out_fn: Out filename
    :type out_fn: str
    :return: None, writes merged fam to specified outfile
    :rtype: None
    """
    n_pheno = len(phenotypes)
    out_lines = []

    with open(fam_fn) as infile:
        for line in infile:
            line = line.strip().split()
            sample_descriptors = line[:5]
            sample_id = sample_descriptors[1]
            
            # -9 is missing phenotype symbol for GEMMA
            sample_phenotypes = pheno_dict.get(sample_id, ["-9"] * n_pheno)
            if sample_id not in pheno_dict:
                print(f"Warning: no phenotypes for {sample_id}")

            out_lines.append(sample_descriptors + sample_phenotypes)
    
    with open(out_fn, "w") as outfile:
        for line in out_lines:
            outfile.write(" ".join(line) + "\n")


def main(argv=None):
    argv = argv or sys.argv
    # Parse cmd args
    if len(argv) != 3:
            print("Usage: python merge_pheno.py <phenotype_file> <fam_file>")
            sys.exit(1)

    pheno_fn = pathlib.Path(argv[1])
    fam_fn = pathlib.Path(argv[2])

    # Construct temporary fam filename
    out_fam_fn = fam_fn.with_name(fam_fn.stem + "_temp" + fam_fn.suffix)
    
    phenotypes, pheno_dict = parse_phenotype_file(pheno_fn)
    modify_fam(fam_fn, pheno_dict, phenotypes, out_fam_fn)

    # Overwrite input fam with temp fam
    shutil.move(out_fam_fn, fam_fn)

    return phenotypes


if __name__ == "__main__":
    main()
    