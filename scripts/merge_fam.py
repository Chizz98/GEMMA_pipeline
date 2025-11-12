#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 12-11-2025

Takes a .fam file and a phenotype file, merges phenotype into fam file. 
"""

def parse_phenotype_file(
    infile_fn: str, 
    delim: str="\t"
    ) -> tuple[list[str], dict[list[str]]]:
    """_summary_

    :param infile_fn: _description_
    :type infile_fn: str
    :param delim: _description_, defaults to "\t"
    :type delim: str, optional
    :raises Exception: _description_
    :return: list containting phenotypes, 
        dictionary with sample ids as keys and phenotypes as values
    :rtype: tuple[list[str], dict[list[str]]]
    """
    out_dict = {}

    with open(infile_fn, "r") as infile:
        header = infile.readline().strip().split(delim)
        pheno_header = header[1:]
        
        for line in infile:
            line = line.strip().split(delim)
            sample_id = line[0]
            phenotypes = line[1:]

            if sample_id in out_dict:
                raise Exception(
                    f"Sample ID {sample_id} exists multiple times in phenotype "
                    f"file"
                )
            out_dict[sample_id] = phenotypes
    return pheno_header, out_dict



def main():
    pheno_fn = r"..\test_files\phenotypes_galore_NA_filtered.txt"
    
    phenotypes, pheno_dict = parse_phenotype_file(pheno_fn)
    print(phenotypes)


if __name__ == "__main__":
    main()
    