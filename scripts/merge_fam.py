#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 12-11-2025

Takes a .fam file and a phenotype file, merges phenotype into fam file. 
"""


def parse_phenotype_file(
    infile_fn: str, 
    delim: str="\t"
    ) -> tuple[list[str], dict[str, list[str]]]:
    """_summary_

    :param infile_fn: _description_
    :type infile_fn: str
    :param delim: _description_, defaults to "\t"
    :type delim: str, optional
    :raises Exception: _description_
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
    phenotypes: list
    ) -> str:
    n_pheno = len(phenotypes)
    out_lines = []
    temp_fam_fn = fam_fn.replace(".fam", "_temp.fam")

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
    
    with open(temp_fam_fn, "w") as outfile:
        for line in out_lines:
            outfile.write(" ".join(line) + "\n")
    return temp_fam_fn


def main():
    pheno_fn = r"..\test_files\phenotypes_galore_NA_filtered_pheno_mismatch.txt"
    fam_fn = r"..\test_files\LKsat_bed_files_het-0.1_v4.fam"
    
    phenotypes, pheno_dict = parse_phenotype_file(pheno_fn)
    temp_fam_fn = modify_fam(fam_fn, pheno_dict, phenotypes)


if __name__ == "__main__":
    main()
    