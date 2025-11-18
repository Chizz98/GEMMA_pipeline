#!/usr/bin/env python3
"""
Author: Chris Dijkstra
Date: 12-11-2025

Takes a .fam file and a phenotype file, merges phenotype into fam file. 
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy


def calc_inflation(p_vals: list, df: int) -> float:
    p_vals = np.array(p_vals)
    p_vals = p_vals[~np.isnan(p_vals)]
    chi2_stats = scipy.stats.chi2.isf(p_vals, df=df)
    expected_median = scipy.stats.chi2.ppf(0.5, df=df)
    genomic_inflation = np.median(chi2_stats) / expected_median
    return genomic_inflation


def parse_assoc(assoc_fn:str) -> dict:
    out_dict = {
        "chrom_list": [],
        "pos_list": [],
        "p_vals": [],
        "actual_pos": []
    }

    # Parse assoc
    with open(assoc_fn) as infile:
        header = infile.readline().strip().split()
        
        chr_max = 0
        cumulative_pos = 0
        last_chr = ""
        for line in infile:
            line = line.strip().split()
            chrom = line[0]
            pos = int(line[2])
            p_val = float(line[11])

            if chrom != last_chr:
                cumulative_pos += chr_max
                chr_max = 0
                last_chr = chrom
            
            pos_cum = pos + cumulative_pos

            if pos > chr_max:
                chr_max = pos
            
            out_dict["chrom_list"].append(chrom)
            out_dict["pos_list"].append(pos_cum)
            out_dict["p_vals"].append(p_val)
            out_dict["actual_pos"].append(pos)
    return out_dict


def manhattan_plot(
        assoc_dict: dict, 
        out_fn: str, 
        lod_plotting_thresh: float
        ) -> None:
    chrom_list = assoc_dict["chrom_list"]
    pos_list = assoc_dict["pos_list"]
    p_vals = assoc_dict["p_vals"]
    actual_pos = np.array(assoc_dict["actual_pos"])
    
    pos_list = np.array(pos_list)
    LOD = -np.log10(np.array(p_vals))

    plotting_chr = np.array(chrom_list)[LOD > lod_plotting_thresh]
    unique_chr = np.unique(plotting_chr)
    plotting_pos = pos_list[LOD > lod_plotting_thresh]
    plotting_lod = LOD[LOD > lod_plotting_thresh]
    plotting_actual_pos = actual_pos[LOD > lod_plotting_thresh]
    bonf_thresh = -np.log10(0.05 / len(p_vals))
    
    colors = ["#1f77b4", "#ff7f0e"]
    xticks = []
    xlabs = []

    plt.figure(figsize=(12, 5))

    for i, contig in enumerate(unique_chr):
        pos_subset = plotting_pos[plotting_chr == contig]
        lod_subset = plotting_lod[plotting_chr == contig]
        actual_pos_subset = plotting_actual_pos[plotting_chr == contig]
        plt.scatter(
            x = pos_subset, 
            y = lod_subset, 
            color = colors[i % 2],
            s=6
            )
        midpoint = (pos_subset.min() + pos_subset.max()) / 2
        xticks.append(midpoint)
        xlabs.append(contig)
        
        # Plot highest lod variant per chr
        max_ind = np.argmax(lod_subset)
        top_pos = pos_subset[max_ind]
        top_actual_pos = actual_pos_subset[max_ind]
        top_snp_lod = lod_subset[max_ind]
        plt.text(top_pos, top_snp_lod + 0.05, str(top_actual_pos), fontsize=6, 
                 rotation=45)

    
    plt.axhline(y=bonf_thresh, color="red", linestyle="dashed")
    plt.xticks(xticks, xlabs, rotation=90)
    plt.tight_layout()
    plt.savefig(out_fn, dpi=300)
    plt.close()


def qqplot(
        assoc_dict: str, 
        out_fn:str, 
        df:int = 1, 
        ) -> None:
    p_vals = np.array(assoc_dict["p_vals"])
    p_vals = p_vals[~np.isnan(p_vals)]

    observed = -np.log10(np.sort(p_vals))
    n = len(observed)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)

    inflation_factor = calc_inflation(p_vals, df=df)

    median_idx = n // 2
    median_expected = expected[median_idx]
    median_observed = observed[median_idx]

    plt.figure(figsize=(5, 5))
    plt.scatter(expected, observed, s=6)

    # Plot refline from 0 to max on both axes
    maxval = max(expected.max(), observed.max())
    plt.plot([0, maxval], [0, maxval], color="red", linestyle="--",
             label="reference line")

    # Plot median lines for genomic inflation factor
    plt.axhline(y=median_observed, color="cyan", linestyle=":", 
                label="median observed")
    plt.axvline(x=median_expected, color="magenta", linestyle=":",
                label="median expected")

    plt.legend(loc="lower right")
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.text(0.05*maxval, 0.9*maxval, f"λGC = {inflation_factor:.3f}", 
             color="black")

    plt.tight_layout()
    plt.savefig(out_fn, dpi=300)
    plt.close()

                     
def main(argv=None):
    argv = argv or sys.argv
    if len(argv) != 4:
            print("Usage: python make_manhattan.py <assoc_fn> <out_prefix> "
            "<lod_plotting_th>")
            sys.exit(1)
    
    assoc_fn = argv[1]
    out_fn = argv[2]
    lod_th = argv[3]

    assoc_dict = parse_assoc(assoc_fn)
    
    manhattan_fn = out_fn + "_manhattan.png"
    manhattan_plot(assoc_dict, manhattan_fn, lod_th)

    qq_fn = out_fn + "_qqplot.png"
    qqplot(assoc_dict, qq_fn)


if __name__ == "__main__":
    main()
        
