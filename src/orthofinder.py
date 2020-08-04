"""
Tools for converting OrthoFinder output to configuration files for I-ADHoRe

Input:
    - A subset of species from the OrthoFinder analysis
    - OrthoFinder output data frame
    - GFF files for each species
    - The GFF features to include for I-ADHoRe
    - The mappings of these features to genes in the GFF

*Note: The last two may differ for different species, so these should be
specified per species, using a multi-option thing in click.

Output = a directory containing:
    - Configuration file for I-ADHoRe
    - Families file for I-ADHoRe
    - Gene lists files for I-ADHoRe

Then some parameters can be changed in the config file and I-ADHoRe can be run
as usual. The eventual CLI could include a full wrapper for I-ADHoRe but that is
no priority.

*Note: GFF utils is crazy slow, and we mostly want only the genes, so a custom
GFF parser might be better after all.
"""
import pandas as pd
import numpy as np
import logging
import os


def get_families_orthofinder(orthofinder_df, species):
    """
    Get the gene families for the species of interest. This is the only
    interaction with the OrthoFinder dataframe I guess. Species identifiers
    should be unique prefixes.
    """
    df = pd.read_csv(orthofinder_df, sep="\t", index_col=0)
    # let's be forgiving when it comes to the species names, allow prefixes
    unique_prefix(df.columns, species)
    f = lambda x, y: x.startswith(y)
    cols = [x for x in df.columns if any([f(x, y) for y in species])]
    df = df[cols]
    df.fillna("", inplace=True)
    return df


def write_families_from_df(df, fname):
    """
    Write families file for I-ADHoRe based on a data frame where every row
    consists of one family.
    """
    genes = {}
    with open(fname, "w") as f:
        for row in df.index:
            for x in df.loc[row].apply(lambda x : x.split(", ")):
                gs = [y for y in x if y != ""]
                strs = ["{}\t{}\n".format(y, row) for y in gs]
                f.write("".join(strs))
                for g in gs:
                    genes[g] = row
    return os.path.abspath(fname), genes


def write_families_from_mcl(df, fname):
    genes = {}
    with open(fname, "w") as o:
        with open(df, "r") as f:
            for i, l in enumerate(f.readlines()):
                gs = l.strip().split()
                strs = ["{}\t{}\n".format(g, i) for g in gs]
                o.write("".join(strs))
                for g in gs:
                    genes[g] = row
    return os.path.abspath(fname), genes


def unique_prefix(strs, prefixes):
    f = lambda x, y: x.startswith(y)
    for y in prefixes:
        if len([x for x in strs if f(x, y) == True]) > 1:
            raise ValueError("Not a unique species prefix: " + y)
    return True


def orthogroup_poisson_filter(df, threshold=3):
    df.fillna("", inplace=True)
    n0 = len(df.index)
    counts = df[df.columns].applymap(
        lambda x: len([y for y in x.split(", ") if y != '']))
    family_sizes = counts.apply(lambda x: sum(x), axis=1)
    idx = poisson_outlier(family_sizes, threshold)
    counts = counts.loc[idx]
    df = df.loc[idx]
    n1 = len(df.index)
    idx = counts.apply(lambda x: np.all(poisson_outlier(x, threshold)), axis=1)
    df = df.loc[idx]
    n2 = len(df.index)
    logging.info("Retained {}/{} families after family-size outlier "
        "filter".format(n1, n0))
    logging.info("Retained {}/{} families after entry outlier filter".format(
        n2, n1))
    return df


def poisson_outlier(vector, threshold):
    """
    Returns a vector with True for non-outliers, False for outliers.
    """
    y = 2*np.sqrt(vector)
    return y < np.median(y) + threshold
