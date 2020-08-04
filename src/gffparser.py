"""
The fully feature `gffutils` package is too slow on typical full GFF files for
the purposes here (database should be created on the fly and destroyed
afterwards).

A command line could look like

--gff file1 file2 --features gene,mRNA gene

for file1 with features gene & mRNA and file2 with features gene (idem for attr)
"""
import gffutils
import tempfile
import logging
import os
import pandas as pd
from collections import defaultdict


def gffs_to_genelists(fnames, features, attributes, genes, outdir):
    confs = []
    genes_data = {}
    for i, fname in enumerate(fnames):
        logging.info("Loading {} into database ... ".format(fname))
        feat = features[i]
        attr = attributes[i]
        db = load_gff(fname, feat)
        genome = os.path.basename(fname) + "_lists"
        gconf, gdata = write_gene_lists(db, genes, genome, outdir,
            attr=attr, features=feat)
        confs.append(gconf)
        genes_data.update(gdata)
    return confs, pd.DataFrame.from_dict(genes_data).transpose()


def load_gff(fn, features=["gene"]):
    """
    Load a gff file in a database, ignoring all features not in `features`.
    Uses tempfiles
    """
    with open(fn, "r") as f:
        lines = [l for l in f.readlines() if
            (l[0] != "#" and l.split("\t")[2] in features)]
    gff = tempfile.NamedTemporaryFile(suffix=".gff", mode="w+t")
    gff.write("\n".join(lines))
    gff.seek(0)
    dbfn = tempfile.NamedTemporaryFile(suffix=".db")
    db = gffutils.create_db(gff.name, dbfn=dbfn.name, force=True,
        merge_strategy='merge')
    # db is an open database connection, as soon as the connection is closed
    # the tempfile should be destroyed.
    return db


def write_gene_lists(db, genes, genome, path, attr="ID", features=["gene"]):
    """
    Write the gene lists for a single gff (single species). This is how the
    results should be organized (and this is also how it's written in the
    configuration file):
        ```
        genome=genome
        chromname path/first_chrom.lst
        chromname path/second_chrom.lst
        chromname path/third_chrom.lst
        ```
    """
    genome_path = os.path.join(path, genome)
    config = {"genome": genome, "lists": []}
    try:
        os.mkdir(genome_path)
    except FileExistsError:
        logging.warning("Directory `{}` already exists, will possibly "
            "overwrite".format(genome_path))
    genomic_elements = defaultdict(list)
    genes_data = {}
    c1, c2 = 0, 0
    for f in features:
        for g in db.features_of_type(f, order_by='start'):
            g_id = g.attributes[attr][0]
            if g_id not in genes.keys():
                logging.debug("Feature {} not found in families".format(g_id))
                c2 += 1
                continue
            c1 += 1
            genes_data[g_id] = {"family": genes[g_id], "feat": f, "chrom": g.chrom,
                                "sp": genome, "strand": g.strand, "start": g.start,
                                "stop": g.stop}
            genomic_elements[g.chrom].append(g_id + g.strand)
    c = c1 + c2
    logging.info("{:.2f}% of genes not found in families".format(100*c2/c))
    for chr, lst in genomic_elements.items():
        p = os.path.join(genome_path, "{}.lst".format(chr))
        with open(p, "w") as f:
            f.write("\n".join(lst))
        config["lists"].append("{} {}".format(chr, os.path.abspath(p)))
    return config, genes_data


def write_karyotype(gdata, fname):
    df = gdata.groupby(["chrom"])[["stop", "sp"]].max()
    df["start"] = 0
    df.to_csv(fname, sep=",")
    return os.path.abspath(fname)
