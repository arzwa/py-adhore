"""
--------------------------------------------------------------------------------
Copyright (C) 2018 Arthur Zwaenepoel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact: arzwa@psb.vib-ugent.be
--------------------------------------------------------------------------------
"""
import click
import logging
import os
import sys
import coloredlogs
import subprocess
from src.orthofinder import *
from src.gffparser import *
from src.utils import *
from src.adhore import *
from src.circos import *


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def cli():
    coloredlogs.install(fmt='%(asctime)s: %(levelname)s\t%(message)s',
        level="INFO", stream=sys.stdout)
    pass


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('data_frame', nargs=1, type=click.Path(exists=True))
@click.argument('species', nargs=1)
@click.argument('gff', nargs=-1, type=click.Path(exists=True))
@click.option("--run", is_flag=True)
@click.option("--mcl", is_flag=True)
@click.option('--features', '-f', default="gene", show_default=True,
    type=str, help='Features to use from gff files for each species,'
    'multiple comma-separated features per species allowed.')
@click.option('--attributes', '-a', default="ID", show_default=True,
    type=str, help='Attributes to use from gff files for each species,'
    'multiple comma-separated attributes per species allowed.')
@click.option('--outlier_filter', '-of', default=-1, show_default=True,
    help='Poisson outlier filtering threshold (< 0: no filtering)')
@click.option('--outdir', '-o', default='py-adhore.out', show_default=True,
    help='output directory')
@click.option('--gap_size', default=30, show_default=True)
@click.option('--cluster_gap', default=35, show_default=True)
@click.option('--q_value', default=0.75, show_default=True)
@click.option('--prob_cutoff', default=0.01, show_default=True)
@click.option('--anchor_points', default=3, show_default=True)
@click.option('--level_2_only', default="false", show_default=True)
@click.option('--alignment_method', default="gg2", show_default=True)
@click.option('--number_of_threads', '-n', default=1, show_default=True)
def of(data_frame, species, run, mcl, gff, features, attributes, outlier_filter,
        outdir, **kwargs):
    """
    Orthofinder/MCL to I-ADHoRE 3.0
    """
    import pandas as pd
    species = species.split(",")
    if len(species) != len(gff):
        logging.error("# of species is different from # of gff files")
        exit()
    try:
        os.mkdir(outdir)
    except FileExistsError:
        logging.warning("Output directory `{}` already exists".format(outdir))
    conf = default_adhore_conf()

    # get features/attributes settings
    feat, attr = parse_feat_attr(features, attributes, species)

    # write families
    logging.info("Writing families file")
    fn = os.path.join(outdir, "families.tsv")
    if not mcl:
        df = get_families_orthofinder(data_frame, species)
        if outlier_filter > 0:
            logging.info("Filtering Poisson outlier families (> {})".format(
                outlier_filter))
            df = orthogroup_poisson_filter(df, threshold=outlier_filter)
        fn, genes = write_families_from_df(df, fn)
    else:
        fn, genes = write_families_from_mcl(data_frame, fn)
    conf["blast_table"] = fn

    # write gene lists
    logging.info("Writing gene lists")
    lconf, gdata = gffs_to_genelists(gff, feat, attr, genes, outdir)
    gdfname = os.path.join(outdir, "genes_data.csv")
    if len(gdata.index) == 0:
        logging.error("No genes found, make sure the --feature and "
            "--attribute specify the correspondence between gff and families")
        exit()
    gdata.to_csv(gdfname)
    conf["lists"] = lconf

    # write karyotype
    logging.info("Writing gene-based karyotype")
    write_karyotype(gdata, os.path.join(outdir, "karyotype.csv"))

    # write configuration file
    logging.info("Writing I-ADHoRe 3.0 configuration file")
    conf.update(kwargs)
    conf["output_path"] = os.path.abspath(os.path.join(outdir, "i-adhore-out"))
    write_adhore_config(conf, os.path.join(outdir, "adhore.conf"))

    if run:
        command = ["i-adhore", os.path.join(outdir, "adhore.conf")]
        logging.info("Running I-ADHoRe [`{}`]".format(" ".join(command)))
        subprocess.run(command)#, capture_output=True)
        summarize_adhore(conf["output_path"], gdfname,
            os.path.join(outdir, "py-adhore.csv"))
        logging.info("These were the parameters for I-ADHoRe: ")
        for k, v in kwargs.items():
            print("--{} '{}' ".format(k, v), end="")
        print()


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('segments', nargs=1, type=click.Path(exists=True))
@click.argument('genesdata', nargs=1)
@click.option("--js", is_flag=True)
@click.option("--all", is_flag=True)
@click.option('--minlen_kt', '-mk', default=5e6, help='Minimum length of '
    'genomic element')
@click.option('--minlen_ch', '-mc', default=1e6, help='Minimum length of '
    'a syntenic block')
@click.option('--min_order', '-mo', default=2, help='Minimum order of '
    'a multiplicon')
@click.option('--outdir', '-o', default=None, help='output directory')
def cc(segments, genesdata, all, js, minlen_kt, minlen_ch, min_order, outdir):
    """
    Circos visualization of I-ADHoRe results
    """
    import pandas as pd
    if not outdir:
        outdir = os.path.join(os.path.dirname(genesdata), "circos")
    try:
        os.mkdir(outdir)
    except FileExistsError:
        logging.warning("Output directory `{}` already exists".format(outdir))
    seg = pd.read_csv(segments, sep="\t", index_col=0)
    seg = segments_filter(seg, min_order)
    gdata = pd.read_csv(genesdata, index_col=0)
    if not js:
        logging.warning("Using --js is recommended")
        kt = write_karyotype(gdata, os.path.join(outdir, "karyotype.txt"))
        ri = write_ribbons(seg, gdata, os.path.join(outdir, "ribbons.txt"))
        cc = write_circos_conf(kt, ri, os.path.join(outdir, "circos.conf"))
    else:
        kt = karyotype_to_json(gdata, minlen=minlen_kt)
        ri = ribbons_to_json(seg, gdata, minlen=minlen_ch)
        if not all: kt = reduce_karyotype(kt, ri)
        get_circosjs_doc(kt, ri, os.path.join(outdir, "circos.html"),
            [minlen_ch, minlen_kt, min_order, all, segments, genesdata])


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('genesdata', nargs=1, type=click.Path(exists=True))
@click.argument('anchorpoints', nargs=1, type=click.Path(exists=True))
@click.option('--outdir', '-o', default="clusters",
              type=click.Path(exists=False), show_default=True)
def cl(genesdata, anchorpoints, outdir):
    """
    Get syntenic clusters from the syntenic network.

    Syntenic clusters are defined as connected components in the
    synteny networs. These are actually the subset of the orthogroups
    that are anchor pairs.
    """
    import pandas as pd
    import src.network
    genesdata = pd.read_csv(genesdata, index_col=0)
    genes, counts = src.network.get_clusters(anchorpoints, genesdata)
    os.mkdir(outdir)
    genes.to_csv(os.path.join(outdir, "clusters.tsv"), sep="\t")
    counts.to_csv(os.path.join(outdir, "profile.csv"), sep=",")

if __name__ == '__main__':
    cli()
