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
from src.orthofinder import *
from src.gffparser import *
from src.utils import *
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
@click.option('--features', '-f', default="gene", show_default=True,
    type=str, help='Features to use from gff files for each species,'
    'multiple comma-separated features per species allowed.')
@click.option('--attributes', '-a', default="ID", show_default=True,
    type=str, help='Attributes to use from gff files for each species,'
    'multiple comma-separated attributes per species allowed.')
@click.option('--outdir', '-o', default='py-adhore.out', show_default=True,
    help='output directory')
def of(data_frame, species, gff, features, attributes, outdir):
    """
    Orthofinder to I-ADHoRE 3.0
    """
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
    logging.info("Writing Families file")
    df = get_families_orthofinder(data_frame, species)
    fn = os.path.join(outdir, "families.tsv")
    fn, genes = write_families_from_df(df, fn)
    conf["blast_table"] = fn

    # write gene lists
    logging.info("Writing gene lists")
    lconf, gdata = gffs_to_genelists(gff, feat, attr, genes, outdir)
    gdata.to_csv(os.path.join(outdir, "genes_data.csv"))
    conf["lists"] = lconf

    # write configuration file
    logging.info("Writing I-ADHoRe 3.0 configuration file")
    conf["output_path"] = os.path.abspath(os.path.join(outdir, "i-adhore-out"))
    write_adhore_config(conf, os.path.join(outdir, "adhore.conf"))


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('segments', nargs=1, type=click.Path(exists=True))
@click.argument('genesdata', nargs=1)
@click.option('--outdir', '-o', default=None, help='output directory')
def cc(segments, genesdata, outdir):
    """
    Circos visualization of I-ADHoRe results
    """
    if not outdir:
        outdir = os.path.join(os.path.dirname(genesdata), "circos")
    try:
        os.mkdir(outdir)
    except FileExistsError:
        logging.warning("Output directory `{}` already exists".format(outdir))
    seg = pd.read_csv(segments, sep="\t", index_col=0)
    gdata = pd.read_csv(genesdata, index_col=0)
    kt = write_karyotype(gdata, os.path.join(outdir, "karyotype.txt"))
    ri = write_ribbons(seg, gdata, os.path.join(outdir, "ribbons.txt"))
    cc = write_circos_conf(kt, ri, os.path.join(outdir, "circos.conf"))


if __name__ == '__main__':
    cli()