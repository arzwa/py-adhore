Copyright (C) 2020 Arthur Zwaenepoel

*This repository is not actively maintained, use at own risk*

# py-adhore

Utilities for doing and visualizing synteny/co-linearity analyses using
I-ADHoRe.  This should eventually replace the tools provided by `wgd syn` in the
[wgd]( https://github.com/arzwa/wgd/) package. It is mainly a wrapper to enable
easier usage of the very fast and sensitive [I-ADHoRe program](http://bioinformatics.psb.ugent.be/webtools/i-adhore/licensing/) for synteny and co-
linearity inference in large comparative genomic data sets; and provides the tools
for associated visualizations one usually desires.

To install, open a terminal session, clone the repository, `cd` into it and
type

```
$ pip install .
```

Works with python3 (you may need to use `pip3` instead of `pip` in the above
command; consider using a [virtualenv](https://virtualenv.pypa.io/en/stable/)).

I-ADHoRe can be downloaded [here](http://bioinformatics.psb.ugent.be/webtools/i-adhore/licensing/) (If the page prompts you to go to the updated website, decline
and go to the bottom of the page to download the package (v3.0)).

## Usage

### Running I-ADHoRe on OrthoFinder gene family data

The `example` directory of this repository contains an example data set,
consisting of gene families (`og.csv`) inferred with OrthoFinder for
*Arabidopsis thaliana* and *Vitis vinifera* and two `.gff` files for the
associated genes of these species (downloaded, along with the proteome fasta
files used in OrthoFinder, from [PLAZA
4.0](https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_dicots/download/index)).

To run I-ADHoRe for this data set you can use the following command

```
$ py-adhore of ./Orthogroups.tsv ath,vvi ./ath.gff ./vvi.gff -f gene -a Name -n 4 --run
```

This will use four threads (`-n 4`) for the I-ADHoRe algorithm. Say you only
wanted to do a within-genome comparison for *Vitis*, you could simply do

```
$ py-adhore of ./Orthogroups.tsv vvi ./vvi.gff -f gene -a Name -o vvi_out -n 4 --run
```

You can find more options and instructions using `py-adhore -of --help`.

To get the associated synteny clusters, use:

```
$ py-adhore cl py-adhore.out/genes_data.csv py-adhore.out/i-adhore-out/anchorpoints.txt
```

We make the following notes:

1. It is very important that the gene IDs in the OrthoFinder gene families
   correspond to the gene IDs found in the features of the gff file that
   correspond with the `-f` setting (here `gene` features) at the attribute
   specified with the `-a` option (here `Name`). If this is not clear have a
   look at the example data.
2. One can use any number and any subset of species from the gene families
   file.
3. One can use a plain MCL gene families file (where each line is whitespace
   separated list of genes) by using the `--mcl` flag. In this case, gene IDs
   must have species IDs as prefixes
4. One can choose to just generate the configuration files for I-ADHoRe without
   running it, simply omit the `--run` flag.
5. Note that if you're using OrthoFinder, you may want to add the singleton genes
   in `Orthogroups.UnassignedGenes.tsv` to the `Orthogroups.tsv` file before
   using with py-adhore, otherwise you will ignore all genes that are in 
   singleton gene families.

### Visualization

One can visualize the results in a circos diagram using Circos or Circos.js as
follows (example for circos.js):

```
$ py-adhore cc --js py-adhore.out/i-adhore-out/segments.txt ./py-adhore.out/genes_data.csv
```

This will create a directory `py-adhore.out/circos` with a HTMl file in the
case of Circos.js. If the `--js` flag is omitted configuration files for Circos
are generated.

## Citation

If you use this code, please do not forget to cite

```
Proost, Sebastian, et al. "i-ADHoRe 3.0—fast and sensitive detection of genomic homology in extremely large data sets."
Nucleic acids research 40.2 (2011): e11-e11.
```
