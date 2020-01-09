Copyright (C) 2020 Arthur Zwaenepoel

VIB/UGent center for plant systems biology - [Bioinformatics & evolutionary
genomics group](http://bioinformatics.psb.ugent.be/beg/)

# py-adhore

Utilities for doing and visualizing synteny/co-linearity analyses using
I-ADHoRe.  This is meant to replace the tools provided by `wgd syn` in the
[wgd]( https://github.com/arzwa/wgd/) package. It is mainly a wrapper to enable
easier usage of I-ADHoRe and to provide the associated visualizations one
usually desires.

To install, clone the repository, `cd` into it and type `pip install .`. Works
with python3.

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
$ py-adhore of ./og.csv ath,vvi ./ath.gff ./vvi.gff -f gene -a Name --run 
```

Say you only wanted to do a within-genome comparison for *Vitis*, you could
simply do

``` 
$ py-adhore of ./og.csv vvi ./vvi.gff -f gene -a Name -o vvi_out --run 
```

You can find more options and instructions using `py-adhore -of --help`.

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

If you use this package, please cite this repository and most importantly:

``` 
Proost, Sebastian, et al. "i-ADHoRe 3.0—fast and sensitive detection of
genomic homology in extremely large data sets." Nucleic acids research 40.2
(2011): e11-e11. 
```
