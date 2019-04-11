# Tools for getting Circos visualizations from I-ADHoRe output
import pandas as pd
import os
pd.set_option('mode.chained_assignment', None)


def write_karyotype(gdata, fname):
    colors = ["green", "blue", "orange"]
    df = gdata.groupby(["chrom"])[["stop", "sp"]].max()
    df["start"] = 0
    df["chrom"] = df["label"] = df.index
    df["color"] = colors[0]
    for (i, sp) in enumerate(df["sp"].unique()):
        df["color"][df["sp"] == sp] = colors[i]
    df["chr"] = "chr"
    df["dash"] = "-"
    colorder = ["chr", "dash", "chrom", "label", "start", "stop", "color"]
    df[colorder].to_csv(fname, sep=" ", index=None, header=None)
    return os.path.abspath(fname)


def write_ribbons(seg, gdata, fname):
    d = seg.groupby("multiplicon")["first"].apply(list)
    e = seg.groupby("multiplicon")["last"].apply(list)
    df = pd.concat([d, e], axis=1)
    lines = []
    for mp in df.index:
        n = len(df.loc[mp]["first"])
        for i in range(n):
            gdata.loc[df.loc[mp]["first"][i]]
            c1 = gdata.loc[df.loc[mp]["first"][i]]["chrom"]
            x1 = gdata.loc[df.loc[mp]["first"][i]]["start"]
            x2 = gdata.loc[df.loc[mp]["last"][i]]["stop"]
            for j in range(i, n):
                c2 = gdata.loc[df.loc[mp]["first"][j]]["chrom"]
                y1 = gdata.loc[df.loc[mp]["first"][j]]["start"]
                y2 = gdata.loc[df.loc[mp]["last"][j]]["stop"]
                lines.append(" ".join([str(x) for x in
                    [c1, x1, x2, c2, y1, y2]]))
    with open(fname, "w") as f:
        f.write("\n".join(lines))
    return os.path.abspath(fname)


def write_circos_conf(karyotype, ribbons, fname):
    conf_str = """
    <<include etc/colors_fonts_patterns.conf>>
    <<include ideogram.conf>>
    <image>
    <<include etc/image.conf>>
    </image>
    chromosomes_units           = 1000000
    chromosomes_display_default = yes
    karyotype = {}
    <links>
    <link>
    ribbon = yes
    file = {}
    bezier_radius = 0.1r
    thickness = 1
    </link>
    </links>
    <<include etc/housekeeping.conf>>
    """.format(karyotype, ribbons)
    with open(fname, "w") as f:
        f.write(conf_str)
    return os.path.abspath(fname)
