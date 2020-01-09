# Tools for getting Circos visualizations from I-ADHoRe output
import pandas as pd
import os
import matplotlib
import random
pd.set_option('mode.chained_assignment', None)


def write_karyotype_circos(gdata, fname):
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


def ribbons_to_json(seg, gdata, minlen=1e6):
    json_list = []
    d = seg.groupby("multiplicon")["first"].apply(list)
    e = seg.groupby("multiplicon")["last"].apply(list)
    df = pd.concat([d, e], axis=1)
    ch_colors = get_chord_colors(list(gdata["sp"].unique()))
    lines = []
    for mp in df.index:
        n = len(df.loc[mp]["first"])
        for i in range(n):
            s1 = gdata.loc[df.loc[mp]["first"][i]]["sp"]
            c1 = gdata.loc[df.loc[mp]["first"][i]]["chrom"]
            x1 = gdata.loc[df.loc[mp]["first"][i]]["start"]
            x2 = gdata.loc[df.loc[mp]["last"][i]]["stop"]
            if x2 - x1 < minlen:
                continue
            for j in range(i+1, n):
                s2 = gdata.loc[df.loc[mp]["first"][j]]["sp"]
                c2 = gdata.loc[df.loc[mp]["first"][j]]["chrom"]
                y1 = gdata.loc[df.loc[mp]["first"][j]]["start"]
                y2 = gdata.loc[df.loc[mp]["last"][j]]["stop"]
                d = {"source": {"id": c1, "start": x1, "end": x2},
                    "target": {"id": c2, "start": y1, "end": y2}}
                d["value"] = ch_colors["_".join(sorted([s1, s2]))]
                json_list.append(d)
    return json_list

def reduce_karyotype(kt, ri):
    chord_elements = set([x["source"]["id"] for x in ri])
    chord_elements = chord_elements | set([x["target"]["id"] for x in ri])
    kt = [x for x in kt if x["id"] in chord_elements]
    return kt


def karyotype_to_json(gdata, minlen=5000000):
    json_list = []
    df = gdata.groupby(["chrom"])[["stop", "sp"]].max()
    df = df[df["stop"] > minlen]
    df = df.sort_values(["sp", "stop"], ascending=False)
    colors = get_colors(list(df["sp"].unique()))
    chrcolors = {}
    for chr in df.index:
        chrcolors[chr] = colors[df.loc[chr]["sp"]]
        json_list.append({"id": chr, "len": df.loc[chr]["stop"],
            "label": df.loc[chr]["sp"], "color": colors[df.loc[chr]["sp"]]})
    return json_list


def nice_colors():
    return ['#5D8AA8', '#FE6F5E', '#177245', '#B22222', '#4F7942',
        '#86608E','#FCC200', '#ffe119', '#911eb4', '#46f0f0', '#f032e6',
        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',  '#808000',
        '#ffd8b1', '#000075', '#808080']


def get_colors(species):
    colors = nice_colors()
    return {sp: colors[i] for i, sp in enumerate(sorted(species))}


def get_chord_colors(species):
    cols = {}
    ns = len(species)
    species = sorted(species)
    cvals = nice_colors()
    for i, sp in enumerate(species):
        cols["{}_{}".format(sp, sp)] = i
    c = len(species)
    for i in range(ns):
        for j in range(i+1, ns):
            cols["_".join(sorted([species[i],species[j]]))] = c
            c += 1
    return cols


def get_circosjs_doc(karyotype, ribbons, fname, pars):
    info = "Based on {0} and {1}".format(pars[-2], pars[-1])
    if pars[3]:
        info += ", showing all elements > {0} and all blocks > {1}".format(
            pars[0], pars[1])
        info += " (all = True). "
    else:
        info += ", showing all elements with blocks > {0} and of total length > {1} all chords > {1}".format(pars[1], pars[0])
        info += " (all = False). "
    info += "Only multiplicons of order > {} are shown.".format(pars[2]-1)
    with open(fname, "w") as f:
        f.write(get_html(get_circosjs(karyotype, ribbons), info))


def get_circosjs(karyotype, ribbons):
    js = """
    var thecolors= ['#5D8AA8', '#FE6F5E', '#177245', '#B22222', '#4F7942',
        '#86608E','#FCC200', '#ffe119', '#911eb4', '#46f0f0', '#f032e6',
        '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',  '#808000',
        '#ffd8b1', '#000075', '#808080']

    var myCircos = new Circos({
        container: '#chart',
        width: 1500,
        height: 1500,
    });

    var configuration = {
        tooltipContent: function (d) {
            return d.label
        },
        innerRadius: 500,
        outerRadius: 520,
        cornerRadius: 0,
        gap: 0.01, // in radian
        labels: {
            display: false,
            position: 'center',
            size: '14',
            color: '#000000',
            radialOffset: 80,
        },
        ticks: {
            display: false,
            color: 'grey',
            spacing: 10000000,
            labels: true,
            labelSpacing: 10,
            labelSuffix: 'Mb',
            labelDenominator: 1000000,
            labelDisplay0: true,
            labelSize: '10px',
            labelColor: '#000000',
            labelFont: 'default',
            majorSpacing: 5,
            size: {
                minor: 2,
                major: 5,
            }
        },
        events: {}
    }
    """
    colorconf = "{'color': function(datum, index) { return thecolors[datum.value] }, opacity: 0.5}"
    #colorconf = "{'color:' 'Spectral', 'opacity': 0.75}"
    #chord_conf = {'color': colorfun, 'opacity': 0.75}
    #chord_conf = {'color': 'Spectral', 'opacity': 0.75}
    js += "\nvar layout_data = " + str(karyotype)
    js += "\nmyCircos.layout(layout_data, configuration);"
    js += "\nmyCircos.chords('synteny', {}, {});".format(
        str(ribbons), str(colorconf))
    js += "\nmyCircos.render();"
    return js


def get_html(js_str, info):
    html = "<!DOCTYPE html><html><head><script src='https://cdn.rawgit.com/"
    html += "nicgirault/circosJS/v2/dist/circos.js'></script>"
    html += """<link rel="stylesheet" href="https://fonts.googleapis.com/icon?family=Material+Icons"> <link rel="stylesheet" href="https://code.getmdl.io/1.3.0/material.indigo-pink.min.css"> <script defer src="https://code.getmdl.io/1.3.0/material.min.js"></script>
    </head><body>"""
    html += "<div><p>{}</p></div>".format(info)
    html += "<div width='100%'><svg width='1500' height='1500' "
    html += "id='chart'></svg></div>"
    html += "<script>{}</script></body></html>".format(js_str)
    return html
