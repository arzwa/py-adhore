import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
from matplotlib.colors import Normalize


def karyotype_axes(ax, kt, g1, g2, minlen=1000000, labelsize=9):
    ax.tick_params(axis='both', which='major', labelsize=labelsize)
    kt = kt[kt["stop"] > minlen]
    ktsp1 = kt[kt["sp"] == g1]
    ktsp2 = kt[kt["sp"] == g2]
    xticks = np.cumsum(ktsp1["stop"].sort_values(ascending=False))
    yticks = np.cumsum(ktsp2["stop"].sort_values(ascending=False))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks.index, rotation=45, ha="right")
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks.index, rotation=45, ha="right")
    ax.grid(linestyle=":")
    ax.set_xlim(0, xticks[-1])
    ax.set_ylim(0, yticks[-1])
    return xticks, yticks


def get_segments(df, g1, g2, xticks, yticks):
    index = {g1: 0, g2:1}
    ticks = [xticks, yticks]
    # g1 corrsponds to x-axis genome, g2 to the y-axis genome
    df = df[df["genome"].isin([g1, g2])]
    df = df[df["list"].isin(list(xticks.index) + list(yticks.index))]
    df = df.groupby("multiplicon").agg(list)
    df = df[df["genome"].apply(set) == set([g1, g2])]
    df = df[df["genome"].apply(lambda x: len(x) > 1)]
    # now we have filtered the data frame down to all we need
    # now it;'s probably easiest to just loop over the segments.
    segs = []
    #return df
    for r in df.index:
        mcon = df.loc[r]
        nseg = len(mcon["genome"])
        for i in range(nseg):
            gi = mcon["genome"][i]
            si = mcon["list"][i]
            indi = index[gi]
            for j in range(i+1,nseg):
                gj = mcon["genome"][j]
                sj = mcon["list"][j]
                x = [[0,0],[0,0]]
                if g1 == g2:
                    indj = 0
                else:
                    indj = index[gj]
                x[indi][0] = ticks[index[gi]][si] + mcon["start"][i]
                x[indi][1] = ticks[index[gi]][si] + mcon["stop"][i]
                x[indj][0] = ticks[index[gj]][sj] + mcon["start"][j]
                x[indj][1] = ticks[index[gj]][sj] + mcon["stop"][j]
                segs.append((x, r))
    return segs


def plot_segments(ax, segs, linewidth=2, color="k"):
    for x, mcon in segs:
        ax.plot(x[0], x[1], color=color, linewidth=linewidth)
        ax.plot(x[1], x[0], color=color, linewidth=linewidth)


def dotplot(ax, kt, df, g1, g2=None, minlen=500000, color="k", linewidth=2,
        labelsize=8):
    if not g2: g2 = g1
    xticks, yticks = karyotype_axes(ax, kt, g1, g2, labelsize=labelsize,
        minlen=minlen);
    segs = get_segments(df, g1, g2, xticks, yticks);
    plot_segments(ax, segs,color=color, linewidth=linewidth)
    return ax


def get_anchor_ks(anchors, ks_distribution):
    an["pid"] = an[['gene_x', 'gene_y']].apply(
        lambda x: '__'.join(sorted(x)), axis=1)
    return an.join(ks, on="pid")


def get_median_ks(df, anchors, ks_distribution):
    anks = get_anchor_ks(anchors, ks_distribution)
    df = df.join(get_anchor_ks(an, ks).groupby(
        "multiplicon")["Ks"].median(), on="multiplicon")
    return df


def colorhack(cmap):
    z = [[0, 0], [0, 0]]
    levels = range(0, 101, 1)
    tmp = plt.contourf(z, levels, cmap=cmap)
    return tmp


def dotplot_ks(ax, kt, df, g1, g2=None, minlen=500000, color="k", linewidth=2,
        labelsize=8, cmap=cm.viridis, outlier_color="firebrick", max_ks=5):
    tmp = colorhack(cmap)
    if not g2: g2 = g1
    xticks, yticks = karyotype_axes(
        ax, kt, g1, g2, labelsize=labelsize, minlen=minlen);
    segs = get_segments(df, g1, g2, xticks, yticks);
    plot_segments_ks(ax, segs, df, cmap=cmap, linewidth=linewidth,
        max_ks=max_ks, outlier_color=outlier_color)
    cbar = plt.colorbar(tmp, fraction=0.02, pad=0.01)
    cbar.ax.set_yticklabels(
        ['{:.2f}'.format(x) for x in np.linspace(0.0, max_ks, 11)])
    return ax


def plot_segments_ks(ax, segs, df, linewidth=2, color="k", max_ks=5,
        cmap=cm.viridis, outlier_color="firebrick"):
    norm = Normalize(vmin=0, vmax=max_ks)
    for x, mcon in segs:
        ks = df[df["multiplicon"] == mcon].iloc[0]["Ks"]
        if ks > max_ks:
            color = outlier_color
        else:
            color = cmap(norm(ks))
        ax.plot(x[0], x[1], color=color, linewidth=linewidth)
        ax.plot(x[1], x[0], color=color, linewidth=linewidth)
    return ax
