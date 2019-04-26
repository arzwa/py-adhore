import pandas as pd
import os


# Parse file sproduced by I-ADHoRe
def summarize_adhore(outdir, gdata, fname):
    """
    Make a one file summary of the I-ADHoRe results. One syntenic block per row
    with a multiplicon ID, anchor points, start coordinate and stop coordinate.
    (Maybe also coordinate for each anchor within the segment. Also each anchor
    point should have an id corresponding to its anchor pairs.)
    """
    seg = pd.read_csv(
        os.path.join(outdir, "segments.txt"), sep="\t", index_col=0)
    gd = pd.read_csv(gdata, index_col=0)
    seg["start"] = list(gd.loc[seg["first"]]["start"])
    seg["stop"] = list(gd.loc[seg["last"]]["stop"])
    seg.to_csv(fname)
