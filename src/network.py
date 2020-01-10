# the anchorpoints output from I-ADHoRe is basically the network as such

# under the transitivity assumption, we do not cluster the network, and we
# simply need to get a profile for each connected component; so we have to get
# the connected components from the anchorpoints file

import networkx as nx

def get_profile(anchorpoints_file, genes_data):
    # anchorpoints should be provided as file, no use to have a df in memory
    sn = get_network(anchorpoints_file)
    cc = nx.connected_components(sn)
    clusters = []
    sps = genes_data.sp.unique()
    for component in cc:
        d = {sp: 0 for sp in sps}
        for sp in genes_data.loc[component].sp:
            d[sp] += 1
        clusters.append(d)
    return pd.DataFrame.from_dict(clusters)  # save with index=None


def get_network(anchorpoints):
    G = nx.Graph()
    with open(anchorpoints, "r") as f:
        for line in f:
            l = line.strip().split()
            if l[0] == "id":
                continue
            gene1, gene2 = l[3], l[4]
            if gene1 not in G:
                G.add_node(gene1)
            if gene2 not in G:
                G.add_node(gene2)
            G.add_edge(gene1, gene2)
    return G
