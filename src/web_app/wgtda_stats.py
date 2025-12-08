# wgtda_stats.py
import ast
import itertools
from typing import Tuple, Dict, Any

import networkx as nx
import pandas as pd


def _parse_list(v):
    if isinstance(v, list):
        return v
    if isinstance(v, str):
        try:
            return ast.literal_eval(v)
        except Exception:
            return [x.strip() for x in v.strip("[]").split(",") if x.strip()]
    return []


def wgtda_network_stats_from_df(
    interactions_df: pd.DataFrame, gene_list_col: str = "geneset"
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Build a co-occurrence graph from WGTDA interactions and compute per-gene + network stats.
    interactions_df must contain a column of gene lists (e.g. 'geneset'); optional 'lifespan'.
    Returns (per_gene_df, summary_dict).
    """
    if gene_list_col not in interactions_df.columns:
        raise ValueError(f"Column '{gene_list_col}' not found in dataframe.")

    df = interactions_df.copy()
    df[gene_list_col] = df[gene_list_col].apply(_parse_list)

    G = nx.Graph()
    lifespan_by_gene = {}

    for _, row in df.iterrows():
        genes = [g for g in row[gene_list_col] if isinstance(g, str) and g]
        genes = list(dict.fromkeys(genes))
        for g in genes:
            G.add_node(g)
        for u, v in itertools.combinations(sorted(genes), 2):
            if G.has_edge(u, v):
                G[u][v]["weight"] = G[u][v].get("weight", 1) + 1
            else:
                G.add_edge(u, v, weight=1)

        ls = row.get("lifespan", None)
        try:
            ls = float(ls)
        except Exception:
            ls = None
        if ls is not None:
            for g in genes:
                lifespan_by_gene.setdefault(g, []).append(ls)

    # Per-gene stats
    degree = pd.Series(dict(G.degree()), name="degree").astype(int)
    num_connections = degree.rename("num_connections")
    deg_cent = pd.Series(nx.degree_centrality(G), name="degree_centrality")
    bet_cent = pd.Series(nx.betweenness_centrality(G, normalized=True), name="betweenness_centrality")
    clo_cent = pd.Series(nx.closeness_centrality(G), name="closeness_centrality")

    avg_sp = {}
    for comp in nx.connected_components(G):
        sub = G.subgraph(comp)
        lengths = dict(nx.all_pairs_shortest_path_length(sub))
        for n in sub.nodes():
            dists = [d for tgt, d in lengths[n].items() if tgt != n]
            avg_sp[n] = float(sum(dists) / len(dists)) if dists else float("nan")
    avg_sp = pd.Series(avg_sp, name="avg_shortest_path_length")

    avg_lifespan = pd.Series(
        {g: (sum(v) / len(v)) for g, v in lifespan_by_gene.items()},
        name="avg_lifespan",
    )

    comp_id = {}
    comp_size = {}
    for i, component in enumerate(nx.connected_components(G), start=1):
        for n in component:
            comp_id[n] = i
            comp_size[n] = len(component)

    comp_id = pd.Series(comp_id, name="component_id")
    comp_size = pd.Series(comp_size, name="component_size")

    per_gene = (
        pd.concat(
            [
                degree,
                num_connections,
                deg_cent,
                bet_cent,
                clo_cent,
                avg_sp,
                avg_lifespan.reindex(G.nodes(), fill_value=float("nan")),
                comp_id,
                comp_size,
            ],
            axis=1,
        )
        .reset_index()
        .rename(columns={"index": "gene"})
        .sort_values(["component_size", "degree"], ascending=[False, False])
        .reset_index(drop=True)
    )

    # Network summary (on largest connected component for APL)
    if G.number_of_nodes() > 0:
        lcc_nodes = max(nx.connected_components(G), key=len)
        L = G.subgraph(lcc_nodes)
        try:
            global_apl = nx.average_shortest_path_length(L)
        except Exception:
            global_apl = float("nan")
        global_clust = nx.average_clustering(G)
    else:
        global_apl, global_clust, lcc_nodes = float("nan"), float("nan"), []

    summary = dict(
        num_nodes=G.number_of_nodes(),
        num_edges=G.number_of_edges(),
        global_avg_shortest_path_length=global_apl,
        global_avg_clustering=global_clust,
        degree_distribution=per_gene["degree"].value_counts().sort_index().to_dict(),
        num_components=nx.number_connected_components(G),
        largest_component_size=len(lcc_nodes) if G.number_of_nodes() else 0,
    )

    return per_gene, summary
