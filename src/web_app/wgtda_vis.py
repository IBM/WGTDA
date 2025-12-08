# wgtda_vis.py
import html
import itertools
import json
from collections import defaultdict
from pathlib import Path
from typing import Optional

import networkx as nx
import pandas as pd
from pyvis.network import Network

from .wgtda_stats import _parse_list  # reuse the same parser


def plot_gene_network_interactive(
    df: pd.DataFrame,
    geneset_col: str = "geneset",
    per_gene_df: Optional[pd.DataFrame] = None,
    output_file: str = "wgtda_network.html",
    title: str = "WGTDA Network (with Topological Stats)",
    highlight_color: str = "#6EC1E4",
    scale_min: int = 10,
    scale_max: int = 30,
) -> str:
    """
    Build a PyVis network with NO hover tooltips. Clicking a node/edge updates a
    right-hand stats panel (implemented in injected JS). Writes and returns the HTML path.
    """
    if geneset_col not in df.columns:
        raise ValueError(f"Column '{geneset_col}' not found in df.")
    df = df.copy()
    df[geneset_col] = df[geneset_col].apply(_parse_list)

    # Edge weights + mean lifespan on edges
    edge_w = defaultdict(int)
    edge_ls = defaultdict(list)
    for _, row in df.iterrows():
        genes = list(dict.fromkeys([g for g in row[geneset_col] if isinstance(g, str) and g]))
        try:
            ls = float(row.get("lifespan", float("nan")))
        except Exception:
            ls = float("nan")
        for u, v in itertools.combinations(sorted(genes), 2):
            edge_w[(u, v)] += 1
            if ls == ls:
                edge_ls[(u, v)].append(ls)

    G = nx.Graph()
    for (u, v), w in edge_w.items():
        m = (sum(edge_ls[(u, v)]) / len(edge_ls[(u, v)])) if edge_ls[(u, v)] else None
        G.add_edge(u, v, weight=w, mean_lifespan=(None if m is None else float(m)))

    node_deg = dict(G.degree())
    max_deg = max(node_deg.values()) if node_deg else 1
    for n in G.nodes():
        G.nodes[n]["size"] = scale_min + (node_deg.get(n, 0) / max_deg) * (scale_max - scale_min)
        G.nodes[n]["color"] = highlight_color  # no 'title' => no hover popups

    # Per-gene stats lookup for JS
    per_gene_stats = {}
    if per_gene_df is not None and not per_gene_df.empty:
        if "gene" not in per_gene_df.columns:
            raise ValueError("per_gene_df must include a 'gene' column.")
        per_gene_stats = (
            per_gene_df.set_index("gene")
            .drop(columns=[c for c in ["gene"] if c in per_gene_df.columns], errors="ignore")
            .to_dict(orient="index")
        )

    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black",
                  notebook=False, cdn_resources="in_line")
    net.from_nx(G)
    net.toggle_physics(True)

    tmp = "_tmp_network.html"
    net.write_html(tmp)

    NODE_STATS_JSON = json.dumps(
        {k: {kk: vv for kk, vv in v.items()
             if (isinstance(vv, (str, int)) or (isinstance(vv, float) and vv == vv))}
         for k, v in per_gene_stats.items()}
    )
    EDGE_STATS = {}
    for (u, v), w in edge_w.items():
        key = f"{u}|{v}"
        m = (sum(edge_ls[(u, v)]) / len(edge_ls[(u, v)])) if edge_ls[(u, v)] else None
        EDGE_STATS[key] = {"from": u, "to": v, "cooccurrence_weight": w,
                           "mean_lifespan": (None if m is None else float(m))}
    EDGE_STATS_JSON = json.dumps(EDGE_STATS)

    html_content = Path(tmp).read_text(encoding="utf-8")

    # Title
    html_content = html_content.replace(
        "<body>",
        "<body>\n"
        "<h1 style='text-align:center;font-family:Arial,sans-serif;margin:16px 0 8px 0;'>"
        + html.escape(title) + "</h1>\n"
    )

    # Side-by-side layout (graph | panel)
    html_content = html_content.replace(
        '<div id="mynetwork"',
        "<div style='display:flex;flex-direction:row;align-items:flex-start;gap:16px;'>\n"
        "  <div id='mynetwork' style='flex:1 1 auto; min-width:0; height:750px;'"
    ).replace(
        "</div>\n</body>",
        "  </div>\n"
        "  <aside id='statsPanel' style='flex:0 0 320px;max-width:320px;border-left:1px solid #eee;"
        "padding:12px;font-family:Arial,sans-serif;overflow-y:auto;'>\n"
        "    <div style='font-weight:bold;margin-bottom:8px;'>Details</div>\n"
        "    <div id='statsBody' style='font-size:14px;color:#333;'>\n"
        "      <em>Click a node (gene) or an edge to see topological statistics.</em>\n"
        "    </div>\n"
        "  </aside>\n"
        "</div>\n"
        "</body>"
    )

    # Minimal JS (concat to avoid f-string brace issues)
    js = (
        "<script>\n"
        "const NODE_STATS = " + NODE_STATS_JSON + ";\n"
        "const EDGE_STATS = " + EDGE_STATS_JSON + ";\n"
        "const fmt = v => (v===null || Number.isNaN(v)) ? '-' : "
        "(typeof v==='number' ? (''+v.toFixed(6)).replace(/0+$/,'').replace(/\\.$/,'') : v);\n"
        "const table = d => { const es = Object.entries(d); if (!es.length) return '<em>No additional stats.</em>';\n"
        "  return '<table style=\"border-collapse:collapse;\">' + es.map(([k,v]) => "
        "'<tr><td style=\"padding:4px 8px 4px 0;color:#666;\">' + k.replace(/_/g,' ').replace(/\\b\\w/g, c=>c.toUpperCase()) + "
        "'</td><td style=\"padding:4px 0;font-weight:600;\">' + fmt(v) + '</td></tr>').join('') + '</table>'; };\n"
        "(function ready(){ if(!window.network||!network.body) return setTimeout(ready,80);\n"
        "  const edges=network.body.data.edges; const panel=document.getElementById('statsBody');\n"
        "  network.on('click',(p)=>{ if(p.nodes&&p.nodes.length){ const id=p.nodes[0];\n"
        "    panel.innerHTML='<div style=\"font-size:16px;font-weight:700;margin-bottom:6px;\">Gene: '+String(id).replace(/</g,'&lt;').replace(/>/g,'&gt;')+'</div>'+table(NODE_STATS[id]||{});\n"
        "  }else if(p.edges&&p.edges.length){ const ed=edges.get(p.edges[0]); if(!ed) return; const key=[ed.from,ed.to].sort().join('|');\n"
        "    panel.innerHTML='<div style=\"font-size:16px;font-weight:700;margin-bottom:6px;\">Edge: '+ed.from+' — '+ed.to+'</div>'+table(EDGE_STATS[key]||{}); }});\n"
        "})();\n"
        "</script>\n"
    )
    html_content = html_content.replace("</body>", js + "</body>")

    Path(output_file).write_text(html_content, encoding="utf-8")
    return str(Path(output_file).resolve())
