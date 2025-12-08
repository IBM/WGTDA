# wgtda_scale_free.py
from collections import Counter
from pathlib import Path
from typing import Tuple, Dict, Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import powerlaw  # <-- new import


def degree_tables(per_gene_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return per-gene ascending table (with degree_normalized) and degree distribution (with ccdf)."""
    if "gene" not in per_gene_df.columns or "degree" not in per_gene_df.columns:
        raise ValueError("per_gene_df must include 'gene' and 'degree'.")

    tbl = per_gene_df[["gene", "degree"]].copy()
    kmax = max(1, int(tbl["degree"].max()))
    tbl["degree_normalized"] = tbl["degree"] / kmax
    tbl = tbl.sort_values("degree", ascending=True).reset_index(drop=True)

    degs = per_gene_df["degree"].astype(int).tolist()
    N = len(degs)
    cnt = Counter(degs)
    dd = pd.DataFrame([(k, cnt[k], cnt[k] / N) for k in sorted(cnt.keys())],
                      columns=["degree", "n_genes", "p_k"]).sort_values("degree")
    dd["ccdf"] = (dd["n_genes"][::-1].cumsum()[::-1]) / N
    return tbl, dd


def powerlaw_tail_fit(degrees) -> Dict[str, Any]:
    """Use powerlaw package to fit discrete power law and return summary stats."""
    if len(degrees) == 0:
        return {"k_min": None, "alpha": None, "ks": None, "n_tail": 0}

    degrees = np.array(degrees)
    fit = powerlaw.Fit(degrees, discrete=True, verbose=False)

    alpha = fit.power_law.alpha
    kmin = fit.power_law.xmin
    ks = fit.power_law.KS()
    n_tail = np.sum(degrees >= kmin)

    # Prepare values for plotting
    vals = np.unique(degrees)

    # ---- Empirical CCDF ----
    sorted_degs = np.sort(degrees)
    ccdf_emp = np.array([np.sum(sorted_degs >= v) / len(sorted_degs) for v in vals])

    # ---- Model CCDF ----
    ccdf_model = fit.power_law.ccdf(vals)

    return {
        "alpha": alpha,
        "k_min": kmin,
        "ks": ks,
        "n_tail": n_tail,
        "vals": vals,
        "ccdf_emp": ccdf_emp,
        "ccdf_model": ccdf_model,
        "fit_obj": fit
    }



def build_scale_free_figure(per_gene_df: pd.DataFrame):
    """Return (plotly Figure, per-gene ascending table, degree distribution, fit dict)."""
    tbl, dd = degree_tables(per_gene_df)
    fit = powerlaw_tail_fit(per_gene_df["degree"].astype(int).tolist())

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Normalized Avg #Connections vs #Genes",
                        f"CCDF (log–log). Tail size n={fit.get('n_tail',0)}, KS={fit.get('ks',float('nan')):.3f}"),
        horizontal_spacing=0.12,
    )

    # ---- Left plot: degree distribution ----
    fig.add_trace(go.Scatter(x=dd["degree"], y=dd["p_k"], mode="lines+markers",
                             name="p(k)", hovertemplate="k=%{x}<br>p(k)=%{y:.4f}<extra></extra>"),
                  row=1, col=1)
    fig.update_xaxes(title_text="Degree k (number of connections)", row=1, col=1)
    fig.update_yaxes(title_text="Normalized frequency p(k)", row=1, col=1)

    # ---- Right plot: empirical vs power-law CCDF ----
    fig.add_trace(go.Scatter(x=fit["vals"], y=fit["ccdf_emp"], mode="markers",
                             name="Empirical CCDF", hovertemplate="k=%{x}<br>P(K≥k)=%{y:.4f}<extra></extra>"),
                  row=1, col=2)

    if fit.get("k_min") is not None:
        fig.add_trace(go.Scatter(x=fit["vals"], y=fit["ccdf_model"], mode="lines",
                                 name=f"Power-law fit (k≥{fit['k_min']}, α≈{fit['alpha']:.2f})",
                                 hoverinfo="skip"),
                      row=1, col=2)

    fig.update_xaxes(title_text="Degree k", type="log", row=1, col=2)
    fig.update_yaxes(title_text="P(K ≥ k)", type="log", row=1, col=2)

    fig.update_layout(title_text="WGTDA Scale-Free Check (Interactive)",
                      showlegend=True, height=520, width=1000,
                      margin=dict(l=40, r=20, t=60, b=40))
    return fig, tbl, dd, fit


def save_plotly_html(fig, path: str) -> str:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(path, include_plotlyjs="cdn", full_html=True)
    return str(Path(path).resolve())
