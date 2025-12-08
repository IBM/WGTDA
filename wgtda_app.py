"""
Module Docstrings to be completed
"""
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from src.web_app.wgtda_stats import wgtda_network_stats_from_df
from src.web_app.wgtda_vis import plot_gene_network_interactive
from src.web_app.wgtda_scale_free import build_scale_free_figure


st.set_page_config(page_title="WGTDA Dashboard", layout="wide")
st.title("WGTDA Dashboard")
st.caption("Render networks, explore topological stats, and check scale-free topology.")

with st.sidebar:
    st.header("Load interactions CSV")
    up = st.file_uploader("Upload interactions CSV", type=["csv"])
    geneset_col = st.text_input("Geneset column name", value="geneset")
    title = st.text_input("Network title", value="Topological Co-occurrence Network")
    run_btn = st.button("Build dashboard")

if not run_btn:
    st.info("⬅️ Upload your interactions CSV and click **Build dashboard**.")
    st.stop()

if not up:
    st.error("Please upload an interactions CSV.")
    st.stop()

interactions = pd.read_csv(up)

# Compute stats
with st.spinner("Computing network stats..."):
    per_gene_df, summary = wgtda_network_stats_from_df(interactions, geneset_col)

# Tabs
tab_network, tab_stats, tab_scalefree = st.tabs(["Network", "Topological Stats", "Scale-Free Check"])

# --- NETWORK ---
with tab_network:
    colA, colB = st.columns([3, 2], gap="large")
    with colA:
        st.subheader("Interactive Network")
        with st.spinner("Rendering WGTDA network..."):
            tmp_html = tempfile.NamedTemporaryFile(delete=False, suffix=".html")
            html_path = plot_gene_network_interactive(
                interactions, geneset_col=geneset_col, per_gene_df=per_gene_df,
                output_file=tmp_html.name, title=title
            )
            html_str = Path(html_path).read_text(encoding="utf-8")
        components.html(html_str, height=800, scrolling=True)

    with colB:
    # --- Top hub genes ---
        st.subheader("Top hub genes (by degree)")
        top_n = st.slider("How many hubs to show?", 5, 50, 10, step=1, key="n_hubs")
        hubs_df = (per_gene_df[["gene", "degree", "betweenness_centrality", "avg_lifespan"]]
                .sort_values("degree", ascending=False)
                .head(top_n)
                .reset_index(drop=True))
        st.dataframe(hubs_df, width="stretch")

        st.download_button("Download hubs CSV",
                        hubs_df.to_csv(index=False).encode("utf-8"),
                        file_name="wgtda_top_hubs.csv", mime="text/csv")
        
        st.subheader("Network Summary")
        st.json(summary)

        st.download_button("Download per-gene stats CSV",
                        per_gene_df.to_csv(index=False).encode("utf-8"),
                        file_name="wgtda_network_stats.csv", mime="text/csv")
        st.download_button("Download network HTML",
                        data=html_str, file_name="wgtda_network.html", mime="text/html")
        
# --- STATS ---
with tab_stats:
    st.subheader("Per-gene Topological Statistics")
    st.dataframe(per_gene_df, width="stretch", height=600)


# --- SCALE-FREE ---
with tab_scalefree:
    st.subheader("Scale-Free Topology (Interactive)")
    fig, tbl_deg_by_gene, ddist, fit = build_scale_free_figure(per_gene_df)
    st.plotly_chart(fig, width="stretch")
    st.markdown(
        f"**Tail fit:** k_min = `{fit.get('k_min')}`, alpha ≈ `{fit.get('alpha')}`, "
        f"KS = `{fit.get('ks')}`, tail size = `{fit.get('n_tail')}`"
    )
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Degree-by-gene (ascending, normalized)**")
        st.dataframe(tbl_deg_by_gene, width="stretch", height=300)
        st.download_button("Download degree-by-gene CSV",
                           tbl_deg_by_gene.to_csv(index=False).encode("utf-8"),
                           file_name="wgtda_degree_by_gene.csv", mime="text/csv")
    with c2:
        st.markdown("**Degree distribution**")
        st.dataframe(ddist, width="stretch", height=300)
        st.download_button("Download degree distribution CSV",
                           ddist.to_csv(index=False).encode("utf-8"),
                           file_name="wgtda_degree_distribution.csv", mime="text/csv")






