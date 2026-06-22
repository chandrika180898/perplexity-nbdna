"""Perplexity Profile page – interactive chart and region table."""

import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Perplexity Profile | Non-B DNA",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("📊 Perplexity Profile")

# ── Guard ─────────────────────────────────────────────────────────────────────
if "analysis_done" not in st.session_state:
    st.warning("No analysis results yet. Please run the **Analysis** step first.", icon="⚠️")
    st.page_link("pages/1_🧬_Analysis.py", label="Go to Analysis →", icon="🧬")
    st.stop()

perps   = st.session_state["perps"]       # np.ndarray
regions = st.session_state["regions"]     # list[(start, end)]
params  = st.session_state["params"]
seq_len = len(st.session_state["seq"])

# ── Display controls ──────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔎 View Controls")
    subsample = st.slider(
        "Subsample every N windows (for speed)",
        min_value=1,
        max_value=max(1, len(perps) // 2000),
        value=max(1, len(perps) // 5000),
        help="Reduce the number of points plotted without affecting the analysis.",
    )
    show_threshold = st.checkbox("Show threshold line", value=True)
    show_regions = st.checkbox("Highlight low-perplexity regions", value=True)

# ── Build chart data ──────────────────────────────────────────────────────────
positions = np.arange(len(perps), dtype=np.int32)

# Subsample for rendering performance
idx = np.arange(0, len(perps), max(1, subsample))
plot_df = pd.DataFrame(
    {
        "Position (nt)": positions[idx],
        "Perplexity": perps[idx].astype(float),
    }
).set_index("Position (nt)")

threshold_val = float(np.percentile(perps, params["threshold_pct"]))

# ── Main perplexity line chart ────────────────────────────────────────────────
st.subheader("Nucleotide Perplexity Landscape")
st.caption(
    f"Window = {params['window']} nt · "
    f"Threshold = {threshold_val:.3f} "
    f"(p{params['threshold_pct']}) · "
    f"Sequence = {seq_len:,} bp"
)

st.line_chart(plot_df, y="Perplexity", height=350, use_container_width=True)

# ── Threshold & region annotations (text fallback) ───────────────────────────
if show_threshold:
    st.caption(
        f"🔵 Perplexity threshold line: **{threshold_val:.4f}** "
        f"(p{params['threshold_pct']} of all windows)"
    )

if show_regions and regions:
    # Highlight bars using st.dataframe with color formatting
    region_df = pd.DataFrame(regions, columns=["Start", "End"])
    region_df["Length (nt)"] = region_df["End"] - region_df["Start"]
    region_df["Mean Perplexity"] = [
        # Map genomic region [s, e) to perplexity-array indices.
        # Window i covers positions [i, i+window), so the last window whose
        # start falls inside [s, e) has index (e - window).  We clamp the
        # slice to at least one element (s+1) to avoid empty slices.
        float(perps[s : max(s + 1, e - params["window"] + 1)].mean())
        if s < len(perps)
        else float("nan")
        for s, e in regions
    ]

    st.divider()
    st.subheader(f"🔴 Low-Perplexity Regions  ({len(regions)} total)")
    st.dataframe(
        region_df.style.format(
            {
                "Mean Perplexity": "{:.4f}",
                "Start": "{:,}",
                "End": "{:,}",
                "Length (nt)": "{:,}",
            }
        ).background_gradient(subset=["Mean Perplexity"], cmap="YlOrRd_r"),
        use_container_width=True,
        hide_index=True,
    )

    # Aggregate stats
    col1, col2, col3 = st.columns(3)
    lengths = region_df["Length (nt)"]
    col1.metric("Total bp in low-perplexity regions", f"{int(lengths.sum()):,}")
    col2.metric("Median region length", f"{int(lengths.median()):,} nt")
    col3.metric("Longest region", f"{int(lengths.max()):,} nt")

elif not regions:
    st.info("No low-perplexity regions detected with the current parameters.")

# ── Perplexity distribution histogram ─────────────────────────────────────────
st.divider()
st.subheader("📈 Perplexity Distribution")

hist_vals, bin_edges = np.histogram(perps, bins=50)
hist_df = pd.DataFrame(
    {
        "Perplexity bin": (bin_edges[:-1] + bin_edges[1:]) / 2,
        "Window count": hist_vals,
    }
).set_index("Perplexity bin")

st.bar_chart(hist_df, y="Window count", height=280, use_container_width=True)

col_a, col_b, col_c, col_d = st.columns(4)
col_a.metric("Min perplexity", f"{float(perps.min()):.4f}")
col_b.metric("Max perplexity", f"{float(perps.max()):.4f}")
col_c.metric("Mean perplexity", f"{float(perps.mean()):.4f}")
col_d.metric("Std dev", f"{float(perps.std()):.4f}")
