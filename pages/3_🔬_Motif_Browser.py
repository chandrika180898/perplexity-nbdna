"""Motif Browser page – filterable table of detected Non-B motifs + CSV download."""

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Motif Browser | Non-B DNA",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🔬 Motif Browser")
st.markdown(
    "Explore Non-B DNA motifs detected inside low-perplexity regions.  "
    "Use the filters below to narrow down results and download a CSV report."
)

# ── Guard ─────────────────────────────────────────────────────────────────────
if "analysis_done" not in st.session_state:
    st.warning("No analysis results yet. Please run the **Analysis** step first.", icon="⚠️")
    st.page_link("pages/1_🧬_Analysis.py", label="Go to Analysis →", icon="🧬")
    st.stop()

motif_results: list[dict] = st.session_state["motif_results"]

if not motif_results:
    st.warning(
        "No Non-B motifs were detected in low-perplexity regions.  "
        "Try reducing the **minimum region length** or adjusting the **threshold percentile** "
        "in the Analysis page.",
        icon="🔍",
    )
    st.stop()

df = pd.DataFrame(motif_results)

# ── Sidebar filters ───────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🗂️ Filters")

    motif_types = sorted(df["Motif"].unique().tolist())
    selected_motifs = st.multiselect(
        "Motif type",
        options=motif_types,
        default=motif_types,
        help="Show only selected motif categories.",
    )

    min_region_len = int(df["Region_Length_nt"].min())
    max_region_len = int(df["Region_Length_nt"].max())
    region_len_range = st.slider(
        "Region length (nt)",
        min_value=min_region_len,
        max_value=max_region_len,
        value=(min_region_len, max_region_len),
        help="Filter by the length of the enclosing low-perplexity region.",
    )

    min_motif_len = int((df["Motif_End"] - df["Motif_Start"]).min())
    max_motif_len = int((df["Motif_End"] - df["Motif_Start"]).max())
    motif_len_range = st.slider(
        "Motif length (nt)",
        min_value=min_motif_len,
        max_value=max_motif_len,
        value=(min_motif_len, max_motif_len),
        help="Filter by the length of the individual motif hit.",
    )

# ── Apply filters ─────────────────────────────────────────────────────────────
df["Motif_Length_nt"] = df["Motif_End"] - df["Motif_Start"]

mask = (
    df["Motif"].isin(selected_motifs)
    & df["Region_Length_nt"].between(*region_len_range)
    & df["Motif_Length_nt"].between(*motif_len_range)
)
filtered = df[mask].reset_index(drop=True)

# ── Summary bar ───────────────────────────────────────────────────────────────
c1, c2, c3 = st.columns(3)
c1.metric("Total hits (all filters)", len(df))
c2.metric("Filtered hits", len(filtered))
c3.metric("Motif types shown", len(selected_motifs))

st.divider()

# ── Motif-type distribution ───────────────────────────────────────────────────
st.subheader("📊 Motif Type Distribution")
counts = (
    filtered.groupby("Motif", as_index=False)
    .size()
    .rename(columns={"size": "Count"})
    .sort_values("Count", ascending=False)
    .set_index("Motif")
)
st.bar_chart(counts, y="Count", height=260, use_container_width=True)

st.divider()

# ── Results table ─────────────────────────────────────────────────────────────
st.subheader(f"🗒️ Motif Hits  ({len(filtered):,} rows)")

display_cols = [
    "Motif",
    "Motif_Start",
    "Motif_End",
    "Motif_Length_nt",
    "Motif_Sequence",
    "LowP_Region_Start",
    "LowP_Region_End",
    "Region_Length_nt",
]

st.dataframe(
    filtered[display_cols].style.format(
        {
            "Motif_Start": "{:,}",
            "Motif_End": "{:,}",
            "LowP_Region_Start": "{:,}",
            "LowP_Region_End": "{:,}",
            "Region_Length_nt": "{:,}",
        }
    ),
    use_container_width=True,
    hide_index=True,
)

# ── Download ──────────────────────────────────────────────────────────────────
st.download_button(
    label="⬇️  Download filtered results (CSV)",
    data=filtered[display_cols].to_csv(index=False),
    file_name="nonb_motifs_filtered.csv",
    mime="text/csv",
    use_container_width=True,
)

st.download_button(
    label="⬇️  Download all results (CSV)",
    data=df[display_cols].to_csv(index=False),
    file_name="nonb_motifs_all.csv",
    mime="text/csv",
    use_container_width=True,
)

# ── Detail view ───────────────────────────────────────────────────────────────
st.divider()
st.subheader("🔍 Inspect a Hit")
if len(filtered) > 0:
    hit_idx = st.number_input(
        "Row index (0-based)",
        min_value=0,
        max_value=max(0, len(filtered) - 1),
        value=0,
        step=1,
    )
    row = filtered.iloc[hit_idx]
    seq_snippet = st.session_state["seq"][row["LowP_Region_Start"]: row["LowP_Region_End"]]
    motif_rel_start = row["Motif_Start"] - row["LowP_Region_Start"]
    motif_rel_end = row["Motif_End"] - row["LowP_Region_Start"]

    st.markdown(f"**Motif type:** `{row['Motif']}`")
    st.markdown(
        f"**Genomic position:** {int(row['Motif_Start']):,} – {int(row['Motif_End']):,} nt"
    )
    st.markdown(
        f"**Enclosing region:** {int(row['LowP_Region_Start']):,} – "
        f"{int(row['LowP_Region_End']):,} nt  ({int(row['Region_Length_nt']):,} nt)"
    )

    # Highlight motif within region sequence
    highlighted = (
        seq_snippet[:motif_rel_start]
        + f"**[{seq_snippet[motif_rel_start:motif_rel_end]}]**"
        + seq_snippet[motif_rel_end:]
    )
    with st.expander("Region sequence (motif in **[bold brackets]**)", expanded=True):
        st.code(highlighted, language=None)
