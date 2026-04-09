import streamlit as st

st.set_page_config(
    page_title="Non-B DNA Perplexity Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Hero ──────────────────────────────────────────────────────────────────────
st.title("🧬 Non-B DNA Perplexity Analyzer")
st.markdown(
    """
    <p style='font-size:1.15rem; color:#555;'>
    A high-performance web tool for detecting <b>low-complexity / Non-B DNA</b>
    structures using nucleotide perplexity and motif analysis.
    </p>
    """,
    unsafe_allow_html=True,
)

st.divider()

# ── Overview cards ────────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.info("### 🧪 Upload\nFASTA / TXT sequences of any length")
with col2:
    st.success("### ⚡ Analyse\nVectorised NumPy perplexity — runs in O(N)")
with col3:
    st.warning("### 📊 Visualise\nInteractive perplexity profile with region highlights")
with col4:
    st.error("### 🔬 Browse\nFilter & download detected Non-B motifs")

st.divider()

# ── Quick-start guide ─────────────────────────────────────────────────────────
st.subheader("🚀 Quick Start")
st.markdown(
    """
    1. **Analysis** – upload your sequence file (FASTA / `.fa` / `.txt`) and
       tune the sliding-window parameters.
    2. **Perplexity Profile** – explore the nucleotide-perplexity landscape and
       all identified low-complexity regions.
    3. **Motif Browser** – filter the detected Non-B motifs, inspect individual
       hits, and download a CSV report.
    4. **About** – read about the algorithm, Non-B DNA biology, and references.

    > Use the **sidebar** (left) to navigate between pages.
    """
)

st.divider()

# ── What is Non-B DNA? ────────────────────────────────────────────────────────
with st.expander("📖 What is Non-B DNA?", expanded=False):
    st.markdown(
        """
        **Non-B DNA** refers to any DNA secondary structure that deviates from
        the canonical right-handed B-form double helix.  These structures arise
        at sequence motifs with special repeat or base-composition properties and
        include:

        | Structure | Typical motif |
        |-----------|--------------|
        | G-Quadruplex (G4) | `G₃⁺NₓG₃⁺NₓG₃⁺NₓG₃⁺` |
        | i-Motif | C-rich counterpart of G4 |
        | Z-DNA | `(CG)ₙ` alternating purine–pyrimidine |
        | Triplex / H-DNA | Homopurine–homopyrimidine tracts |
        | Inverted / Mirror / Direct Repeats | Various repeat units |
        | Short Tandem Repeats (STR) | `(unit)₄⁺` |
        | PolyA / PolyT | `A₇⁺` or `T₇⁺` |

        Non-B DNA regions are enriched at:
        * Replication origins and fragile sites
        * Promoters and regulatory elements
        * Cancer-associated genomic rearrangement break-points

        **Low nucleotide perplexity** (i.e., low sequence complexity) is a
        fast proxy for locating these regions before applying motif-level regex
        detection.
        """
    )

# ── How perplexity is calculated ──────────────────────────────────────────────
with st.expander("🧮 How is perplexity calculated?", expanded=False):
    st.markdown(
        r"""
        For a DNA window of length $W$, the **nucleotide perplexity** is
        defined as the exponentiated Shannon entropy over the four-nucleotide
        alphabet $\{A, C, G, T\}$:

        $$
        H = -\sum_{b \in \{A,C,G,T\}} p_b \log_2 p_b
        \qquad
        \text{Perplexity} = 2^{H}
        $$

        * **Perplexity = 4** → perfectly uniform composition (maximum complexity)
        * **Perplexity → 1** → single nucleotide dominates (minimum complexity)

        The tool computes perplexity for every overlapping window using
        **vectorised cumulative-sum counting** (NumPy), reducing the
        complexity from O(N·W) to **O(N)**.  Low-perplexity regions are then
        identified with a **Kadane-variant algorithm** and filtered by minimum
        length before motif scanning.
        """
    )

st.info(
    "👈  Navigate to **Analysis** in the sidebar to upload your sequence and start.",
    icon="▶️",
)
