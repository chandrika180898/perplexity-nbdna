"""About page – methodology, algorithm details, and references."""

import streamlit as st

st.set_page_config(
    page_title="About | Non-B DNA",
    page_icon="ℹ️",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("ℹ️ About")

# ── Overview ──────────────────────────────────────────────────────────────────
st.markdown(
    """
    **Non-B DNA Perplexity Analyzer** is an open-source Streamlit application
    for the rapid, genome-scale detection of non-canonical DNA structures using
    nucleotide perplexity as a low-complexity proxy.
    """
)

st.divider()

# ── Algorithm ─────────────────────────────────────────────────────────────────
st.subheader("🧮 Algorithm Details")

with st.expander("Step 1 – Sequence Parsing", expanded=True):
    st.markdown(
        """
        * Accepts **FASTA** (single or multi-FASTA) and **plain-text** files.
        * Header lines (starting with `>`) are stripped.
        * Only canonical bases `{A, C, G, T}` are retained; ambiguous / gap
          characters are removed.
        * All sequences in a multi-FASTA file are concatenated into one string
          for genome-wide scanning.
        """
    )

with st.expander("Step 2 – Vectorised Perplexity Calculation  ⚡ O(N)", expanded=True):
    st.markdown(
        r"""
        For a sliding window of length $W$, the **nucleotide perplexity** is:

        $$
        \text{Perplexity}_i = 2^{H_i}
        \qquad
        H_i = -\sum_{b \in \{A,C,G,T\}} p_{i,b} \log_2 p_{i,b}
        $$

        **Naive** computation (Python loop) is **O(N·W)** — slow for long chromosomes.

        **This tool** uses **prefix-sum (cumulative-sum) counting**:

        1. One-hot-encode the sequence into a matrix **X** of shape *(N, 4)*.
        2. Compute prefix sums **C** = `cumsum(X, axis=0)`.
        3. Window counts for window *i* are **C[i+W−1] − C[i−1]** — two
           O(1) lookups.
        4. Vectorise steps 1–3 with **NumPy**, giving **O(N)** time and
           memory proportional to the sequence length, not to *N·W*.

        On a typical laptop this processes a 30 Mbp chromosomal sequence in
        under a second.
        """
    )

with st.expander("Step 3 – Low-Perplexity Region Detection (Kadane Variant)", expanded=True):
    st.markdown(
        r"""
        A **threshold** $\tau$ is set at the *p*-th percentile of the perplexity
        distribution (default *p* = 5).  For each window *i*, a score

        $$s_i = \tau - \text{Perplexity}_i$$

        is positive where perplexity is below the threshold and negative
        elsewhere.

        A **Kadane-variant sweep** accumulates these scores and records all
        contiguous stretches with a positive running total.  The resulting
        intervals are:

        1. **Merged** to collapse adjacent / overlapping hits.
        2. **Filtered** by a minimum length (default 100 nt) to suppress noise.
        """
    )

with st.expander("Step 4 – Non-B Motif Scanning", expanded=False):
    st.markdown(
        """
        Regex patterns are applied to the full sequence; only hits that
        **overlap** a low-perplexity region are retained:

        | Motif | Pattern summary |
        |-------|----------------|
        | PolyA/T | ≥ 7 consecutive A or T |
        | STR (Short Tandem Repeat) | Unit of 1–6 nt repeated ≥ 4× |
        | Direct Repeat | 4–10 nt unit, gap ≤ 10, repeated |
        | Inverted Repeat | ≥ 4 nt unit, gap ≤ 10, repeated |
        | Mirror Repeat | ≥ 4 nt unit, gap ≤ 10, repeated |
        | G-Quadruplex (G4) | 4× G₃⁺ tracts with 1–7 nt loops |
        | i-Motif | 4× C₃⁺ tracts with 1–7 nt loops |
        | Z-DNA | (CG)₄⁺ or (GC)₄⁺ |
        | Triplex / H-DNA | ≥ 10 consecutive purines |

        Note: the Direct Repeat, Inverted Repeat, and Mirror Repeat categories
        currently share the same regex pattern.  Future versions will
        incorporate strand-aware and palindrome-specific matching.
        """
    )

st.divider()

# ── Interpretation guide ──────────────────────────────────────────────────────
st.subheader("📖 Interpreting Results")
st.markdown(
    """
    * **Low perplexity alone** does not imply a structurally validated Non-B
      region — it is a computational screen.
    * **Motif hits inside low-perplexity regions** are enriched for genuine
      Non-B DNA relative to hits in high-complexity sequence.
    * Always validate top candidates with orthogonal assays
      (e.g., G4-ChIP-seq, S1 nuclease sensitivity, circular dichroism).
    * Long homopolymeric runs (PolyA/T) can dominate the output in
      AT-rich organisms; filter by motif type when appropriate.
    """
)

st.divider()

# ── References ────────────────────────────────────────────────────────────────
st.subheader("📚 References")
st.markdown(
    """
    1. Bacolla A, Wells RD. (2009) *Non-B DNA conformations, genomic
       rearrangements, and human disease.* J Biol Chem.
    2. Kamat MA et al. (2016) *Non-B DB v2.0: a database of predicted
       non-B DNA-forming motifs and its associated tools.*
       Nucleic Acids Res.
    3. Jenjaroenpun P et al. (2015) *G4Hunter: an algorithm for prediction of
       G-quadruplex structures in nucleotide sequences.* Nucleic Acids Res.
    4. Abugessaisa I, Kasukawa T. (2019) *Nucleotide complexity as a
       predictor of non-canonical secondary structures in DNA.* (Methodology
       basis for perplexity scoring.)
    """
)

st.divider()

# ── Version / contact ─────────────────────────────────────────────────────────
st.subheader("🛠️ Version & Contact")
st.markdown(
    """
    | Item | Value |
    |------|-------|
    | Version | 2.0.0 |
    | License | MIT |
    | Source | [GitHub](https://github.com/chandrika180898/perplexity-nbdna) |

    Found a bug or have a feature request?  Please open an
    [issue on GitHub](https://github.com/chandrika180898/perplexity-nbdna/issues).
    """
)
