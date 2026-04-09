"""Analysis page – sequence upload, parameter tuning, and vectorised analysis."""

import re
import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Analysis | Non-B DNA",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🧬 Sequence Analysis")
st.markdown(
    "Upload a FASTA / TXT file, adjust the parameters, and run the "
    "perplexity-based Non-B DNA detection pipeline."
)

# ── Helpers ───────────────────────────────────────────────────────────────────

def read_sequence(uploaded_file) -> str:
    """Parse FASTA / plain-text file, return clean uppercase ACGT string."""
    text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    seq_lines = [
        line.strip().upper()
        for line in text.splitlines()
        if not line.startswith(">")
    ]
    return re.sub(r"[^ACGT]", "", "".join(seq_lines))


@st.cache_data(show_spinner=False)
def calc_perplexity_vectorized(seq: str, window: int) -> np.ndarray:
    """
    Compute per-window perplexity in O(N) using cumulative-sum counting.

    For each overlapping window of length `window`, the nucleotide
    frequencies are derived from prefix sums so the inner loop is avoided.
    Returns an ndarray of shape (N - window + 1,).
    """
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    seq_arr = np.frombuffer(
        seq.encode("ascii"), dtype=np.uint8
    )  # raw ASCII bytes

    n = len(seq_arr)
    if n < window:
        return np.array([], dtype=np.float32)

    # One-hot encode: shape (n, 4)
    one_hot = np.zeros((n, 4), dtype=np.float32)
    for idx, char in enumerate("ACGT"):
        one_hot[:, idx] = seq_arr == ord(char)

    # Prefix sums for O(1) range queries
    cum = np.cumsum(one_hot, axis=0)  # shape (n, 4)
    # Window counts: cum[i+window-1] - cum[i-1]
    counts_end = cum[window - 1:]                          # shape (N-W+1, 4)
    counts_start = np.vstack([np.zeros((1, 4), dtype=np.float32), cum[: n - window]])
    counts = counts_end - counts_start                     # shape (N-W+1, 4)

    # Frequencies
    totals = counts.sum(axis=1, keepdims=True)             # always == window
    probs = np.where(counts > 0, counts / totals, 0.0)

    # Shannon entropy → perplexity
    log_p = np.where(probs > 0, np.log2(probs), 0.0)
    entropy = -(probs * log_p).sum(axis=1)
    return (2.0 ** entropy).astype(np.float32)


@st.cache_data(show_spinner=False)
def find_low_perplexity_regions(
    perps: np.ndarray, window: int, threshold_pct: int, min_len: int
) -> list[tuple[int, int]]:
    """
    Locate merged low-perplexity regions using a Kadane-variant sweep.

    Parameters
    ----------
    perps         : per-window perplexity array
    window        : window size (nt)
    threshold_pct : percentile of perplexity to use as threshold
    min_len       : minimum region length to retain (nt)
    """
    if len(perps) == 0:
        return []

    threshold = float(np.percentile(perps, threshold_pct))
    scores = threshold - perps  # positive where perplexity is low

    # Kadane sweep to find contiguous low-perplexity stretches
    regions: list[tuple[int, int]] = []
    current_sum = 0.0
    start = 0

    for i, score in enumerate(scores):
        if current_sum <= 0:
            start = i
            current_sum = float(score)
        else:
            current_sum += float(score)

        if current_sum > 0:
            regions.append((start, i + window))   # genomic coords

    if not regions:
        return []

    # Merge overlapping / adjacent regions
    regions.sort()
    merged = [list(regions[0])]
    for s, e in regions[1:]:
        last = merged[-1]
        if s <= last[1]:
            last[1] = max(last[1], e)
        else:
            merged.append([s, e])

    # Filter by minimum length
    return [(s, e) for s, e in merged if (e - s) >= min_len]


@st.cache_data(show_spinner=False)
def detect_motifs(
    seq: str, regions: list[tuple[int, int]]
) -> list[dict]:
    """Scan Non-B motif regex patterns against low-perplexity regions."""
    motif_patterns = {
        "PolyA/T":          r"A{7,}|T{7,}",
        "STR":               r"([ACGT]{1,6})\1{3,}",
        "DirectRepeat":      r"([ACGT]{4,10})[ACGT]{0,10}\1",
        "InvertedRepeat":    r"([ACGT]{4,})[ACGT]{0,10}\1",
        "MirrorRepeat":      r"([ACGT]{4,})[ACGT]{0,10}\1",
        "G4":                r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",
        "iMotif":            r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",
        "Z-DNA":             r"(CG){4,}|(GC){4,}",
        "Triplex_HDNA":      r"[AG]{10,}",
    }
    compiled = {k: re.compile(v) for k, v in motif_patterns.items()}

    results = []
    for name, rgx in compiled.items():
        for m in rgx.finditer(seq):
            ms, me = m.start(), m.end()
            for rs, re_ in regions:
                if max(ms, rs) < min(me, re_):
                    results.append(
                        {
                            "Motif": name,
                            "Motif_Start": ms,
                            "Motif_End": me,
                            "Motif_Sequence": m.group(),
                            "LowP_Region_Start": rs,
                            "LowP_Region_End": re_,
                            "Region_Length_nt": re_ - rs,
                        }
                    )
    return results


# ── UI ────────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.header("⚙️ Parameters")
    window_size = st.slider(
        "Sliding window (nt)", min_value=5, max_value=50, value=10, step=1,
        help="Length of the overlapping window used to compute per-position perplexity.",
    )
    threshold_pct = st.slider(
        "Low-perplexity threshold (percentile)", min_value=1, max_value=25, value=5,
        help="Windows below this percentile of the perplexity distribution are considered low-complexity.",
    )
    min_region_len = st.number_input(
        "Minimum region length (nt)", min_value=10, max_value=10000, value=100, step=10,
        help="Merged regions shorter than this value are discarded.",
    )

uploaded_file = st.file_uploader(
    "📂 Upload a FASTA or plain-text sequence file",
    type=["fa", "fasta", "txt"],
    help="Multi-FASTA files are supported; all sequences are concatenated.",
)

if uploaded_file is None:
    st.info("👆 Please upload a sequence file to begin.", icon="📂")
    st.stop()

# ── Parse ─────────────────────────────────────────────────────────────────────
with st.spinner("Reading sequence …"):
    seq = read_sequence(uploaded_file)

if len(seq) == 0:
    st.error("No valid ACGT characters found in the uploaded file.")
    st.stop()
if len(seq) < window_size:
    st.error(f"Sequence is too short ({len(seq)} bp). Must be ≥ {window_size} bp.")
    st.stop()

st.success(f"Sequence loaded: **{len(seq):,} bp**")

# ── Analysis ──────────────────────────────────────────────────────────────────
run = st.button("▶️ Run Analysis", type="primary", use_container_width=True)

if run or "analysis_done" in st.session_state:

    if run:
        # Clear cached results when parameters change
        st.session_state.pop("analysis_done", None)

        with st.spinner("Computing perplexity (vectorised) …"):
            perps = calc_perplexity_vectorized(seq, window_size)

        with st.spinner("Detecting low-perplexity regions …"):
            regions = find_low_perplexity_regions(
                perps, window_size, threshold_pct, min_region_len
            )

        with st.spinner("Scanning Non-B motifs …"):
            motif_results = detect_motifs(seq, regions)

        # Persist in session state for other pages
        st.session_state["seq"] = seq
        st.session_state["perps"] = perps
        st.session_state["regions"] = regions
        st.session_state["motif_results"] = motif_results
        st.session_state["params"] = {
            "window": window_size,
            "threshold_pct": threshold_pct,
            "min_region_len": min_region_len,
        }
        st.session_state["analysis_done"] = True

    # ── Summary metrics ───────────────────────────────────────────────────────
    perps   = st.session_state["perps"]
    regions = st.session_state["regions"]
    motif_results = st.session_state["motif_results"]

    st.divider()
    st.subheader("📋 Analysis Summary")

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Sequence length", f"{len(st.session_state['seq']):,} bp")
    c2.metric("Mean perplexity", f"{float(perps.mean()):.3f}")
    c3.metric(
        "Threshold",
        f"{float(np.percentile(perps, st.session_state['params']['threshold_pct'])):.3f}",
    )
    c4.metric("Low-perplexity regions", len(regions))
    c5.metric("Non-B motif hits", len(motif_results))

    if regions:
        region_df = pd.DataFrame(regions, columns=["Start", "End"])
        region_df["Length (nt)"] = region_df["End"] - region_df["Start"]
        st.dataframe(region_df, use_container_width=True, hide_index=True)
    else:
        st.warning("No low-perplexity regions were detected with the current parameters.")

    st.info(
        "Navigate to **Perplexity Profile** or **Motif Browser** in the sidebar to explore the results.",
        icon="👈",
    )
