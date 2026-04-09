import streamlit as st
import pandas as pd
import numpy as np
import re
from collections import Counter

st.title("Low Perplexity Non-B DNA Detector (Optimized)")

uploaded_file = st.file_uploader(
    "Upload FASTA or TXT sequence",
    type=["txt","fa","fasta"]
)

# ---------------- READ SEQUENCE (Optimized) ----------------
def read_sequence(uploaded_file):
    text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    
    # Faster sequence extraction using list comprehension
    seq_lines = [line.strip().upper() for line in text.splitlines() if not line.startswith(">")]
    sequence = "".join(seq_lines)
    sequence = re.sub("[^ACGT]", "", sequence)
    
    return sequence


# ---------------- PERPLEXITY (Vectorized) ----------------
def calc_perplexity_vectorized(seq_array):
    """Vectorized perplexity calculation for multiple sequences"""
    results = np.zeros(len(seq_array))
    
    for i, seq in enumerate(seq_array):
        # Use Counter for faster counting
        counts = Counter(seq)
        total = len(seq)
        
        if total == 0:
            results[i] = 0
            continue
        
        # Calculate entropy efficiently
        probs = [counts[n]/total for n in "ACGT" if counts[n] > 0]
        entropy = -sum(p * np.log2(p) for p in probs)
        results[i] = 2 ** entropy
    
    return results


# ---------------- SLIDING WINDOWS (Optimized with numpy) ----------------
def sliding_windows_fast(seq, window=10):
    """Fast sliding window using string slicing"""
    if len(seq) < window:
        return [], np.array([])
    
    # Pre-allocate windows list
    n_windows = len(seq) - window + 1
    windows = []
    sequences = []
    
    # Extract all windows in one pass
    for i in range(n_windows):
        sub = seq[i:i+window]
        windows.append((i, i+window))
        sequences.append(sub)
    
    # Calculate perplexities in batch
    perplexities = calc_perplexity_vectorized(sequences)
    
    return windows, perplexities


# ---------------- PERCENTILE (Optimized with numpy) ----------------
def percentile(values, percent):
    if not values:
        return None
    return np.percentile(values, percent)


# ---------------- MERGE REGIONS (Optimized) ----------------
def merge_regions(regions):
    if not regions:
        return []
    
    # Use numpy for sorting if possible
    regions = sorted(regions)
    merged = [list(regions[0])]
    
    for s, e in regions[1:]:
        last = merged[-1]
        if s <= last[1]:
            last[1] = max(last[1], e)
        else:
            merged.append([s, e])
    
    return merged


# ---------------- KADANE LOW PERPLEXITY DETECTION (Optimized) ----------------
def kadane_low_perplexity(windows, perps, threshold):
    scores = threshold - perps
    
    regions = []
    current_sum = 0
    start = 0
    
    for i, score in enumerate(scores):
        if current_sum <= 0:
            start = i
            current_sum = score
        else:
            current_sum += score
        
        if current_sum > 0:
            w_start = windows[start][0]
            w_end = windows[i][1]
            regions.append((w_start, w_end))
    
    return merge_regions(regions)


# ---------------- FILTER MINIMUM REGION LENGTH ----------------
def filter_min_length(regions, min_len=100):
    return [(s, e) for s, e in regions if (e - s) >= min_len]


# ---------------- MOTIFS (Pre-compiled, Optimized) ----------------
def build_regex():
    motifs = {
        "PolyA/T": re.compile(r"A{7,}|T{7,}"),
        "STR": re.compile(r"([ACGT]{1,6})\1{3,}"),
        "DirectRepeat_DR": re.compile(r"([ACGT]{4,10})[ACGT]{0,10}\1"),
        "InvertedRepeat_IR": re.compile(r"([ACGT]{4,})[ACGT]{0,10}\1"),
        "MirrorRepeat_MR": re.compile(r"([ACGT]{4,})[ACGT]{0,10}\1"),
        "G4": re.compile(r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"),
        "iMotif": re.compile(r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}"),
        "ZDNA": re.compile(r"(CG){4,}|(GC){4,}"),
        "Triplex_HDNA": re.compile(r"[AG]{10,}")
    }
    return motifs


# ---------------- OPTIMIZED MOTIF SEARCH ----------------
def find_motifs_in_regions(seq, motifs, regions):
    """Search for motifs only within low perplexity regions"""
    results = []
    
    # Pre-convert regions to list of tuples if not already
    region_list = list(regions)
    
    # For each motif, search in the sequence
    for name, pattern in motifs.items():
        # Find all motif occurrences
        for match in pattern.finditer(seq):
            ms, me = match.start(), match.end()
            motif_seq = match.group()
            
            # Check which regions this motif overlaps with
            for rs, re_ in region_list:
                if max(ms, rs) < min(me, re_):
                    results.append({
                        "Motif": name,
                        "Motif_Start": ms,
                        "Motif_End": me,
                        "Motif_Sequence": motif_seq,
                        "LowP_Start": rs,
                        "LowP_End": re_,
                        "Region_Length": re_ - rs
                    })
    
    return results


# ---------------- CACHED COMPUTATIONS ----------------
@st.cache_data
def process_sequence(seq, window_size=10, min_region_len=100, percentile_thresh=5):
    """Cached main processing function"""
    # Sliding windows
    windows, perps = sliding_windows_fast(seq, window_size)
    
    if not perps:
        return None
    
    # Threshold
    threshold = percentile(perps, percentile_thresh)
    
    # Kadane detection
    regions = kadane_low_perplexity(windows, perps, threshold)
    
    # Filter by length
    merged = filter_min_length(regions, min_region_len)
    
    # Motifs
    motifs = build_regex()
    
    # Find motifs
    results = find_motifs_in_regions(seq, motifs, merged)
    
    return results


# ---------------- MAIN ----------------
if uploaded_file:
    seq = read_sequence(uploaded_file)
    
    if len(seq) == 0:
        st.error("No valid DNA sequence found in file.")
        st.stop()
    
    if len(seq) < 100:
        st.error("Sequence must be at least 100 bp.")
        st.stop()
    
    # Add progress indicator
    with st.spinner("Processing sequence... This may take a moment for large sequences."):
        # Process with caching
        results = process_sequence(seq, window_size=10, min_region_len=100, percentile_thresh=5)
    
    if results is None or not results:
        st.warning("No motifs detected in low perplexity regions.")
    else:
        df = pd.DataFrame(results)
        
        # Display statistics
        st.success(f"Found {len(df)} motif-region overlaps")
        
        # Show summary stats
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Unique Motifs", df['Motif'].nunique())
        with col2:
            st.metric("Total Regions", df['LowP_Start'].nunique())
        with col3:
            st.metric("Avg Region Length", f"{df['Region_Length'].mean():.0f} bp")
        
        # Display dataframe
        st.dataframe(df)
        
        # Download button
        st.download_button(
            "Download CSV",
            df.to_csv(index=False),
            "low_perplexity_nonB_results.csv",
            "text/csv"
        )
        
        # Optional: Show motif distribution
        if st.checkbox("Show motif distribution"):
            motif_counts = df['Motif'].value_counts()
            st.bar_chart(motif_counts)
