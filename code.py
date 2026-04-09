import streamlit as st
import pandas as pd
import math
import re


# ---------------- READ SEQUENCE ----------------
def read_sequence(uploaded_file):

    text = uploaded_file.getvalue().decode("utf-8")

    lines = text.splitlines()

    seq = []

    for line in lines:

        if line.startswith(">"):
            continue

        seq.append(line.strip().upper())

    sequence = "".join(seq)

    sequence = re.sub("[^ACGT]", "", sequence)

    return sequence


# ---------------- PERPLEXITY ----------------
def calculate_perplexity(seq):

    counts = {n: seq.count(n) for n in "ACGT"}

    total = sum(counts.values())

    if total == 0:
        return 0

    probs = [c/total for c in counts.values() if c > 0]

    entropy = -sum(p * math.log2(p) for p in probs)

    return 2 ** entropy


# ---------------- SLIDING WINDOW ----------------
def sliding_windows(seq, window=100):

    windows = []
    perplexities = []

    if len(seq) < window:
        return windows, perplexities

    for i in range(len(seq) - window + 1):

        sub = seq[i:i+window]

        p = calculate_perplexity(sub)

        windows.append((i, i+window, p))

        perplexities.append(p)

    return windows, perplexities


# ---------------- PERCENTILE ----------------
def percentile(values, percent):

    if len(values) == 0:
        return 0

    values = sorted(values)

    index = int(len(values) * percent / 100)

    index = min(index, len(values)-1)

    return values[index]


# ---------------- LOW PERPLEXITY ----------------
def bottom_percentile_windows(windows, perplexities, percent=5):

    threshold = percentile(perplexities, percent)

    regions = []

    for start, end, p in windows:

        if p <= threshold:
            regions.append((start, end))

    return regions, threshold


# ---------------- MERGE REGIONS ----------------
def merge_regions(regions):

    if not regions:
        return []

    regions = sorted(regions)

    merged = [list(regions[0])]

    for s, e in regions[1:]:

        last = merged[-1]

        if s <= last[1]:

            last[1] = max(last[1], e)

        else:

            merged.append([s, e])

    return merged


# ---------------- MOTIFS ----------------
def build_nonb_regex():

    motifs = {

        "polyA_polyT": r"A{7,}|T{7,}",

        "STR_repeat": r"([ACGT]{1,6})\1{4,}",

        "G_quadruplex":
        r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",

        "i_motif":
        r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",

        "Z_DNA":
        r"(CG){4,}|(GC){4,}"
    }

    return {k: re.compile(v) for k,v in motifs.items()}


# ---------------- SCAN MOTIFS ----------------
def scan_motifs(seq, regex_dict):

    hits = []

    for name, regex in regex_dict.items():

        for m in regex.finditer(seq):

            hits.append((name, m.start(), m.end(), m.group()))

    return hits


# ---------------- OVERLAP ----------------
def overlap(a_start, a_end, b_start, b_end):

    return max(a_start, b_start) < min(a_end, b_end)


# ---------------- INTERSECTION ----------------
def intersect_motifs_lowP(motifs, regions):

    results = []

    for name, ms, me, motif_seq in motifs:

        for rs, re in regions:

            if overlap(ms, me, rs, re):

                results.append({

                    "Motif": name,
                    "Motif_Start": ms,
                    "Motif_End": me,
                    "Sequence": motif_seq,
                    "LowP_Start": rs,
                    "LowP_End": re
                })

    return results


# ---------------- STREAMLIT APP ----------------

st.title("Low Perplexity Non-B DNA Detector")

uploaded_file = st.file_uploader(
    "Upload FASTA or TXT sequence",
    type=["txt","fa","fasta"]
)


if uploaded_file:

    seq = read_sequence(uploaded_file)

    if len(seq) < 100:

        st.error("Sequence must be at least 100 bp")

    else:

        windows, perplexities = sliding_windows(seq,100)

        low_regions, threshold = bottom_percentile_windows(
            windows,
            perplexities,
            5
        )

        merged = merge_regions(low_regions)

        regex_dict = build_nonb_regex()

        motifs = scan_motifs(seq, regex_dict)

        overlaps = intersect_motifs_lowP(motifs, merged)

        st.write("Perplexity Threshold:", threshold)

        if overlaps:

            df = pd.DataFrame(overlaps)

            st.dataframe(df)

            csv = df.to_csv(index=False)

            st.download_button(
                "Download CSV",
                csv,
                "results.csv",
                "text/csv"
            )

        else:

            st.warning("No motifs detected in low-perplexity regions.")
