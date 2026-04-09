import streamlit as st
import pandas as pd
import math
import re

st.title("Low Perplexity Non-B DNA Detector")

uploaded_file = st.file_uploader(
    "Upload FASTA or TXT sequence",
    type=["txt","fa","fasta"]
)

# ---------------- READ SEQUENCE ----------------
def read_sequence(uploaded_file):
    text = uploaded_file.getvalue().decode("utf-8", errors="ignore")
    seq_lines = []

    for line in text.splitlines():
        if line.startswith(">"):
            continue
        seq_lines.append(line.strip().upper())

    sequence = "".join(seq_lines)
    sequence = re.sub("[^ACGT]", "", sequence)

    return sequence


# ---------------- PERPLEXITY ----------------
def calc_perplexity(seq):
    counts = {n: seq.count(n) for n in "ACGT"}
    total = sum(counts.values())

    if total == 0:
        return 0

    probs = [c/total for c in counts.values() if c > 0]
    entropy = -sum(p * math.log2(p) for p in probs)

    return 2 ** entropy


# ---------------- SLIDING WINDOWS ----------------
def sliding_windows(seq, window=100):

    windows = []
    perplexities = []

    if len(seq) < window:
        return windows, perplexities

    for i in range(len(seq) - window + 1):
        sub = seq[i:i+window]
        p = calc_perplexity(sub)

        windows.append((i, i+window, p))
        perplexities.append(p)

    return windows, perplexities


# ---------------- PERCENTILE ----------------
def percentile(values, percent):

    if not values:
        return None

    values = sorted(values)

    index = int(len(values) * percent / 100)
    index = min(index, len(values)-1)

    return values[index]


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
def build_regex():

    motifs = {
        "PolyA/T": r"A{7,}|T{7,}",
        "STR": r"([ACGT]{1,6})\1{4,}",
        "G4": r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",
        "iMotif": r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",
        "ZDNA": r"(CG){4,}|(GC){4,}"
    }

    return {k: re.compile(v) for k, v in motifs.items()}


# ---------------- MAIN ----------------
if uploaded_file:

    seq = read_sequence(uploaded_file)

    if len(seq) == 0:
        st.error("No valid DNA sequence found in file.")
        st.stop()

    if len(seq) < 100:
        st.error("Sequence must be at least 100 bp.")
        st.stop()

    windows, perps = sliding_windows(seq, 100)

    if not perps:
        st.error("Perplexity calculation failed.")
        st.stop()

    threshold = percentile(perps, 5)

    low_regions = [(s, e) for s, e, p in windows if p <= threshold]

    merged = merge_regions(low_regions)

    regex_dict = build_regex()

    results = []

    for name, rgx in regex_dict.items():

        for m in rgx.finditer(seq):

            ms, me = m.start(), m.end()

            for rs, re_ in merged:

                if max(ms, rs) < min(me, re_):

                    results.append({
                        "Motif": name,
                        "Motif_Start": ms,
                        "Motif_End": me,
                        "Motif_Sequence": m.group(),
                        "LowP_Start": rs,
                        "LowP_End": re_
                    })

    if not results:
        st.warning("No motifs detected in low perplexity regions.")
    else:
        df = pd.DataFrame(results)

        st.dataframe(df)

        st.download_button(
            "Download CSV",
            df.to_csv(index=False),
            "low_perplexity_nonB_results.csv",
            "text/csv"
        )
