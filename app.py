import streamlit as st
import pandas as pd
import math
import re

st.title("Low Perplexity Non-B DNA Detector")

uploaded_file = st.file_uploader(
    "Upload FASTA or TXT file",
    type=["txt","fa","fasta"]
)

# -------- Read Sequence --------
def read_sequence(file):

    text = file.getvalue().decode("utf-8")

    seq = []

    for line in text.splitlines():

        if line.startswith(">"):
            continue

        seq.append(line.strip().upper())

    sequence = "".join(seq)

    sequence = re.sub("[^ACGT]", "", sequence)

    return sequence


# -------- Perplexity --------
def perplexity(seq):

    counts = {n: seq.count(n) for n in "ACGT"}

    total = sum(counts.values())

    if total == 0:
        return 0

    probs = [c/total for c in counts.values() if c > 0]

    entropy = -sum(p * math.log2(p) for p in probs)

    return 2 ** entropy


# -------- Sliding Window --------
def sliding_windows(seq, window=100):

    windows=[]
    perps=[]

    if len(seq) < window:
        return windows,perps

    for i in range(len(seq)-window+1):

        sub = seq[i:i+window]

        p = perplexity(sub)

        windows.append((i,i+window,p))

        perps.append(p)

    return windows,perps


# -------- Merge Regions --------
def merge_regions(regions):

    if not regions:
        return []

    regions = sorted(regions)

    merged=[list(regions[0])]

    for s,e in regions[1:]:

        last=merged[-1]

        if s<=last[1]:

            last[1]=max(last[1],e)

        else:

            merged.append([s,e])

    return merged


# -------- Motifs --------
def build_regex():

    motifs = {

        "PolyA_T": r"A{7,}|T{7,}",

        "STR": r"([ACGT]{1,6})\1{4,}",

        "G4": r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",

        "iMotif": r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",

        "ZDNA": r"(CG){4,}|(GC){4,}"
    }

    return {k:re.compile(v) for k,v in motifs.items()}


# -------- App Logic --------
if uploaded_file:

    seq = read_sequence(uploaded_file)

    if len(seq) < 100:

        st.error("Sequence must be at least 100 bp")

    else:

        windows,perps = sliding_windows(seq,100)

        if len(perps) == 0:

            st.error("Unable to calculate perplexity")

        else:

            threshold = sorted(perps)[int(len(perps)*0.05)]

            low_regions=[(s,e) for s,e,p in windows if p<=threshold]

            merged = merge_regions(low_regions)

            regex_dict = build_regex()

            results=[]

            for name,rgx in regex_dict.items():

                for m in rgx.finditer(seq):

                    ms,me = m.start(),m.end()

                    for rs,re in merged:

                        if max(ms,rs) < min(me,re):

                            results.append({

                                "Motif":name,
                                "Motif_Start":ms,
                                "Motif_End":me,
                                "Sequence":m.group(),
                                "LowP_Start":rs,
                                "LowP_End":re
                            })

            if len(results) == 0:

                st.warning("No motifs detected")

            else:

                df = pd.DataFrame(results)

                st.dataframe(df)

                st.download_button(
                    "Download CSV",
                    df.to_csv(index=False),
                    "results.csv"
                )
