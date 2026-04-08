# file: low_perplexity_nonB_folder_scan.py

import math
import re
import csv
import os


def read_sequence(file):

    seq = []

    with open(file) as f:

        for line in f:

            if line.startswith(">"):
                continue

            seq.append(line.strip().upper())

    sequence = "".join(seq)

    sequence = re.sub("[^ACGT]", "", sequence)

    return sequence


def calculate_perplexity(seq):

    counts = {n: seq.count(n) for n in "ACGT"}

    total = sum(counts.values())

    probabilities = []

    for c in counts.values():

        if c > 0:
            probabilities.append(c / total)

    entropy = -sum(p * math.log2(p) for p in probabilities)

    perplexity = 2 ** entropy

    return perplexity


def sliding_windows(seq, window=100):

    windows = []
    perplexities = []

    for i in range(len(seq) - window + 1):

        sub = seq[i:i + window]

        p = calculate_perplexity(sub)

        windows.append((i, i + window, p))
        perplexities.append(p)

    return windows, perplexities


def percentile(values, percent):

    values = sorted(values)

    index = int(len(values) * percent / 100)

    return values[index]


def bottom_percentile_windows(windows, perplexities, percent=5):

    threshold = percentile(perplexities, percent)

    regions = []

    for start, end, p in windows:

        if p <= threshold:
            regions.append((start, end))

    return regions, threshold


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


def build_nonb_regex():

    motifs = {

        "polyA_polyT": r"A{7,}|T{7,}",

        "STR_repeat": r"([ACGT]{1,6})\1{4,}",

        "G_quadruplex":
        r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",

        "i_motif":
        r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",

        "Z_DNA":
        r"(CG){4,}|(GC){4,}",

        "eGZ":
        r"(CGG){4,}|(GGC){4,}",

        "Sticky_triplex":
        r"(GAA){50,}|(TTC){50,}"
    }

    return {k: re.compile(v) for k, v in motifs.items()}


def scan_motifs(seq, regex_dict):

    hits = []

    for name, regex in regex_dict.items():

        for m in regex.finditer(seq):

            hits.append((name, m.start(), m.end(), m.group()))

    return hits


def overlap(a_start, a_end, b_start, b_end):

    return max(a_start, b_start) < min(a_end, b_end)


def intersect_motifs_lowP(motifs, regions):

    results = []

    for name, ms, me, seq in motifs:

        for rs, re in regions:

            if overlap(ms, me, rs, re):

                results.append((name, ms, me, seq, rs, re))

    return results


def write_csv(results, output_file):

    with open(output_file, "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow([
            "Motif",
            "Motif_start",
            "Motif_end",
            "Motif_sequence",
            "LowP_start",
            "LowP_end"
        ])

        for r in results:
            writer.writerow(r)


def process_file(file_path):

    print("Processing:", file_path)

    seq = read_sequence(file_path)

    windows, perplexities = sliding_windows(seq, 100)

    low_regions, threshold = bottom_percentile_windows(
        windows,
        perplexities,
        percent=5
    )

    merged_regions = merge_regions(low_regions)

    regex_dict = build_nonb_regex()

    motifs = scan_motifs(seq, regex_dict)

    overlaps = intersect_motifs_lowP(motifs, merged_regions)

    output = os.path.splitext(file_path)[0] + "_results.csv"

    write_csv(overlaps, output)

    print("Output:", output)
    print("Perplexity threshold:", threshold)
    print()


# -------------------------
# SET YOUR FOLDER PATH HERE
# -------------------------

folder = r"C:\Users\user\Desktop\OT"

for file in os.listdir(folder):

    if file.endswith((".txt", ".fa", ".fasta")):

        full_path = os.path.join(folder, file)

        process_file(full_path)
