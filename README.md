# Non-B DNA Perplexity Analyzer

A high-performance, multi-page **Streamlit** web application for detecting
**Non-B DNA** structures in genomic sequences using nucleotide perplexity and
regex-based motif scanning.

---

## ✨ Features

| Feature | Details |
|---------|---------|
| ⚡ Fast perplexity | Vectorised O(N) calculation with NumPy prefix sums |
| 📂 FASTA support | Single / multi-FASTA and plain-text files |
| 🔍 9 Non-B motifs | G4, i-Motif, Z-DNA, Triplex, STR, PolyA/T, Direct/Inverted/Mirror Repeats |
| 📊 Interactive charts | Perplexity landscape, region table, motif distribution bar chart |
| 🗂️ Filterable results | Filter by motif type, region length, motif length |
| ⬇️ CSV export | Download filtered or full result sets |
| 🧭 5-page navigation | Home · Analysis · Perplexity Profile · Motif Browser · About |

---

## 🚀 Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Run the app

```bash
streamlit run app.py
```

Then open <http://localhost:8501> in your browser.

---

## 📁 Project Structure

```
perplexity-nbdna/
├── app.py                          # Home page (landing / intro)
├── pages/
│   ├── 1_🧬_Analysis.py           # Upload & run analysis
│   ├── 2_📊_Perplexity_Profile.py # Perplexity chart & region table
│   ├── 3_🔬_Motif_Browser.py      # Filterable motif results + CSV
│   └── 4_ℹ️_About.py             # Methodology, references
├── requirements.txt
└── README.md
```

---

## 🧮 Algorithm

### Perplexity (O(N) vectorised)

For each overlapping window of length *W*, the nucleotide perplexity is:

```
Perplexity = 2^H   where   H = -Σ p_b log2(p_b)
```

* **Perplexity = 4** → perfectly uniform base composition (high complexity)
* **Perplexity → 1** → single nucleotide dominates (low complexity / Non-B prone)

Counts are computed via **cumulative-sum (prefix-sum)** arrays so the inner
loop over window characters is eliminated, giving **O(N)** time instead of
the naïve O(N·W).

### Low-perplexity region detection

1. Compute a threshold τ at the *p*-th percentile (default p = 5).
2. Assign scores `s_i = τ − Perplexity_i` (positive where low-complexity).
3. Run a **Kadane-variant sweep** to find contiguous stretches of positive
   score.
4. Merge overlapping intervals and discard regions shorter than `min_len`
   (default 100 nt).

### Non-B motif scanning

Nine compiled regex patterns are applied to the full sequence; only hits
that overlap a detected low-perplexity region are reported.

---

## ⚙️ Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Sliding window | 10 nt | Window size for per-position perplexity |
| Threshold percentile | 5 | Percentile of perplexity used as low-complexity cut-off |
| Minimum region length | 100 nt | Shortest low-perplexity region to retain |

All parameters are adjustable via the sidebar in the **Analysis** page.

---

## 📦 Dependencies

| Package | Minimum version | Purpose |
|---------|----------------|---------|
| `streamlit` | 1.34 | Web UI framework |
| `pandas` | 2.0 | Tabular data handling |
| `numpy` | 1.26 | Vectorised perplexity computation |

---

## 📄 License

This project is released under the [MIT License](LICENSE).
