# Strand-Asymmetry-SNP-Analysis-Tool

A Python script for analyzing synonymous single-nucleotide polymorphism (SNP) frequencies
in coding sequences (CDS) from multi-sequence FASTA alignments.

Developed for the study:
> **"Replication and gene expression have different impacts on synonymous nucleotide
> polymorphism contributing to strand inequality in the *Escherichia coli* chromosome"**
> Deka N. et al. (2025)

---

## What Does This Program Do?

This tool takes a set of aligned gene sequences (in FASTA format) from multiple
bacterial strains and calculates how often each of the 12 possible nucleotide
substitutions (e.g., C→T, G→A, A→T, etc.) occurs in coding regions.

It separates mutations into:
- **Synonymous** (silent) — mutations that do not change the amino acid
- **Non-synonymous** — mutations that change the amino acid

For each mutation type, it computes:
- **Expected counts** — how many synonymous/non-synonymous changes are theoretically
  possible given the codon composition of the reference sequence
- **Observed counts** — how many unique mutations are actually seen across all strains
- **Observed/Expected ratios** — to identify over- or under-represented mutation types
- **Transition/Transversion (ti/tv) ratios** — for each mutation type and overall
- **Codon degeneracy** — percentage of four-fold degenerate (FFD) and two-fold
  degenerate (TFD) codons in the reference
- **Mutation positions** — exact nucleotide position, codon, and amino acid change
  for every observed mutation, classified as synonymous/non-synonymous and
  transition/transversion

The results are saved as a multi-sheet **Excel file** (one file per gene/input FASTA).

---

## Who Should Use This Tool?

This tool is designed for **biologists and bioinformaticians** studying:
- Mutation bias in bacterial genomes
- Strand asymmetry in chromosomal SNPs
- Impact of gene expression on nucleotide substitution patterns
- Codon usage and degeneracy analysis

No advanced programming knowledge is required to run it.

---

## Input Format

- One or more **FASTA files** (`.fasta`, `.fa`, or `.txt` extension)
- Each FASTA file should contain **aligned sequences** of one gene from multiple strains
- All sequences must use only the bases: **A, C, G, T** (or U — it is handled automatically)
- Sequences with ambiguous characters (N, ?, etc.) are automatically skipped

For example a sample fasta file named sugE.txt may be used
---

## Output Format

One Excel file is generated per input FASTA, named `<gene_name>_mutation_analysis.xlsx`.

It contains these sheets:

| Sheet Name | Contents |
|---|---|
| `Syn_exp` | Expected synonymous mutation counts per codon and mutation type |
| `NonSyn_exp` | Expected non-synonymous mutation counts |
| `Syn_obs` | Observed synonymous mutation counts |
| `NonSyn_obs` | Observed non-synonymous mutation counts |
| `Ratios` | Observed/Expected ratios for all 12 mutation types |
| `Observed_Total` | Total observed mutations per type (syn + non-syn) |
| `Summary` | ti/tv counts, ratios, codon degeneracy (FFD%, TFD%) |
| `Mutation_positions` | Position, codon, amino acid, and class of every observed mutation |

The sample output file may be seen as sugE_mutation_analysis.xlsx
---

## How to Run the Program

### Step 1 — Install Python

Download Python 3.8 or newer from: https://www.python.org/downloads/

During installation on Windows, check ✅ **"Add Python to PATH"**.

### Step 2 — Install required packages

Open **Command Prompt** (Windows) or **Terminal** (Mac/Linux) and type:

```bash
pip install biopython pandas openpyxl
```

Wait for installation to complete.

### Step 3 — Prepare your files

Place your FASTA files (`.fasta`, `.fa`, or `.txt`) in a folder together with
the script `titv_Mar.py`.

### Step 4 — Run the program

Open Command Prompt / Terminal, navigate to your folder, and type:

```bash
python titv_Mar.py
```

The program will automatically find all FASTA files in the same folder and
analyze them one by one.

**Example (Windows):**

cd C:\Users\YourName\snp_analysis
python titv_Mar.py

### Step 5 — View results

After the program finishes, one `.xlsx` Excel file will appear in the same folder
for each FASTA file analyzed. Open it with Microsoft Excel or LibreOffice Calc.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `ModuleNotFoundError: No module named 'Bio'` | Run `pip install biopython` |
| `No FASTA files found` | Make sure your files end in `.fasta`, `.fa`, or `.txt` and are in the same folder as the script |
| Sequences with characters like N, ? are skipped | This is expected — only A/C/G/T/U sequences are used |
| Empty output sheet | Your gene may have fewer than 3 sequences; try with more strains |

---

## Citation

If you use this tool in your research, please cite:

> Deka N., et al. (2025). Replication and gene expression have different impacts on
> synonymous nucleotide polymorphism contributing to strand inequality in the
> *Escherichia coli* chromosome.

---

## License

MIT License — free to use, modify, and distribute with attribution.
