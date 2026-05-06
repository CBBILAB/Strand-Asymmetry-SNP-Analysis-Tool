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

**Example FASTA input:**
