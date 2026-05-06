import os
from collections import defaultdict, Counter
from Bio import SeqIO
import pandas as pd

# ---------------- Codon Table ---------------- #
SynonymousCodons = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'],
    'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG']
}

bases = ["A", "C", "G", "T"]
STOP_CODONS = {"TAA", "TAG", "TGA"}

def is_stop_codon(codon):
    return codon in STOP_CODONS

# --- Mutation Types in Required Order --- #
mutation_types = [
    "A→G", "C→T", "G→A", "T→C",
    "A→C", "A→T", "C→A", "C→G",
    "G→C", "G→T", "T→A", "T→G"
]

# --- Fixed Codon Order (U instead of T for display) ---
codon_order = [
    "UUU","UUC","UUA","UUG","CUU","CUC","CUA","CUG",
    "AUU","AUC","AUA","AUG","GUU","GUC","GUA","GUG",
    "UCU","UCC","UCA","UCG","CCU","CCC","CCA","CCG",
    "ACU","ACC","ACA","ACG","GCU","GCC","GCA","GCG",
    "UAU","UAC","UAA","UAG","CAU","CAC","CAA","CAG",
    "AAU","AAC","AAA","AAG","GAU","GAC","GAA","GAG",
    "UGU","UGC","UGG","UGA","CGU","CGC","CGA","CGG","AGU",
    "AGC","AGA","AGG","GGU","GGC","GGA","GGG"
]

# --- Degeneracy classification --- #
four_fold_codons = {
    "CCU","CCC","CCA","CCG",  # Proline
    "GCU","GCC","GCA","GCG",  # Alanine
    "ACU","ACC","ACA","ACG",  # Threonine
    "GGU","GGC","GGA","GGG",  # Glycine
    "GUU","GUC","GUA","GUG"   # Valine
}

two_fold_codons = {
    "UUU","UUC",    # Phe
    "UAU","UAC",    # Tyr
    "CAU","CAC",    # His
    "CAA","CAG",    # Gln
    "AAU","AAC",    # Asn
    "AAA","AAG",    # Lys
    "GAU","GAC",    # Asp
    "GAA","GAG",    # Glu
    "UGU","UGC"     # Cys
}

# ---------------- Helper Functions ---------------- #
def translate_codon(codon):
    for aa, codons in SynonymousCodons.items():
        if codon in codons:
            return aa
    return None


def calculate_codon(sequence):
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    counts = defaultdict(int)
    for codon in codons:
        if len(codon) == 3:
            counts[codon] += 1
    return counts


def generate_expected_dicts(CDict):
    Synonym_exp = {codon: {mut: 0 for mut in mutation_types} for codon in CDict}
    NonSynm_exp = {codon: {mut: 0 for mut in mutation_types} for codon in CDict}

    for codon, count in CDict.items():
        if count == 0:
            continue
        original_aa = translate_codon(codon)
        for pos in range(3):
            for b in bases:
                if b == codon[pos]:
                    continue
                mutated_codon = codon[:pos] + b + codon[pos+1:]

                # Exclude mutations that produce stop codons
                if is_stop_codon(mutated_codon):
                    continue

                mutated_aa = translate_codon(mutated_codon)
                mutation = f"{codon[pos]}→{b}"

                if mutated_aa is None:
                    continue
                if mutated_aa == original_aa:
                    Synonym_exp[codon][mutation] += count
                else:
                    NonSynm_exp[codon][mutation] += count
    return Synonym_exp, NonSynm_exp


def sum_mutations(exp_dict):
    """Sum all codon-level counts for each mutation type."""
    totals = {mut: 0 for mut in mutation_types}
    for codon, mut_dict in exp_dict.items():
        for mut, val in mut_dict.items():
            if mut in totals:
                totals[mut] += val
    return totals


# ---------------- Observed counting (per-nucleotide unique) ---------------- #
def get_observed_mutations_unique(fasta_file):
    # Read sequences but only keep 'clean' sequences (only A/C/G/T/U allowed)
    raw_sequences = []
    skipped = 0
    for rec in SeqIO.parse(fasta_file, "fasta"):
        s = str(rec.seq).upper()
        # allowed characters: A C G T U
        allowed = set("ACGTU")
        if set(s) <= allowed:
            # convert U to T for internal processing
            raw_sequences.append(s.replace("U", "T"))
        else:
            skipped += 1

    if not raw_sequences:
        raise ValueError(f"No clean sequences found in {fasta_file}. {skipped} sequence(s) skipped due to invalid characters.")

    # Trim to minimum length (multiple of 3)
    min_len = min(len(s) for s in raw_sequences)
    min_len = min_len - (min_len % 3)
    sequences = [s[:min_len] for s in raw_sequences]

    # Build reference sequence (majority rule)
    ref_nuc = []
    for pos in range(min_len):
        bases_here = [s[pos] for s in sequences if s[pos] in "ACGT"]
        if bases_here:
            ref_nuc.append(Counter(bases_here).most_common(1)[0][0])
        else:
            ref_nuc.append("N")
    ref_seq = "".join(ref_nuc)

    Synonym_obs = defaultdict(lambda: defaultdict(int))
    NonSynm_obs = defaultdict(lambda: defaultdict(int))
    Observed_total = defaultdict(int)

    seen = set()  # keep track of unique codon mutations (per position + ref_base + alt)

    for i in range(0, min_len, 3):
        codon_ref = ref_seq[i:i+3]
        if len(codon_ref) < 3 or any(b not in "ACGT" for b in codon_ref):
            continue
        original_aa = translate_codon(codon_ref)

        for s in sequences:
            codon_alt = s[i:i+3]
            if len(codon_alt) < 3 or any(b not in "ACGT" for b in codon_alt):
                continue

            diffs = [(pos, codon_ref[pos], codon_alt[pos]) for pos in range(3) if codon_ref[pos] != codon_alt[pos]]

            # only allow single-nucleotide codon mutations per sequence (we count them unique across sequences)
            if len(diffs) == 1:
                offset, ref_base, alt = diffs[0]
                mut_key = f"{ref_base}→{alt}"
                key = (i, codon_ref, ref_base, alt)
                if key not in seen:
                    seen.add(key)
                    Observed_total[mut_key] += 1
                    mutated_codon = codon_ref[:offset] + alt + codon_ref[offset+1:]

                    # Exclude observed mutations forming stop codons
                    if is_stop_codon(mutated_codon):
                        continue

                    mutated_aa = translate_codon(mutated_codon)

                    if mutated_aa is None or original_aa is None:
                        NonSynm_obs[codon_ref][mut_key] += 1
                    else:
                        if mutated_aa == original_aa:
                            Synonym_obs[codon_ref][mut_key] += 1
                        else:
                            NonSynm_obs[codon_ref][mut_key] += 1
            # else ignore multi-nucleotide differences in a single codon for observed unique counting

    print(f"[{os.path.basename(fasta_file)}] Read {len(raw_sequences)+skipped} sequences: {len(raw_sequences)} clean used, {skipped} skipped.")
    return ref_seq, Synonym_obs, NonSynm_obs, Observed_total, sequences


# ---------------- Mutation Position Finder ---------------- #
# --- Transition vs Transversion check ---
def mutation_type(base_from, base_to):
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    if (base_from, base_to) in transitions:
        return "ti"  # transition
    else:
        return "tv"  # transversion


def get_mutation_positions(ref_seq, sequences):
    mutation_positions = []
    min_len = len(ref_seq)

    for i in range(0, min_len, 3):
        codon_ref = ref_seq[i:i+3]
        if len(codon_ref) < 3 or any(b not in "ACGT" for b in codon_ref):
            continue
        ref_AA = translate_codon(codon_ref)

        # For this codon, collect offsets observed across sequences, but only from sequences
        # that have a single-nucleotide difference in this codon. Ignore sequences with
        # multi-nucleotide differences within this codon.
        offset_to_alts = defaultdict(set)
        valid_single_diff_sequences = 0

        for s in sequences:
            codon_alt = s[i:i+3]
            if len(codon_alt) < 3 or any(b not in "ACGT" for b in codon_alt):
                continue
            diffs = [(pos, codon_ref[pos], codon_alt[pos]) for pos in range(3) if codon_ref[pos] != codon_alt[pos]]

            if len(diffs) == 1:
                valid_single_diff_sequences += 1
                offset, ref_base, alt = diffs[0]
                offset_to_alts[offset].add((ref_base, alt, codon_alt))

        # If mutations are observed at more than one offset within this codon, skip the codon entirely
        if len(offset_to_alts) != 1:
            continue

        # Only one offset present; record each alternative base observed at this offset
        only_offset = next(iter(offset_to_alts.keys()))
        ref_base = codon_ref[only_offset]
        nucleotide_pos = i + only_offset + 1  # 1-based

        # For uniqueness, collect alt bases
        alts = {alt for (_, alt, _) in offset_to_alts[only_offset]}

        for alt in alts:
            mut_codon = codon_ref[:only_offset] + alt + codon_ref[only_offset+1:]

            # Skip stop codon mutations
            if is_stop_codon(mut_codon):
                continue

            mut_AA = translate_codon(mut_codon)
            mut_class = mutation_type(ref_base, alt)
            if ref_AA == mut_AA:
                label = "Syn" + ("ti" if mut_class == "ti" else "tv")
            else:
                label = "NSyn" + ("ti" if mut_class == "ti" else "tv")

            mutation_positions.append({
                "Codon_range": f"{i+1}-{i+3}",
                "Nucleotide_pos": nucleotide_pos,
                "Reference_codon": codon_ref,
                "Reference_base": ref_base,
                "Mutations": alt,
                "Ref_AA": ref_AA,
                "Mut_AA": mut_AA,
                "Class": label
            })

    return mutation_positions


# ---------------- Processing ---------------- #
def process(ref_seq, Synonym_exp, NonSynm_exp, Synonym_obs, NonSynm_obs, Observed_total, CDict):
    total_Synonym_exp = sum_mutations(Synonym_exp)
    total_NonSynm_exp = sum_mutations(NonSynm_exp)
    total_Synonym_obs = sum_mutations(Synonym_obs)
    total_NonSynm_obs = sum_mutations(NonSynm_obs)

    ratio_syn = {}
    ratio_nonsyn = {}
    for mut in mutation_types:
        exp_syn = total_Synonym_exp.get(mut, 0)
        obs_syn = total_Synonym_obs.get(mut, 0)
        exp_nonsyn = total_NonSynm_exp.get(mut, 0)
        obs_nonsyn = total_NonSynm_obs.get(mut, 0)

        ratio_syn[mut] = (obs_syn / exp_syn) if exp_syn > 0 else 0
        ratio_nonsyn[mut] = (obs_nonsyn / exp_nonsyn) if exp_nonsyn > 0 else 0

    def dict_to_df(codon_dict, CDict):
        # Create dataframe from mutation dictionary
        df = pd.DataFrame.from_dict(codon_dict, orient="index").fillna(0)

        # Ensure all mutation columns exist
        df = df.reindex(columns=mutation_types, fill_value=0)

        # Convert index to RNA format
        df.index = [c.replace("T", "U") for c in df.index]

        # Reindex to full codon list
        df = df.reindex(codon_order, fill_value=0)

        # Insert correct Codon_Count from reference dictionary
        df.insert(0, "Codon_Count", [CDict.get(c.replace("U", "T"), 0) for c in df.index])
        
        # --- Add TOTAL row --- #
        total_row = df.sum(numeric_only=True)
        total_row.name = "TOTAL"

        df = pd.concat([df, total_row.to_frame().T])

        return df

    df_syn_exp = dict_to_df(Synonym_exp, CDict)
    df_nonsyn_exp = dict_to_df(NonSynm_exp, CDict)
    df_syn_obs = dict_to_df(Synonym_obs, CDict)
    df_nonsyn_obs = dict_to_df(NonSynm_obs, CDict)

    df_ratios = pd.DataFrame({
        "Total_Synonym_exp": [total_Synonym_exp.get(m, 0) for m in mutation_types],
        "Total_Synonym_obs": [total_Synonym_obs.get(m, 0) for m in mutation_types],
        "Synonym_ratio": [ratio_syn[m] for m in mutation_types],
        "Total_NonSynm_exp": [total_NonSynm_exp.get(m, 0) for m in mutation_types],
        "Total_NonSynm_obs": [total_NonSynm_obs.get(m, 0) for m in mutation_types],
        "NonSynm_ratio": [ratio_nonsyn[m] for m in mutation_types]
    }, index=mutation_types)

    df_observed_total = pd.DataFrame({
        "Observed_Total": [Observed_total.get(m, 0) for m in mutation_types]
    }, index=mutation_types)

    return (df_syn_exp, df_nonsyn_exp, df_syn_obs, df_nonsyn_obs,
            df_ratios, df_observed_total,
            total_Synonym_exp, total_NonSynm_exp,
            total_Synonym_obs, total_NonSynm_obs)


def make_summary_sheet(total_Synonym_exp, total_NonSynm_exp,
                       total_Synonym_obs, total_NonSynm_obs,
                       CDict):

    # --- Expected sums --- #
    exp_syn_ti = sum(val for mut, val in total_Synonym_exp.items()
                     if mutation_type(mut[0], mut[-1]) == "ti")
    exp_syn_tv = sum(val for mut, val in total_Synonym_exp.items()
                     if mutation_type(mut[0], mut[-1]) == "tv")
    exp_nonsyn_ti = sum(val for mut, val in total_NonSynm_exp.items()
                        if mutation_type(mut[0], mut[-1]) == "ti")
    exp_nonsyn_tv = sum(val for mut, val in total_NonSynm_exp.items()
                        if mutation_type(mut[0], mut[-1]) == "tv")

    # --- Observed sums --- #
    obs_syn_ti = sum(val for mut, val in total_Synonym_obs.items()
                     if mutation_type(mut[0], mut[-1]) == "ti")
    obs_syn_tv = sum(val for mut, val in total_Synonym_obs.items()
                     if mutation_type(mut[0], mut[-1]) == "tv")
    obs_nonsyn_ti = sum(val for mut, val in total_NonSynm_obs.items()
                        if mutation_type(mut[0], mut[-1]) == "ti")
    obs_nonsyn_tv = sum(val for mut, val in total_NonSynm_obs.items()
                        if mutation_type(mut[0], mut[-1]) == "tv")

    # --- Obs/Exp ratios for NON-synonymous --- #
    nonsyn_ti_ratio = obs_nonsyn_ti / exp_nonsyn_ti if exp_nonsyn_ti > 0 else 0
    nonsyn_tv_ratio = obs_nonsyn_tv / exp_nonsyn_tv if exp_nonsyn_tv > 0 else 0
    
    # --- nti and ntv --- #
    nti = obs_syn_ti + obs_nonsyn_ti
    ntv = obs_syn_tv + obs_nonsyn_tv
    nti_ntv_ratio = nti / ntv if ntv > 0 else 0
    
    
    # --- Ratios --- #
    exp_syn_ratio = exp_syn_ti / exp_syn_tv if exp_syn_tv > 0 else 0
    exp_nonsyn_ratio = exp_nonsyn_ti / exp_nonsyn_tv if exp_nonsyn_tv > 0 else 0
    exp_ti_tv = ((exp_syn_ti + exp_nonsyn_ti) /
                 (exp_syn_tv + exp_nonsyn_tv)
                 if (exp_syn_tv + exp_nonsyn_tv) > 0 else 0)

    obs_syn_ratio = obs_syn_ti / obs_syn_tv if obs_syn_tv > 0 else 0
    obs_nonsyn_ratio = obs_nonsyn_ti / obs_nonsyn_tv if obs_nonsyn_tv > 0 else 0
    obs_ti_tv = ((obs_syn_ti + obs_nonsyn_ti) /
                 (obs_syn_tv + obs_nonsyn_tv)
                 if (obs_syn_tv + obs_nonsyn_tv) > 0 else 0)

    # --- Degeneracy counts (from reference only) --- #
    four_count = sum(
        CDict.get(c.replace("U", "T"), 0)
        for c in four_fold_codons
    )

    two_count = sum(
        CDict.get(c.replace("U", "T"), 0)
        for c in two_fold_codons
    )

    total_codons = sum(CDict.values())
    
    # --- Degeneracy percentages --- #
    FFD_percent = (four_count / total_codons) * 100 if total_codons > 0 else 0
    TFD_percent = (two_count / total_codons) * 100 if total_codons > 0 else 0
    TFD_FFD_ratio = two_count / four_count if four_count > 0 else 0

    data = {
        
        "Total_Exp_Syn_ti": [exp_syn_ti],
        "Total_Exp_Syn_tv": [exp_syn_tv],
        "Total_Exp_NonSyn_ti": [exp_nonsyn_ti],
        "Total_Exp_NonSyn_tv": [exp_nonsyn_tv],

        "Exp_Syn_ti/tv": [exp_syn_ratio],
        "Exp_NonSyn_ti/tv": [exp_nonsyn_ratio],
        "Exp_ti/tv": [exp_ti_tv],

        "Total_Obs_Syn_ti": [obs_syn_ti],
        "Total_Obs_Syn_tv": [obs_syn_tv],
        "Total_Obs_NonSyn_ti": [obs_nonsyn_ti],
        "Total_Obs_NonSyn_tv": [obs_nonsyn_tv],

        "Obs_Syn_ti/tv": [obs_syn_ratio],
        "Obs_NonSyn_ti/tv": [obs_nonsyn_ratio],
        "Obs_ti/tv": [obs_ti_tv],
        
              
        "NSyn_ti_obs/exp": [nonsyn_ti_ratio],
        "NSyn_tv_obs/exp": [nonsyn_tv_ratio],

        "nti": [nti],
        "ntv": [ntv],
        "nti/ntv": [nti_ntv_ratio],
        
        
        "Total_Codons": [total_codons],

        "Four_fold_count": [four_count],
        "Two_fold_count": [two_count],
        
        "FFD%": [FFD_percent],
        "TFD%": [TFD_percent],
        "TFD:FFD": [TFD_FFD_ratio],

    }

    return pd.DataFrame(data)


def process_multiple_genes(fasta_files):
    for fasta_file in fasta_files:
        gene_name = os.path.splitext(os.path.basename(fasta_file))[0]
        ref_seq, Synonym_obs, NonSynm_obs, Observed_total, sequences = get_observed_mutations_unique(fasta_file)
        CDict = calculate_codon(ref_seq)
        Synonym_exp, NonSynm_exp = generate_expected_dicts(CDict)

        (df_syn_exp, df_nonsyn_exp, df_syn_obs, df_nonsyn_obs,
         df_ratios, df_observed_total,
         total_Synonym_exp, total_NonSynm_exp,
         total_Synonym_obs, total_NonSynm_obs) = process(
            ref_seq, Synonym_exp, NonSynm_exp, Synonym_obs, NonSynm_obs,
            Observed_total, CDict
        )

        df_summary = make_summary_sheet(total_Synonym_exp, total_NonSynm_exp,
                                        total_Synonym_obs, total_NonSynm_obs,
                                        CDict)

        mutation_positions = get_mutation_positions(ref_seq, sequences)
        df_mutpos = pd.DataFrame(mutation_positions)

        with pd.ExcelWriter(f"{gene_name}_mutation_analysis.xlsx", engine="openpyxl") as writer:
            df_syn_exp.to_excel(writer, sheet_name="Syn_exp")
            df_nonsyn_exp.to_excel(writer, sheet_name="NonSyn_exp")
            df_syn_obs.to_excel(writer, sheet_name="Syn_obs")
            df_nonsyn_obs.to_excel(writer, sheet_name="NonSyn_obs")
            df_ratios.to_excel(writer, sheet_name="Ratios")
            df_observed_total.to_excel(writer, sheet_name="Observed_Total")
            df_summary.to_excel(writer, sheet_name="Summary")
            df_mutpos.to_excel(writer, sheet_name="Mutation_positions", index=False)


# ---------------- Main ---------------- #
if __name__ == "__main__":
    fasta_files = [f for f in os.listdir() if f.endswith(".fasta") or f.endswith(".fa") or f.endswith(".txt") ]
    if not fasta_files:
        print("No FASTA files found in current directory.")
    else:
        process_multiple_genes(fasta_files)
        print("Analysis complete. Excel files generated.")
