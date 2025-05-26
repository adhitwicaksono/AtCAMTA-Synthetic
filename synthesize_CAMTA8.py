# synthesize_CAMTA8.py
# Author: Adhityo Wicaksono
# Synthetic generation of AtCAMTA8 using MSA-derived consensus strategy

import random
import pandas as pd
from pathlib import Path

# Codon table for back-translation
reverse_codon_table = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA']
}

# Classification logic
THRESHOLDS = {
    'Conserved': 80.0,
    'Mutable': (50.0, 79.9),
    'Synthetic': 0.0
}

def classify(conservation):
    if conservation >= THRESHOLDS['Conserved']:
        return 'Conserved'
    elif THRESHOLDS['Mutable'][0] <= conservation <= THRESHOLDS['Mutable'][1]:
        return 'Mutable'
    else:
        return 'Synthetic'

def synthesize_from_consensus(consensus_df):
    aa_seq = ""
    for _, row in consensus_df.iterrows():
        state = classify(row['Conservation (%)'])
        if state == 'Conserved':
            aa_seq += row['Consensus_AA']
        elif state == 'Mutable':
            variants = list(set(row['Column']))
            variants = [v for v in variants if v != row['Consensus_AA']]
            aa_seq += random.choice(variants) if variants else row['Consensus_AA']
        else:  # Synthetic
            aa_seq += random.choice("ACDEFGHIKLMNPQRSTVWY")
    return aa_seq

def back_translate(protein_seq):
    return ''.join([random.choice(reverse_codon_table[aa]) for aa in protein_seq if aa in reverse_codon_table])

if __name__ == "__main__":
    input_csv = Path("data/AtCAMTA_consensus_profile.csv")
    output_path = Path("outputs/AtCAMTA8_synthetic")
    output_path.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)
    syn_protein = synthesize_from_consensus(df)
    syn_cds = back_translate(syn_protein)

    with open(output_path / "AtCAMTA8_protein.fasta", "w") as f:
        f.write(f">AtCAMTA8_synthetic_protein\n{syn_protein}\n")

    with open(output_path / "AtCAMTA8_CDS.fasta", "w") as f:
        f.write(f">AtCAMTA8_synthetic_CDS\n{syn_cds}\n")

    print("AtCAMTA8 synthesis complete. Files saved in outputs/AtCAMTA8_synthetic/")
