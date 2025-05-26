# synthesize_CAMTA7.py
# Author: Adhityo Wicaksono
# Synthetic generation of AtCAMTA7 by codon-level divergence from AtCAMTA3

import random

# Reverse codon table (AA -> list of synonymous codons)
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

# Forward codon table (simplified)
codon_table = {codon: aa for aa, codons in reverse_codon_table.items() for codon in codons}

def translate_dna(dna_seq):
    protein_seq = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        aa = codon_table.get(codon, 'X')
        if aa == '*': break
        protein_seq += aa
    return protein_seq

def synthesize_camta7(original_dna):
    aa_seq = translate_dna(original_dna)
    synthetic_cds = ''.join([random.choice(reverse_codon_table[aa]) for aa in aa_seq])
    return aa_seq, synthetic_cds

if __name__ == "__main__":
    # Example input: AtCAMTA3 CDS
    from pathlib import Path

    input_path = Path("data/AtCAMTA-CDS.fasta")
    output_path = Path("outputs/AtCAMTA7_synthetic")
    output_path.mkdir(parents=True, exist_ok=True)

    from Bio import SeqIO
    records = list(SeqIO.parse(input_path, "fasta"))
    camta3 = [r for r in records if "At2g22300" in r.id][0]

    aa_seq, syn_cds = synthesize_camta7(str(camta3.seq))

    # Write outputs
    with open(output_path / "AtCAMTA7_protein.fasta", "w") as f:
        f.write(f">AtCAMTA7_synthetic_protein\n{aa_seq}\n")

    with open(output_path / "AtCAMTA7_CDS.fasta", "w") as f:
        f.write(f">AtCAMTA7_synthetic_CDS\n{syn_cds}\n")

    print("AtCAMTA7 synthesis complete. Files saved in outputs/AtCAMTA7_synthetic/")
