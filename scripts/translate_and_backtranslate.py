# translate_and_backtranslate.py
# Helper functions for translating DNA to protein and back-translating protein to DNA

from codon_tables import codon_table, reverse_codon_table
import random


def translate_dna(dna_seq):
    """
    Translate a DNA sequence into a protein sequence.
    Translation stops at the first stop codon.
    """
    protein_seq = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        aa = codon_table.get(codon, 'X')
        if aa == '*':
            break
        protein_seq += aa
    return protein_seq


def back_translate(protein_seq):
    """
    Generate a DNA sequence from a protein sequence
    using randomly selected synonymous codons.
    """
    return ''.join(
        [random.choice(reverse_codon_table[aa]) for aa in protein_seq if aa in reverse_codon_table]
    )


if __name__ == "__main__":
    test_dna = "ATGGCTAGTCCC"
    print("Translated:", translate_dna(test_dna))

    test_protein = "MAST"
    print("Back-translated:", back_translate(test_protein))
