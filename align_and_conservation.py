# align_and_conservation.py
# Manual MSA + consensus conservation profile from AtCAMTA protein sequences

from pathlib import Path
import pandas as pd
from itertools import zip_longest


def load_fasta(file_path):
    """Load protein sequences from a FASTA file."""
    entries = Path(file_path).read_text().strip().split("\n>")
    records = []
    for entry in entries:
        lines = entry.strip().splitlines()
        header = lines[0].replace('>', '') if lines[0].startswith('>') else lines[0]
        seq = ''.join(lines[1:]).replace(" ", "").replace("\n", "")
        records.append((header, seq))
    return records


def pad_sequences(records):
    """Naively align sequences by padding with '-' to the max length."""
    max_len = max(len(seq) for _, seq in records)
    return [(hdr, seq.ljust(max_len, '-')) for hdr, seq in records]


def compute_conservation(aligned_records):
    """Compute consensus amino acid and conservation per column."""
    aligned_seqs = [seq for _, seq in aligned_records]
    consensus_profile = []

    for i, column in enumerate(zip_longest(*aligned_seqs, fillvalue='-')):
        aa_counts = {aa: column.count(aa) for aa in set(column)}
        most_common = max(aa_counts, key=aa_counts.get)
        conservation = round(aa_counts[most_common] / len(column) * 100, 2)
        consensus_profile.append({
            "Position": i + 1,
            "Consensus_AA": most_common,
            "Conservation (%)": conservation,
            "Column": ''.join(column)
        })
    return pd.DataFrame(consensus_profile)


if __name__ == "__main__":
    input_fasta = "data/AtCAMTA-Prot.fasta"
    output_csv = "data/AtCAMTA_consensus_profile.csv"

    records = load_fasta(input_fasta)
    aligned_records = pad_sequences(records)
    consensus_df = compute_conservation(aligned_records)
    consensus_df.to_csv(output_csv, index=False)

    print(f"Consensus conservation profile saved to {output_csv}")
