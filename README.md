# 🌱 AtCAMTA-Synthetic

**AI-Guided Design of Synthetic CAMTA Transcription Factors in Plants**  
By [Adhityo Wicaksono](https://github.com/adhityow) | Code support by ChatGPT (OpenAI, May 2025)

---

## 🧬 Overview

This repository contains the scripts and logic used to generate **AtCAMTA7** and **AtCAMTA8**, two synthetic transcription factors derived from the *Arabidopsis thaliana* CAMTA family (AtCAMTA1–6). These synthetic genes are designed using a combination of:

- **Codon-level divergence** from AtCAMTA3 (AtCAMTA7)
- **Consensus- and AI-guided de novo synthesis** based on full AtCAMTA family alignment (AtCAMTA8)

The project demonstrates a bottom-up synthetic genome strategy, enabling gradual functional divergence and innovation using AI in plant molecular biology.

---

## 🧪 What’s Included

- `synthesize_CAMTA7.py` – Codon-reshuffling engine to generate AtCAMTA7
- `synthesize_CAMTA8.py` – De novo generator from multi-CAMTA consensus logic
- `codon_tables.py` – Reversible codon table for back-translation
- `translate_and_backtranslate.py` – Helpers to convert between CDS and protein
- `align_and_conservation.py` – (Optional) aligner for AtCAMTA1–6 proteins
- `outputs/` – FASTA files of AtCAMTA7 and AtCAMTA8 (protein and CDS)
- `data/` – Raw protein, CDS, and genomic sequences for AtCAMTA1–6

---

## 📊 Key Features

- **Functional Conservation Aware**: Preserves CG-1, ANK, and IQ domain structures
- **Synthetic Divergence Logic**: Position-by-position redesign based on conservation profile
- **Codon Customization**: Designed for plant chassis expression (e.g., *Nicotiana*, *Arabidopsis*)
- **AI + Domain Biology**: Inspired by Evo-2, ProGen, and CAMTA domain modularity

---

## 🛠 Requirements

- Python 3.8+
- `pandas`
- `Bio` (Biopython, optional)
- `random`, `collections`, `subprocess`

---

## 🚀 How to Run

```bash
# Synthesize AtCAMTA7
python synthesize_CAMTA7.py

# Synthesize AtCAMTA8 from consensus scaffold
python synthesize_CAMTA8.py

# Output: FASTA files for both CDS and protein
