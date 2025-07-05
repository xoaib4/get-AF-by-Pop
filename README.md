# get-AF-by-Pop
Calculate allele frequency of all variants by population in a VCF file
# Allele Frequency Table by Population

This repository provides a Bash script that computes **allele frequencies (AF)** for every variant in a multi‑sample VCF file:

* **AF_Total** – allele frequency across *all* samples  
* **AF_<POP>** – allele frequency within each population, where the population code is extracted from each sample’s ID

The final result is `allele_freq_table.tsv`, a tab‑delimited table ready for downstream analysis or Excel.

---

## Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| **bash** | ≥ 3.0 | Shell interpreter (works on macOS/HPC) |
| **bcftools** | ≥ 1.9 | VCF processing, AF computation |

Your VCF must be **bgzipped** (`.vcf.gz`) and **indexed** (`.vcf.gz.tbi`) and must contain per‑sample **GT** genotypes.

---

## Sample Naming Convention

The script derives the population tag from **characters 8–10** of each sample ID (0‑based indexing).  
Example:

| Sample ID | Extracted population |
|-----------|----------------------|
| `HG_ADG_001` | `ADG` |
| `HG_BPG_022` | `BPG` |

You can change that rule by editing these two variables at the top of the script:

```bash
POP_TAG_START=7   # 0‑based position of first character
POP_TAG_LEN=3     # length of the population tag
How the Script Works
Read sample IDs from the VCF header.

Group sample IDs into population lists based on the tag.

Run bcftools +fill-tags to calculate INFO/AF

for all samples (Total)

for each population (subset)

Merge every AF column into one table:

sql
Copy
Edit
CHROM  POS  REF  ALT  AF_Total  AF_ADG  AF_BPG  …  AF_RYK
Missing frequencies are denoted . (e.g., if all genotypes in that population are missing).

Quick Start
bash
Copy
Edit
# 1. Make the script executable
chmod +x get_AF_by_pop.sh

# 2. Run the script
./get_AF_by_pop.sh

# 3. Output
ls -lh allele_freq_table.tsv
Output Columns
Column	Description
CHROM	Chromosome
POS	1‑based position
REF	Reference allele
ALT	Alternate allele
AF_Total	Allele frequency in the entire dataset
AF_<POP>	Allele frequency in each detected population

Tips & Extensions
Monomorphic sites are retained (via bcftools view -a), so AF can be 0 or 1.

To split multiallelic sites beforehand:

bash
Copy
Edit
bcftools norm -m -any -Oz -o split.vcf.gz input.vcf.gz
tabix -p vcf split.vcf.gz
Adapt the script to read a population map file instead of parsing IDs.

Add more tags (e.g., AC, AN, missing rate) by modifying the bcftools query format string.
