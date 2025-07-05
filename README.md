# get-AF-by-Pop
Calculate allele frequency of all variants by population in a VCF file
# Allele Frequency Table by Population

This script computes **allele frequencies (AF)** for all variants in a bgzipped VCF file, both:

- ✅ Across the **entire dataset** ("AF_Total")
- ✅ Separately for **each population**, inferred from sample names (e.g., "AF_BPG", "AF_ADG", etc.)

The output is a tab-delimited table:  
`allele_freq_table.tsv`

---

## 🔧 Requirements

- `bash` (v3+)
- [`bcftools`](https://samtools.github.io/bcftools/) (v1.9 or higher)
  - Must be in your `$PATH`
- Your VCF must:
  - Be **bgzipped** and **indexed** (`.vcf.gz` + `.vcf.gz.tbi`)
  - Contain **GT fields** (genotypes for individuals)
  - Use sample names with **embedded population codes**

---

## 📁 Input Files

### 1. `VCF`

- Set in the script via:

  ```bash
  VCF="your_file.vcf.gz"
2. Sample naming convention
The script extracts population names from characters 8–10 of each sample ID.
For example:

nginx
Copy
Edit
HG_ADG_001  →  ADG
HG_BPG_022  →  BPG
You can adjust this logic using these parameters at the top of the script:

bash
Copy
Edit
POP_TAG_START=7     # 0-based index
POP_TAG_LEN=3
🧪 What the Script Does
Reads sample IDs from the VCF

Groups them into populations based on the embedded tag (e.g., ADG, RUS)

Computes allele frequency for:

All samples → AF_Total

Each population subset → AF_POPNAME

Merges all frequencies into a single table:

sql
Copy
Edit
CHROM  POS  REF  ALT  AF_Total  AF_ADG  AF_BPG  …  AF_RYK
🚀 Running the Script
Make executable:

bash
Copy
Edit
chmod +x get_AF_by_pop.sh
Run it:

bash
Copy
Edit
./get_AF_by_pop.sh
Result:

Copy
Edit
allele_freq_table.tsv
📌 Output Columns
Column	Description
CHROM	Chromosome
POS	1-based position
REF	Reference allele
ALT	Alternate allele
AF_Total	Allele frequency in entire dataset
AF_<POP>	Allele frequency in each population

🧠 Notes
If some variants are missing for a population (e.g., due to missing genotypes), their AF will be ".".

The script keeps all sites, including monomorphic ones (via bcftools view -a).

Make sure your sample names are unique and consistent in the VCF.

🧼 Cleaning Up
Temporary .samples and .af files are auto-deleted after the script finishes.
