#!/usr/bin/env bash
set -euo pipefail

###############################################################################
VCF="or1_198_bi.vcf.gz"      # bgzipped & indexed VCF
THREADS=4
POP_TAG_START=7              # 0‑based index of pop code in sample ID
POP_TAG_LEN=3
###############################################################################

command -v bcftools >/dev/null 2>&1 || { echo "Need bcftools"; exit 1; }

TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

echo "➤ Splitting sample IDs by population tag …"
bcftools query -l "${VCF}" | while read -r id; do
    pop=${id:${POP_TAG_START}:${POP_TAG_LEN}}
    echo "${id}" >> "${TMPDIR}/${pop}.samples"
done

mapfile -t POP_NAMES < <(ls "${TMPDIR}"/*.samples | xargs -n1 basename | \
                         sed 's/.samples//' | sort)
echo "   Populations detected: ${POP_NAMES[*]}"

get_af () {
    local label=$1   # pop code or Total
    local samplist=$2
    local outfile=$3

    if [[ -z "${samplist}" ]]; then
        # ---- AF for all samples ------------------------------------------------
        bcftools +fill-tags "${VCF}" -Ou -- -t AF |
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > "${outfile}"
    else
        # ---- AF for population subset -----------------------------------------
        # -a  keeps sites even if monomorphic; -S  selects samples
        bcftools view -Ou -a -S "${samplist}" "${VCF}" | \
        bcftools +fill-tags -Ou -- -t AF | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > "${outfile}"
    fi
}

echo "➤ Computing AF (Total) …"
get_af "Total" "" "${TMPDIR}/Total.af"

for pop in "${POP_NAMES[@]}"; do
    echo "   ↳ ${pop}"
    get_af "${pop}" "${TMPDIR}/${pop}.samples" "${TMPDIR}/${pop}.af"
done

echo "➤ Merging tables …"
paste_files=( "${TMPDIR}/Total.af" )
for pop in "${POP_NAMES[@]}"; do paste_files+=( "${TMPDIR}/${pop}.af" ); done

paste "${paste_files[@]}" | \
awk -v pops="Total ${POP_NAMES[*]}" '
BEGIN {
    split(pops, a, " ");
    printf "CHROM\tPOS\tREF\tALT";
    for (i = 1; i <= length(a); i++) printf "\tAF_%s", a[i];
    printf "\n";
}
{
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4;
    col = 5;
    while (col <= NF) {
        if ($(col) == "")
            printf "\t.";
        else
            printf "\t%s", $(col);
        col += 5;
    }
    printf "\n";
}' > allele_freq_table.tsv
