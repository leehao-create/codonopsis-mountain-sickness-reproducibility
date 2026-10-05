#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module_dir="$(cd -- "${script_dir}/.." && pwd)"
output="${NCBI_GENE_MAP:-${module_dir}/data/reference/NCBI_GeneID_to_symbol.txt}"
url="https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"

mkdir -p "$(dirname -- "${output}")"
tmp_file="$(mktemp)"
trap 'rm -f "${tmp_file}"' EXIT

curl --fail --location --retry 3 "${url}" --output "${tmp_file}"
{
  printf 'NCBI gene (formerly Entrezgene) ID\tGene name\n'
  gzip -dc "${tmp_file}" | awk -F '\t' 'NR > 1 {print $2 "\t" $3}'
} > "${output}"

printf 'Wrote current NCBI Homo_sapiens.gene_info mapping to %s\n' "${output}"
printf 'Record the download date and SHA-256 before analysis; current NCBI content may differ from the archived mapping.\n'
