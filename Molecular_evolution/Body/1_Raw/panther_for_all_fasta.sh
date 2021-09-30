#!/bin/bash

ls omm_NT_fasta.v10b_116tax_CDS_final/ | sed 's/.fasta//' > gene_names_list.txt

for name in $(cat gene_names_list.txt)
do
inp="~/Raw/omm_NT_fasta.v10b_116tax_CDS_final/${name}.fasta"
out="~/Derived/Panther_results/${name}_panther_output.txt"
echo "$inp"
echo "$out"
done

