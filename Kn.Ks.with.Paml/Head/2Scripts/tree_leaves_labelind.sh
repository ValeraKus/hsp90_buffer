#!/bin/bash

#ls ../../Body/1_Raw/omm_RooTree.v10b_116tax_CDS_final/ | sed 's/.rootree//' > trees_names_list.txt

ls ../../Body/2_Derived/nonclients/tree/ | sed 's/.fasta.nw//' > trees_names_list.txt


for name in $(cat trees_names_list.txt)
do
inp="../../Body/2_Derived/nonclients/tree/${name}.fasta.nw"
out="../../Body/2_Derived/nonclients/labeled_trees/${name}_leaves_labeled.nw"



cat $inp | sed 's/\([A-Z][a-z]*_[a-z]*_*[a-z]*:[0-9]*\.*[0-9]*\)/ \1 #/g' | perl -pe 'BEGIN{$A=1;} s/#/"#".$A++/ge' > $out

done

rm trees_names_list.txt
