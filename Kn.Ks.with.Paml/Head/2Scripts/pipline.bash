!#/bin/bash

#1 Downloading data
#wget -P ../../Body/1_Raw http://www.orthomam.univ-montp2.fr/orthomam_v10b/cds/archives/omm_NT_fasta.v10c_116tax_CDS_final.tar.gz
wget -P ../../Body/1_Raw http://www.orthomam.univ-montp2.fr/orthomam_v10b/cds/archives/omm_RooTree.v10b_116tax_CDS_final.tar.gz
wget -P ../../Body/1_Raw http://www.orthomam.univ-montp2.fr/orthomam_v10b/cds/archives/omm_AA_fasta.v10c_116tax_CDS_final.tar.gz

#2 Decompress data
tar -xzvf ../../Body/1_Raw/*.tar.gz

#3 Remoove archives
rm ../../Body/1_Raw/*.tar.gz
