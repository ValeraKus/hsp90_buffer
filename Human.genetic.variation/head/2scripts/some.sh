#!/bin/zsh

##hsp90 coordinates in human genome hg19 
ensID='ENSG00000096384'
start=$(awk -v id=$ensID 'BEGIN { FS="\t" } {if ($0 !~ "#" && $3 == "gene" && $9 ~ id) print $4}' ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/1raw/gencode.v19.annotation.gtf) 
#> ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/2derived/HSP90AB1.hg19.coordinates.txt)

chr=$(awk -v id=$ensID 'BEGIN { FS="\t" } {if ($0 !~ "#" && $3 == "gene" && $9 ~ id) print $1}' ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/1raw/gencode.v19.annotation.gtf | sed 's/chr//')

#start=$(cut -f1 ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/2derived/HSP90AB1.hg19.coordinates.txt)
#end=$(cut -f2 ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/2derived/HSP90AB1.hg19.coordinates.txt)
echo $start
echo $chr
#echo $end

awk -v st=$start -v c=$chr '{if ($0 !~ "#" && $2==c && $3>=st-1000000 && $3 <=st+1000000) print $0}' ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/1raw/PGS000713.txt > ~/Desktop/Work/hsp90_buffer/Human.genetic.variation/body/2derived/$ensID.PGS000713.snps.txt


