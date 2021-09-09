#!/bin/bash

cat Body/1_Raw/human.hsp.genes.txt | cut -f 10 | tail -97 > HSP/Body/2_Derived/human.hsp.ensemblid.list.txt 
