#!/bin/bash
set -e
cd /path/to/output/directory/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Bacteroides_vulgatus/assembly_summary.txt
genomeNum=$(grep -c "." assembly_summary.txt)
if [ $genomeNum -gt 700 ]; then
grep "Complete\|Chromosome" assembly_summary.txt | cut -f20 > var.txt
else
cut -f20 assembly_summary.txt | sed '1,2d' > var.txt
fi

for f in `cat var.txt`
do
name=$(grep -w "$f" assembly_summary.txt | cut -f9 | cut -f2 -d'=' | sed 's/ /_/g' | sed 's/\//_/g' | sed 's/\:/_/g' | sed 's/)/_/g' | sed 's/(/_/g')
xx=$(grep -w "$f" assembly_summary.txt | cut -f20 | cut -f10 -d'/')
wget $f/$xx\_genomic.fna.gz
gunzip *.gz
rename $xx\_genomic.fna $name.fna $xx\_genomic.fna
done
rm var.txt
rm assembly_summary.txt
