#!/bin/bash

# wget http://mendel.stanford.edu/sidowlab/downloads/gerp/hg19.GERP_scores.tar.gz
# tar xzvf hg19.GERP_scores.tar.gz

for infile in chr*.maf.rates; do

	IFS='.' read -r -a array <<< "$infile"
	bit1="${array[0]}"
	chrom=${bit1:3}
	outfile=gerp.chr"${chrom}".maf.rates.bed
	echo 'processing' $infile 'to produce' $outfile
	awk -v chrom="$chrom" 'BEGIN {FS="\t";OFS="\t"} {print chrom, NR, $0}' "${infile}" > "${outfile}"
done

sed -e 's/^M\t/MT\t/' gerp.chrM.maf.rates.bed > gerp.chrMT.maf.rates.bed
rm gerp.chrM.maf.rates.bed

