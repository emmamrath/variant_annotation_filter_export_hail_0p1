#!/bin/bash

echo '##fileformat=VCFv4.2' > temp_hdr.txt
echo '##INFO=<ID=CONDEL_STRAND,Number=1,Type=String,Description="CONDEL STRAND">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_TRANSCRIPT,Number=1,Type=String,Description="CONDEL TRANSCRIPT">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_PROTEIN,Number=1,Type=String,Description="CONDEL PROTEIN">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_AA_POS,Number=1,Type=Integer,Description="CONDEL AA_POS">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_AA_REF,Number=1,Type=String,Description="CONDEL AA_REF">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_AA_ALT,Number=1,Type=String,Description="CONDEL AA_ALT">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_SIFT,Number=1,Type=Float,Description="CONDEL SIFT">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_PPH2,Number=1,Type=Float,Description="CONDEL PPH2">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_MA,Number=1,Type=Float,Description="CONDEL MA">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_FATHMM,Number=1,Type=Float,Description="CONDEL FATHMM">' >> temp_hdr.txt
echo '##INFO=<ID=CONDEL_SCORE,Number=1,Type=Float,Description="CONDEL SCORE">' >> temp_hdr.txt
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> temp_hdr.txt

cat fannsdb.tsv | awk '{FS="\t";OFS=""} {print $1, "\t", $2, "\t.\t", $4, "\t", $5, "\t60\tPASS\tCONDEL_STRAND=", $3, ";CONDEL_TRANSCRIPT=", $6, ";CONDEL_PROTEIN=", $7, ";CONDEL_AA_POS=", $8, ";CONDEL_REF=", $9, ";CONDEL_ALT=", $10, ";CONDEL_SIFT=", $11, ";CONDEL_PPH2=", $12, ";CONDEL_MA=", $13, ";CONDEL_FATHMM=", $14, ";CONDEL_SCORE=", $15}' | sort | uniq > temp_body.txt

cat temp_hdr.txt temp_body.txt > condel_scores.vcf
