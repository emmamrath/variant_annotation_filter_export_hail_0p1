#!/bin/bash

# cd data/clinvar
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
# gunzip variant_summary.txt.gz
# mv variant_summary.txt clinvar_variant_summary.txt
# cd ../../code
# ./create_vcf_from_clinvar_file.sh

python create_VCF_from_clinvar_variant_summary_TSV_file_using_UCSC_genes.py -i ../data/clinvar/clinvar_variant_summary.txt -r ../data/reference_genome_hs37d5x/hs37d5x.fa -g ../data/reference_RefSeq_genes/UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq.txt -o ../data/clinvar/clinvar_variant_summary.vcf
grep '^#' ../data/clinvar/clinvar_variant_summary.vcf > ../data/clinvar/temp_clinvar_hdr.vcf
grep -v '^#'  ../data/clinvar/clinvar_variant_summary.vcf | sort -k1,1 -k2,2n >  ../data/clinvar/temp_clinvar_body.vcf
cat  ../data/clinvar/temp_clinvar_hdr.vcf  ../data/clinvar/temp_clinvar_body.vcf >  ../data/clinvar/clinvar_variant_summary_sorted.vcf
rm  ../data/clinvar/temp_clinvar_hdr.vcf  ../data/clinvar/temp_clinvar_body.vcf
