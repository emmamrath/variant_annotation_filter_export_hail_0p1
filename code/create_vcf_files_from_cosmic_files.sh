#!/bin/bash

# cd data/cosmic
# sftp -o User=<my_email@my_email.org> sftp-cancer.sanger.ac.uk
# get /files/grch37/cosmic/v84/CosmicMutantExport.tsv.gz
# get /files/grch37/cosmic/v84/CosmicResistanceMutations.tsv.gz
# get /files/grch37/cosmic/v84/VCF/CosmicCodingMuts.vcf.gz
# get /files/grch37/cosmic/v84/VCF/CosmicNonCodingVariants.vcf.gz
# get /files/grch37/cosmic/v84/CosmicHGNC.tsv.gz
# cd ../../code
# ./create_vcf_files_from_cosmic_files.sh

grep '^#' ../data/cosmic/CosmicCodingMuts.vcf > ../data/cosmic/temp_CosmicCodingMuts_hdr.vcf
grep 'GENE=' ../data/cosmic/CosmicCodingMuts.vcf > ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_body.vcf

cat ../data/cosmic/temp_CosmicCodingMuts_hdr.vcf | sed 's/##INFO=<ID=GENE/##INFO=<ID=COSMICCODINGMUTS_GENE/g' | sed 's/##INFO=<ID=STRAND/##INFO=<ID=COSMICCODINGMUTS_STRAND/g' | sed 's/##INFO=<ID=CDS/##INFO=<ID=COSMICCODINGMUTS_CDS/g' | sed 's/##INFO=<ID=AA/##INFO=<ID=COSMICCODINGMUTS_AA/g' | sed 's/##INFO=<ID=CNT/##INFO=<ID=COSMICCODINGMUTS_CNT/g' > temp_CosmicCodingMuts_hdr2.vcf
cat ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_body.vcf | sed 's/GENE=/COSMICCODINGMUTS_GENE=/g' | sed 's/STRAND=/COSMICCODINGMUTS_STRAND=/g' | sed 's/CDS=/COSMICCODINGMUTS_CDS=/g' | sed 's/AA=/COSMICCODINGMUTS_AA=/g' | sed 's/CNT=/COSMICCODINGMUTS_CNT=/g' > ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_body2.vcf
cat ../data/cosmic/temp_CosmicCodingMuts_hdr2.vcf ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_body2.vcf > ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_tooManyInfoFields.vcf

python filter_info_fields_in_vcf.py -i ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_tooManyInfoFields.vcf -f 'COSMICCODINGMUTS_GENE,COSMICCODINGMUTS_STRAND,COSMICCODINGMUTS_CDS,COSMICCODINGM
UTS_AA,COSMICCODINGMUTS_CNT' -o ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail.vcf
vt normalize ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail.vcf -n -r ../data/reference_genome_hs37d5x/hs37d5x.fa -o ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize_withDuplicates.vcf

grep '^#' ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize_withDuplicates.vcf > ../data/cosmic/temp_hdr.vcf
grep -v '^#' ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize_withDuplicates.vcf > ../data/cosmic/temp_body.vcf
cat ../data/cosmic/temp_body.vcf | sort -k1,1 -k2,2n | uniq > ../data/cosmic/temp_body_uniq.vcf
cat ../data/cosmic/temp_hdr.vcf temp_body_uniq.vcf > ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf
bgzip -c ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf > ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf.gz
tabix -f -p vcf ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf.gz
ln -s ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf.gz ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf.bgz
ln -s ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf.gz.tbi ../data/cosmic/CosmicCodingMuts_genesOfInterest_for_Hail_vtnormalize.vcf.bgz.tbi

# We have more than 1 entry per SNP.
# They are all GRCh37. They are different transcripts. 2 have long gene_name. Get only chopped gene_name.
# They are from different samples. Some are same SNP, some are not.
# I will convert the CosmicMutantExport_genesOfInterest.tsv Mutation_CDS to VCF.REF and VCF.ALT and when annotating take only the first one that matches.
# Half the reverse-strand genes have the reverse strand in the Mutation_genome_position, and thus my program will reverse-complement the reverse-strand nucleotides.
# The other half of the reverse-strand genes have the forward nucleotides.

head -n 1 ../data/cosmic/CosmicMutantExport.tsv > ../data/cosmic/temp_CosmicMutantExport_hdr.tsv
awk 'FNR != 1 {print}' ../data/cosmic/CosmicMutantExport.tsv > ../data/cosmic/temp_CosmicMutantExport_body.tsv
cat ../data/cosmic/temp_CosmicMutantExport_hdr.tsv ../data/cosmic/temp_CosmicMutantExport_body.tsv > ../data/cosmic/temp_CosmicMutantExport_genesOfInterest_for_Hail.tsv

python create_vcf_from_CosmicMutantExport_tsv_file.py -i ../data/cosmic/temp_CosmicMutantExport_genesOfInterest_for_Hail.tsv -r ../data/reference_genome_hs37d5x/hs37d5x.fa -g ../data/UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq.txt -o ../data/cosmic/temp_CosmicMutantExport_genesOfInterest_for_Hail.vcf

grep '^#' ../data/cosmic/temp_CosmicMutantExport_genesOfInterest_for_Hail.vcf > ../data/cosmic/temp_hdr.vcf
grep -v '^#' ../data/cosmic/temp_CosmicMutantExport_genesOfInterest_for_Hail.vcf | sort -k1,1 -k2,2n > ../data/cosmic/temp_body.vcf
cat ../data/cosmic/temp_hdr.vcf ../data/cosmic/temp_body.vcf > ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_sorted.vcf
bgzip -c ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_sorted.vcf > ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_sorted.vcf.gz
tabix -f -p vcf ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_sorted.vcf.gz

vt normalize ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_sorted.vcf.gz -n -r ../data/reference_genome_hs37d5x/hs37d5x.fa -o ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize_withDuplicates.vcf

grep '^#' ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize_withDuplicates.vcf > ../data/cosmic/temp_hdr.vcf
grep -v '^#' ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize_withDuplicates.vcf > ../data/cosmic/temp_body.vcf
cat ../data/cosmic/temp_body.vcf | sort -k1,1 -k2,2n | uniq > ../data/cosmic/temp_body_uniq.vcf
cat ../data/cosmic/temp_hdr.vcf ../data/cosmic/temp_body_uniq.vcf > ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf

mv ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize_someRefEqualsAlt.vcf

grep '^#' ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize_someRefEqualsAlt.vcf > ../data/cosmic/temp_hdr.vcf
grep -v '^#' ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize_someRefEqualsAlt.vcf | grep -vP '\tR\t' | grep -vP '\tM\t' > ../data/cosmic/temp_body_someRefEqualsAlt.vcf
awk '{FS="\t";OFS="\t";} ($4 != $5) {print $0}' ../data/cosmic/temp_body_someRefEqualsAlt.vcf > ../data/cosmic/temp_body.vcf
cat ../data/cosmic/temp_hdr.vcf ../data/cosmic/temp_body.vcf > ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf

bgzip -c ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf > ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf.gz
tabix -f -p vcf ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf.gz

ln -s ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf.gz ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf.bgz
ln -s ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf.gz.tbi ../data/cosmic/CosmicMutantExport_genesOfInterest_for_Hail_vtnormalize.vcf.bgz.tbi

rm ../data/cosmic/temp_CosmicCodingMuts_hdr.vcf ../data/cosmic/temp_CosmicCodingMuts_genesOfInterest_body.vcf ../data/cosmic/temp_hdr.vcf../data/cosmic/temp_body.vcf ../data/cosmic/temp_body_someRefEqualsAlt.vcf

