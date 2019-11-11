#!/bin/bash

inlist=$1 # "list_shards_for_hail.txt"
indir=$2 # "/my/cohort/data_extract_split_sorted_vep_reformatted" # this indir is specified in inlist file too
outdir=$3 # "/my/cohort/data_extract_split_sorted_vep_reformatted_hail" # this outdir is specified in inlist file too
caddlist=$4 # "list_cadd_shards.txt"
gerplist=$5 # "list_gerp_shards.txt"
gnomadlist=$6 # "list_gnomad_shards.txt"

mkdir -p "${outdir}"

# cat list_shards_for_hail.txt | head
# shard0001	/my/cohort/data_extract_split_sorted_vep_reformatted/my_cohort.shard0001.chrom1.pos1-13077999.split.sorted.vep_reformatted.vcf	/my/cohort/data_extract_split_sorted_vep_reformatted_hail/my_cohort.shard0001.chrom1.pos1-13077999.split.sorted.vep_reformatted_hail
# shard0002	/my/cohort/data_extract_split_sorted_vep_reformatted/my_cohort.shard0002.chrom1.pos13078000-17150659.split.sorted.vep_reformatted.vcf	/my/cohort/data_extract_split_sorted_vep_reformatted_hail/my_cohort.shard0002.chrom1.pos13078000-17150659.split.sorted.vep_reformatted_hail

# cat list_cadd_shards.txt | head
# shard0001	/my/hail_databases/CADD_1_4_GRCh37/cadd_vcf_bgz_files/cadd.shard0001.chrom1.pos1-13077999.exon_intervals.vcf.bgz
# shard0002	/my/hail_databases/CADD_1_4_GRCh37/cadd_vcf_bgz_files/cadd.shard0002.chrom1.pos13078000-17150659.exon_intervals.vcf.bgz

# cat list_gerp_shards.txt | head
# shard0001	/my/hail_databases/gerp_db/gerp.shard0001.chrom1.pos1-13077999.vcf.bgz
# shard0002	/my/hail_databases/gerp_db/gerp.shard0002.chrom1.pos13078000-17150659.vcf.bgz

python annotate_vep_sample_variants_in_hail_in_subsets_with_vep_fields_already_filled_for_all_shards.py -invcf_list "${inlist}" -cadd_list "${caddlist}" -gerp_list "${gerplist}"

