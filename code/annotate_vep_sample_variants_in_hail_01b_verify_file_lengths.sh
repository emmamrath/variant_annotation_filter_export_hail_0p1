#!/bin/bash

indir=$3 # /my/data/directory/my_cohort_data_extract_split
outdir=$3 # /my/data/directory/my_cohort_data_extract_split_sorted_vep_reformatted

echo 'For each shard, this script will display the number of variants in shards_exons_split and shards_exons_split_sorted_vep_reformatted. They should be the same.'

for shard in {0001..0145}; do
	echo $shard
	infile_length=$(grep -v "^#" "${indir}"/*.shard"${shard}".*.split.vcf | wc -l | cut -d" " -f1)
	outfile_length=$(grep -v "^#" "${outdir}"/*.shard"${shard}".*.sorted.vep_reformatted.vcf | wc -l | cut -d" " -f1)
	echo $shard $infile_length $outfile_length
done

echo 'For each shard, this script has displayed the number of variants in shards_exons_split and shards_exons_split_sorted_vep_reformatted. They should be the same.'


