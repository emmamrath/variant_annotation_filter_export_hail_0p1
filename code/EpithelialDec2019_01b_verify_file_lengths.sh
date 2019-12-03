#!/bin/bash

#indir=../MGRB_ISKS_RegulatoryRegionsOct2019_extract_gnomadAFnfe_lt_0p0001_and_gnomadAF_lt_0p005_split
#outdir=../MGRB_ISKS_RegulatoryRegionsOct2019_extract_gnomadAFnfe_lt_0p0001_and_gnomadAF_lt_0p005_split_sorted_vep_reformatted
#indir=../MGRB3_ISKS2_RISCsarcoma_LionsKidsSarcoma.shards_exons_split
#outdir=../MGRB3_ISKS2_RISCsarcoma_LionsKidsSarcoma.shards_exons_split_sorted_vep_reformatted
indir=../EpithelialDec2019_01_run_vep.sh.shards_exons_split
outdir=../EpithelialDec2019_01_run_vep.sh.shards_exons_split_sorted_vep_reformatted

echo 'For each shard, this script will display the number of variants in shards_exons_split and shards_exons_split_sorted_vep_reformatted. They should be the same.'

for shard in {0001..0145}; do
	echo $shard
	infile_length=$(grep -v "^#" "${indir}"/*.shard"${shard}".*.split.vcf | wc -l | cut -d" " -f1)
	outfile_length=$(grep -v "^#" "${outdir}"/*.shard"${shard}".*.sorted.vep_reformatted.vcf | wc -l | cut -d" " -f1)
	echo $shard $infile_length $outfile_length
done

echo 'For each shard, this script has displayed the number of variants in shards_exons_split and shards_exons_split_sorted_vep_reformatted. They should be the same.'


