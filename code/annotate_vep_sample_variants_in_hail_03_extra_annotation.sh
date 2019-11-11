#!/bin/bash

shard1=$1 # "0001"
shard2=$2 # "0145"
indir=$3 # "/my/cohort/data_extract_split_sorted_vep_reformatted_hail"
outdir=$4 # "/my/cohort/data_extract_split_sorted_vep_reformatted_hail_extra"

mkdir -d "${outdir}"

gene_file="../data/UCSC_tables_GRCh37_RefSeq_genes_gene_length_for_canonical_and_longest_transcripts.txt"

for shard in $(eval echo {$shard1..$shard2}); do

	num_infiles=`ls -1 "${indir}"/*.shard"${shard}".*.tsv 2>/dev/null | wc -l`
	echo ''; echo 'num_infiles' $num_infiles

	if [ "$num_infiles" -gt 0 ]; then

		infile=`ls -1 "${indir}"/*.shard"${shard}".*.tsv`

		infile_basename=$(basename $infile)
		infile_prefix=${infile_basename::-4}
                outfile="${outdir}"/"${infile_prefix}".extra.tsv

		echo 'python add_gene_annotation_to_tab_delimited_file.py -i '"${infile}"' -o '"${outfile}" -num_aa "${gene_file}"
		python add_gene_annotation_to_tab_delimited_file.py -i "${infile}" -o "${outfile}" -num_aa "${gene_file}"

		echo 'finished this shard' $outfile; echo ''
	fi
done
