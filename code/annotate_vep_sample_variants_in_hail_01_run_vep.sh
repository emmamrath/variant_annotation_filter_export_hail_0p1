#!/bin/bash

shard1=$1 # 0001
shard2=$2 # 0145
indir=$3 # /my/data/directory/my_cohort_data_extract_split
tmpdir=$4 # ./tmp

vep_executable="${5}" # /my/vep/installation/ensembl-vep/vep
vep_cache_data="${6}" # /my/vep/installation/vep_data
vep_fasta="${7}" # /my/vep/installation/vep_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
vep_plugins="${8}" # /my/dir/.vep/Plugins
bcftools_executable="${9}" # /my/bcftools/installation/bcftools/install/bin/bcftools
vep_loftee_plugin_string="${10}" # LoF,loftee_path:/my/loftee/installation/LOFTEE/loftee,human_ancestor_fa:/my/loftee/installation/human_ancestor_data/human_ancestor.fa.gz,conservation_file:/my/loftee/installation/conservation_data/phylocsf_gerp.sql,gerp_file:/my/loftee/installation/gerp_scores_data_recompress_bgz/GERP_scores.final.sorted.txt.bgz
vep_loftool_plugin_string="${11}" # LoFtool,/my/vep/installation/vep_plugins/download_plugin_data/LoFtool_scores.txt
vep_revel_plugin_string="${12}" # REVEL,/my/vep/installation/vep_plugins/revel_data/new_tabbed_revel.tsv.gz

outdir1="${indir}"_sorted
outdir2="${indir}"_sorted_vep"
outdir3="${indir}"_sorted_vep_reformatted"

mkdir -p "${outdir1}"
mkdir -p "${outdir2}"
mkdir -p "${outdir3}"
mkdir -p "${tmpdir}"

for shard in $(eval echo {$shard1..$shard2}); do

	num_infiles=`ls -1 "${indir}"/*.shard"${shard}".*.vcf 2>/dev/null | wc -l`
	echo ''; echo 'num_infiles' $num_infiles

	if [ "$num_infiles" -gt 0 ]; then

		infile=`ls -1 "${indir}"/*.shard"${shard}".*.vcf`

		echo ''; echo 'sort' $infile 'before vep'

		infile_basename=$(basename $infile)
		infile_prefix=${infile_basename::-13}
		outfile1="${outdir1}"/"${infile_prefix}".sorted.vcf

		tmpdir_bcftools="${tmpdir}"/"${shard}" # bcftools deletes this tmpdir after using it, thus each shard needs its own temporary bcftools directory
		mkdir -p "${tmpdir_bcftools}"
		"${bcftools_executable}" sort -o "${outfile1}" "${infile}" -T "${tmpdir_bcftools}"

		echo ''; echo 'run vep on' $outfile1

		outfile2="${outdir2}"/"${infile_prefix}".sorted.vep.vcf
		tmpfile_hdr_vep="${tmpdir}"/"${infile_prefix}"_hdr_vep.txt

		num_variants=`grep -v '^#' "${outfile1}" | wc -l`
		num_loop=$(($num_variants/1000 + 1))

		tmpfile_sorted_hdr="${tmpdir}"/"${infile_prefix}"_sorted_hdr.txt
		grep '^#' "${outfile1}" > "${tmpfile_sorted_hdr}"

		from_i=0
		to_i=1
		while [ $to_i -le $num_loop ]; do
			from_var=$(($from_i * 1000 + 1))
			to_var=$(($to_i * 1000))
			#print 'from_i' $from_i 'to_i' $to_i 'from_var' $from_var 'to_var' $to_var
			tmpfile_subset="${tmpdir}"/"${infile_prefix}"_subset_"${from_var}"_"${to_var}".before_vep.vcf
			tmpfile_subset_vep="${tmpdir}"/"${infile_prefix}"_subset_"${from_var}"_"${to_var}".after_vep.vcf
			if [ $to_i -eq 1 ]; then
				grep -v '^#' "${outfile1}" | head -n 1000 | cat "${tmpfile_sorted_hdr}" - > "${tmpfile_subset}"
			else
				if [ $to_i -eq $num_loop ]; then
					to_var=$num_variants
					tmpfile_subset="${tmpdir}"/"${infile_prefix}"_subset_"${from_var}"_"${to_var}".before_vep.vcf
					tmpfile_subset_vep="${tmpdir}"/"${infile_prefix}"_subset_"${from_var}"_"${to_var}".after_vep.vcf
					tail_amt=$(($to_var - $from_var + 1))
					grep -v '^#' "${outfile1}" | tail -n "${tail_amt}" | cat "${tmpfile_sorted_hdr}" - > "${tmpfile_subset}"
				else
					grep -v '^#' "${outfile1}" | head -n "${to_var}" | tail -n 1000 | cat "${tmpfile_sorted_hdr}" - > "${tmpfile_subset}"
				fi
			fi

			echo ''; echo 'call vep for' $tmpfile_subset 'to produce' $tmpfile_subset_vep
			"${vep_executable}" --vcf --offline --cache \
				-dir_cache "${vep_cache_data}" \
				--fasta "${vep_fasta}" \
				-i "${tmpfile_subset}" -o "${tmpfile_subset_vep}" --everything \
				--plugin "${vep_loftee_plugin_string}" \
				--plugin "${vep_loftool_plugin_string}" \
				--plugin "${vep_revel_plugin_string}" \
				--hgvs --clin_sig_allele --mane --plugin LOVD --plugin CSN --plugin SpliceRegion --plugin TSSDistance --plugin miRNA --force_overwrite --dir_plugins "${vep_plugins}"

			((from_i++))
			((to_i++))
		done

		echo ''; echo 'cat' "${tmpdir}"'/'"${infile_prefix}"'_subset_*_*.after_vep.vcf, sort, produce' $outfile2
		grep '^#' "${tmpdir}"/"${infile_prefix}"_subset_1_*.after_vep.vcf > "${tmpfile_hdr_vep}"
		cat "${tmpdir}"/"${infile_prefix}"_subset_*_*.after_vep.vcf | grep -v '^#' | sort -k1,1 -k2,2n -k4,4 -k5,5 | uniq | cat "${tmpfile_hdr_vep}" - > "${outfile2}"

		echo ''; echo 'run reformat vep file' $outfile2

		outfile3="${outdir3}"/"${infile_prefix}".sorted.vep_reformatted.vcf

		echo 'python convert_vcf_vep_info_into_standard_info_fields.py -i '"${outfile2}"' -o '"${outfile3}"
		python convert_vcf_vep_info_into_standard_info_fields.py -i "${outfile2}" -o "${outfile3}" "${tmpdir}" "${tmpdir}"/hail.log

		echo ''; echo 'finished this shard' $outfile3
	fi
done

