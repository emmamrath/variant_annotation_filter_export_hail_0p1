#!/bin/bash
set -euo pipefail

indir=$1 # "/my/cohort/shard_data"
outdir=$2 # "/my/cohort/shard_data_extract"
tmpdir_base=$3 # "./tmp"

in_exon_list="../data/UCSC_tables_GRCh37_RefSeq_genes_canonical_exons_bedtools_merged.txt"

mkdir -p "${outdir}"
tmpdir="${tmpdir_base}"/00_02_grep_tabix
mkdir -p "${tmpdir}"

for shard in {0001..0145}; do

  infile=`ls -1 "${indir}"/*.shard"${shard}".*.vcf.gz`
  infile_basename="$(basename $infile)"

  IFS='.' read -r -a array <<< "$infile_basename"
  for element in "${array[@]}"; do
    if  [[ $element == shard* ]]; then
        shard=`echo $element | cut -c6-`
    fi
    if  [[ $element == chrom* ]]; then
        chrom=`echo $element | cut -c6-`
    fi
    if  [[ $element == pos* ]]; then
        start_end=`echo $element | cut -c4-`
        IFS='-' read -r -a array2 <<< "$start_end"
        start_pos="${array2[0]}"
        end_pos="${array2[1]}"
    fi
  done

  infile_base=${infile_basename%???????}
  outfile_basename="${infile_base}.extract.vcf"
  outfile="${outdir}/${outfile_basename}"

  in_vcf_header="${tmpdir}"/temp_header_shard"${shard}".vcf
  tabix -H "${infile}" > "${in_vcf_header}"

  chrom_for_grep="^${CHROM}\t"

  tmpfile="${tmpdir}"/"${infile_base}.exon_intervals.txt"
  rm -rf $tmpfile
  grep -P "${chrom_for_grep}" "$in_exon_list" | \
    awk -v chrom="$chrom" -v startpos="$startpos" -v endpos="$endpos" 'BEGIN {FS="\t";OFS=""} (((startpos<=$2)&&(endpos>=$3))||((endpos>=$2)&&(endpos<=$3))||((startpos>=$2)&&(startpos<=$3))||((startpos>=$2)&&(endpos<=$3))) {print $1, ":", $2, "-", $3}' >> "${tmpfile}"

  :>"${outfile}"
  while IFS= read -r inline; do
    echo $inline $infile_basename
    tabix $infile $inline >> $outfile
  done < "$tmpfile"

  sort -k1,1 -k2,2n -k4,4 -k5,5 "${outfile}" | uniq | cat "${in_vcf_header}" - | bgzip > "${outfile}".gz
  tabix -p vcf "${outfile}".gz

done

