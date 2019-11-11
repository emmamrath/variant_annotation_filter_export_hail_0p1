#!/bin/bash
set -euo pipefail

indir=$1 # "/my/cohort/shard_data_extract"
outdir=$2 # "/my/cohort/data_extract_split"
tmpdir_base=$3 # "./tmp"

reference="../data/reference_genome_hs37d5x/hs37d5x.fa"

mkdir -p "${outdir}"
tmpdir="${tmpdir_base}"/00_03_split
mkdir -p "${tmpdir}"

for shard in {0001..0145}; do

  infile=`ls -1 "${indir}"/*.shard"${shard}".*.vcf.gz`

  infile_basename="$(basename $infile)"
  len=${#infile_basename}
  outfile_prefix=${infile_basename::len-7}
  outfile="${outdir}"/"${outfile_prefix}".split.vcf

  qsub -z -v SHARD="${shard}",INFILE="${infile}",OUTFILE="${outfile}",OUTDIR="${outdir}",REFERENCE="${reference}" -N spl${shard} 04_split_multiallelics.pbs

  outfile_basename=$(basename $OUTFILE)
  tmp_outfile="${tmpdir}"/"${outfile_basename}"
  tmpdir_bcftools="${tmpdir}"/bcftools_"${outfile_basename}" # each run of bcftools needs its own tmpdir because bcftools creates it before running and deletes it after running
  mkdir -p "${tmpdir_bcftools}"

  java -Xmx1800M -Djava.io.tmpdir="${tmpdir}" -jar bin/GenomeAnalysisTK.jar \
    -T LeftAlignAndTrimVariants \
    -R "${reference}" \
    --variant "${infile}" \
    -o "${tmp_outfile}" \
    --splitMultiallelics

  bcftools sort -o "${outfile}" "${tmp_outfile}" -T "${tmpdir_bcftools}"

done

