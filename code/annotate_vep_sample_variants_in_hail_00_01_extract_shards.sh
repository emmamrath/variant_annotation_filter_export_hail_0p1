#!/bin/bash
set -euo pipefail

indir=$1 # /my/cohort/before_shards
outdir=$2 # /my/cohort/shard_data
tmpdir_base=$3 # ./tmp

infile="../data/hs37d5x_shards_4cols.txt"
reference="../data/reference_genome_hs37d5x/hs37d5x.fa"

mkdir -p "${outdir}"
tmpdir="${tmpdir_base}"/00_01_selectvariants
mkdir -p "${tmpdir}"

outstring=""
while IFS= read -r inline; do

  IFS=$'\t' read -r -a array <<< "$inline"
  shard="${array[0]}"
  chrom="${array[1]}"
  start="${array[2]}"
  end="${array[3]}"
  interval="${chrom}":"${start}"-"${end}"

  echo 'Processing shard' $shard

  for infile in "${indir}/"*.vcf.gz; do

    echo 'Processing sample' $infile 'shard' $shard

    infile_basename=$(basename $infile)
    IFS=$'.' read -r -a array2 <<< "$infile_basename"
    sampleid="${array2[0]}"
    outfile="${outdir}"/"${sampleid}".WholeGenome.shard"${shard}".chrom"${chrom}".pos"${start}"-"${end}".vcf.gz

    java -Xmx1800M -Djava.io.tmpdir="${tmpdir}" -jar bin/GenomeAnalysisTK.jar \
       -T SelectVariants \
       -R "${reference}" \
       -V "${infile}" \
       -o "${outfile}" \
       -L "${interval}"

  done
done < "$infile"

