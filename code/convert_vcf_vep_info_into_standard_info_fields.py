#!/usr/bin/python
# python convert_vcf_vep_info_into_standard_info_fields.py -i infile -o output_vcf
# python convert_vcf_vep_info_into_standard_info_fields.py -i 02_ISKS_MGRB_HC.shard0140.chromY.pos9266324-20168886.whole_exome.vep.vcf -o 02_ISKS_MGRB_HC.shard0140.chromY.pos9266324-20168886.whole_exome.vep_reformatted.vcf

# Input vcf contains vep info from command:
#
# ./vep --vcf --offline --cache -dir_cache ./vep_data -i 02_ISKS_MGRB_HC.shard0001.chrom1.pos1-13077999.whole_exome.vcf.gz -o whole_exome_extract_vep/02_ISKS_MGRB_HC.shard0001.chrom1.pos1-13077999.whole_exome.vep.vcf --everything --force_overwrite
#
# vcf header contains the follow header line, as one line:
# ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format:
# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|
# SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|
# gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
# 
# vcf header contains data records in the following format, as one line:
# X	152481538	.	T	C	521.34	.	
# AC=1;AF=0.0001358;AN=7330;BaseQRankSum=0.287;ClippingRankSum=0;DP=106041;ExcessHet=3.0103;FS=0;InbreedingCoeff=-0.0002;MLEAC=1;MLEAF=0.0001358;MQ=60;MQRankSum=0;QD=17.98;ReadPosRankSum=0.147;SOR=0.804;
# CSQ=C|3_prime_UTR_variant|MODIFIER|MAGEA1|ENSG00000198681|Transcript|ENST00000356661|protein_coding|3/3||||1692|||||||-1||SNV|HGNC|6796|YES|||CCDS14720.1|
# ENSP00000349085|P43355|A8IF97|# UPI0000035FCB|||||||||||||||||||||||||||||||	GT:AD:DP:GQ:PL

# A newer version of the vep call:
#
# /my/vep/installation/ensembl-vep/vep --vcf --offline --cache -dir_cache /my/vep/installation/vep_data --fasta /my/vep/installation/vep_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz -i /my/data/my_cohort_extract/my_cohort.shard0004.chrom1.pos29953084-103888907.vcf.gz -o /my/data/my_cohort_extract_vep/my_cohort.shard0004.chrom1.pos29953084-103888907.vep.vcf --everything --plugin LoFtool,/my/vep/installation/vep_plugins/download_plugin_data/LoFtool_scores.txt --plugin REVEL,/my/vep/installation/vep_plugins/revel_data/new_tabbed_revel.tsv.gz --hgvs --clin_sig_allele --mane --plugin LOVD --plugin CSN --plugin SpliceRegion --plugin TSSDistance --plugin miRNA --force_overwrite
#
# ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoFtool|REVEL|LOVD|CSN|SpliceRegion|TSSDistance|miRNA">

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import random
import commands
import argparse
import re
# import subprocess
# from multiprocessing import Pool

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def merge_dictionaries( *dict_args ):
	# Given any number of dicts, shallow copy and merge into a new dict,
	# precedence goes to key value pairs in latter dicts.
	result = {}
	for dictionary in dict_args:
		result.update(dictionary)
	return result

######################################################
def write_new_vep_headers( inline, outfile ):

	bits1 = inline.split('"')
	bit1 = bits1[1]
	bits2 = bit1.split('Format: ')
	hdr_string = bits2[1]
	orig_vep_field_hdrs = hdr_string.split('|')
	vep_field_hdrs = []
	for this_hdr in orig_vep_field_hdrs:
		new_hdr = "VEP_" + this_hdr
		vep_field_hdrs.append( new_hdr )

	extra_vep_field_hdrs = []
	extra_vep_field_hdrs.append( 'VEP_2_Consequence' )
	extra_vep_field_hdrs.append( 'VEP_2_IMPACT' )
	extra_vep_field_hdrs.append( 'VEP_2_SYMBOL' )
	extra_vep_field_hdrs.append( 'VEP_2_Gene' )
	extra_vep_field_hdrs.append( 'VEP_2_Feature' )
	extra_vep_field_hdrs.append( 'VEP_2_CANONICAL' )

	all_vep_field_hdrs = vep_field_hdrs + extra_vep_field_hdrs
	for this_hdr in all_vep_field_hdrs:
		outline = '##INFO=<ID=' + this_hdr + ',Number=1,Type=String,Description="VEP field ' + this_hdr + '">' + "\n"
		outfile.write( outline )

	return vep_field_hdrs, extra_vep_field_hdrs

######################################################
def reformat_vep( old_vep_info, vep_field_hdrs, extra_vep_field_hdrs ):

	# For CDKN2A, fill information for ARF p14 (p53 negative regulator) and 1 for p16 INC4A (retinoblastoma negative regulator).
	# VEP SYMBOLS for the CDKN2A region: C9orf53, CDKN2A, CDKN2B-AS1, RP11-145E5.5, RP11-149I2.4
	# The transcript for CDKN2A_ARF is ENST00000579755
	# The transcript for CKN2A_INK4A is ENST00000498124

	chosen_vep_fields = {}
	chosen_extra_vep_fields = {}
	for this_key in vep_field_hdrs:
		chosen_vep_fields[this_key] = ''
	for this_key in extra_vep_field_hdrs:
		chosen_extra_vep_fields[this_key] = ''

	vep_alts = old_vep_info.split(',')
	vep_alts.reverse() # Process multiple vep transcripts in reverse order so that the first transcript will be chosen, as it will be processed last.
	for this_vep_string in vep_alts:
		vep_fields = {}
		vep_values = this_vep_string.split('|')
		for i in range( 0, len(vep_field_hdrs) ):
			this_key = str(vep_field_hdrs[i])
			this_value = str(vep_values[i])
			vep_fields[this_key] = this_value

		# Make sure CDKN2A is the chosen vep string in the case there are multiple transcripts for SNPs hitting that gene, with ARF in the main fields and INK4A in the secondary fields.
		# As with all genes, the canonical is the preferred transcript, if the canonical is present. Otherwise, another transcript will be used.
		if (vep_fields['VEP_SYMBOL'] in ['C9orf53', 'CDKN2A', 'CDKN2B-AS1', 'RP11-145E5.5', 'RP11-149I2.4']):
			if (vep_fields['VEP_SYMBOL'] == 'CDKN2A'):
				if (vep_fields['VEP_Feature'] == 'ENST00000579755'):
					if (vep_fields['VEP_CANONICAL'] == 'YES'):
						chosen_vep_fields = vep_fields
					elif (chosen_vep_fields['VEP_CANONICAL'] == 'YES'):
						do_nothing = 1
					else:
						chosen_vep_fields = vep_fields
				elif (vep_fields['VEP_Feature'] == 'ENST00000498124'):
					if (vep_fields['VEP_CANONICAL'] == 'YES'):
						chosen_extra_vep_fields['VEP_2_Consequence'] = vep_fields['VEP_Consequence']
						chosen_extra_vep_fields['VEP_2_IMPACT'] = vep_fields['VEP_IMPACT']
						chosen_extra_vep_fields['VEP_2_SYMBOL'] = vep_fields['VEP_SYMBOL']
						chosen_extra_vep_fields['VEP_2_Gene'] = vep_fields['VEP_Gene']
						chosen_extra_vep_fields['VEP_2_Feature'] = vep_fields['VEP_Feature']
						chosen_extra_vep_fields['VEP_2_CANONICAL'] = vep_fields['VEP_CANONICAL']
					elif (chosen_extra_vep_fields['VEP_2_CANONICAL'] == 'YES'):
						do_nothing = 1
					else:
						chosen_extra_vep_fields['VEP_2_Consequence'] = vep_fields['VEP_Consequence']
						chosen_extra_vep_fields['VEP_2_IMPACT'] = vep_fields['VEP_IMPACT']
						chosen_extra_vep_fields['VEP_2_SYMBOL'] = vep_fields['VEP_SYMBOL']
						chosen_extra_vep_fields['VEP_2_Gene'] = vep_fields['VEP_Gene']
						chosen_extra_vep_fields['VEP_2_Feature'] = vep_fields['VEP_Feature']
						chosen_extra_vep_fields['VEP_2_CANONICAL'] = vep_fields['VEP_CANONICAL']

			if ((vep_fields['VEP_SYMBOL'] == 'CDKN2A') and (chosen_vep_fields['VEP_SYMBOL'] != 'CDKN2A')):
				save_second_transcript = {}
				save_second_transcript['VEP_2_Consequence'] = vep_fields['VEP_Consequence']
				save_second_transcript['VEP_2_IMPACT'] = vep_fields['VEP_IMPACT']
				save_second_transcript['VEP_2_SYMBOL'] = vep_fields['VEP_SYMBOL']
				save_second_transcript['VEP_2_Gene'] = vep_fields['VEP_Gene']
				save_second_transcript['VEP_2_Feature'] = vep_fields['VEP_Feature']
				save_second_transcript['VEP_2_CANONICAL'] = vep_fields['VEP_CANONICAL']

		# The canonical is the preferred transcript, if the canonical is present. Otherwise, another transcript will be used.
		else:
			if (vep_fields['VEP_CANONICAL'] == 'YES'):
				chosen_vep_fields = vep_fields
			elif (chosen_vep_fields['VEP_CANONICAL'] == 'YES'):
				do_nothing = 1
			else:
				chosen_vep_fields = vep_fields

	new_vep_info = ''
	final_vep_fields = merge_dictionaries( chosen_vep_fields, chosen_extra_vep_fields )
	for this_key in final_vep_fields:
		this_value = str(final_vep_fields[this_key])
		if (new_vep_info == ''):
			new_vep_info = str(this_key) + '=' + this_value
		else:
			new_vep_info = new_vep_info + ';' + str(this_key) + '=' + this_value
	if (new_vep_info == ''):
		new_vep_info = '.'

	return new_vep_info

######################################################
def build_new_info_with_reformatted_vep( this_info, vep_field_hdrs, extra_vep_field_hdrs ):

	new_info = '.'
	if (this_info != '.'):
		info_fields = this_info.split(';')
		for info_pair in info_fields:
			key_and_value = info_pair.split('=')
			new_vep_info = info_pair
			this_key = str(key_and_value[0])
			if (this_key == "CSQ"):
				idx = info_pair.find('=')
				old_vep_info = info_pair[idx+1:]
				new_vep_info = reformat_vep( old_vep_info, vep_field_hdrs, extra_vep_field_hdrs )
			if (new_info == '.'):
				new_info = new_vep_info
			else:
				new_info = new_info + ';' + new_vep_info

	return new_info

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Reformat the vep fields of a vcf record.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input VCF')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output VCF')
	args = parser.parse_args()

	outfile = open(args.outfile, 'w')

	# Open the input VCF and read its first record

	vep_field_hdrs = []
	extra_vep_field_hdrs = []
	infile = open(args.infile)
	for inline in infile:

		inline = inline.strip()
		is_header = False
		if (len(inline) >= 1):
			if (inline[0:1] == "#"):
				is_header = True
				is_vep_hdr = False
				is_chrom_hdr = False
				if (len(inline) >= 15):
					if (inline[0:15] == "##INFO=<ID=CSQ,"):
						is_vep_hdr = True
				if (is_vep_hdr): # Write new vep header instead of the old one.
					vep_field_hdrs, extra_vep_field_hdrs = write_new_vep_headers( inline, outfile )
				else: # Write out all other headers as is.
					outfile.write( inline + "\n" )

		if (is_header == False):

			infields = inline.split("\t")
			this_chrom = str(infields[0])
			this_pos = str(infields[1])
			this_id = str(infields[2])
			this_ref = str(infields[3])
			this_alt = str(infields[4])
			this_qual = str(infields[5])
			this_filter = str(infields[6])
			this_info = str(infields[7])
			this_format = ''
			this_samples = ''
			if (len(infields) >= 9):
				this_format = str(infields[8])
				if (len(infields) >= 10):
					this_samples = str(infields[9])
					if (len(infields) >= 11):
						for i in range( 10, len(infields) ):
							this_samples = this_samples + "\t" + str(infields[i])

			new_info = build_new_info_with_reformatted_vep( this_info, vep_field_hdrs, extra_vep_field_hdrs )

			outline = this_chrom + "\t" + this_pos + "\t" + this_id + "\t" + this_ref + "\t" + this_alt + "\t" + this_qual + "\t" + this_filter + "\t" + new_info
			if (this_format != ''):
				outline = outline + "\t" + this_format
			if (this_samples != ''):
				outline = outline + "\t" + this_samples
			outline = outline + "\n"
			outfile.write( outline )

	outfile.close()

if __name__=='__main__':
    main()


