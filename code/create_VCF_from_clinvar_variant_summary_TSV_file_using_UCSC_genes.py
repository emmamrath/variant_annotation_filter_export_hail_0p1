#!/usr/bin/python
# python create_VCF_from_clinvar_variant_summary_TSV_file_using_UCSC_genes.py -i infile -o outfile -r reference_genome_fasta_file -g file_of_gene_names_with_strand_info
# python create_VCF_from_clinvar_variant_summary_TSV_file_using_UCSC_genes.py -i temp_variant_summary.txt -r /my/reference_genome_for_MGRB/hs37d5x/hs37d5x.fa -g /my/reference_data/UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_20170928.bed -o clinvar_variant_summary.vcf


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import commands
import argparse
import re

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def get_idx_of_last_digit_in_string(s):

	idx = -1
	new_str = s[::-1]
	match_str = re.search('\d', new_str)
	idx = len(new_str) - match_str.start(0) - 1
	return idx

######################################################
def compl(char):

	rtn_char = ''
	if (char == 'A'):
		rtn_char = 'T'
	elif (char == 'T'):
		rtn_char = 'A'
	elif (char == 'C'):
		rtn_char = 'G'
	elif (char == 'G'):
		rtn_char = 'C'
	elif (char == 'N'):
		rtn_char = 'N'
	elif (char == '.'):
		rtn_char = '.'
	else:
		print 'Unexpected nucleotide', char
	return rtn_char

######################################################
def reverse_compl(seq):

	rtn_seq = ''
	for i in range( len(seq)-1, -1, -1 ):
		char = seq[ i : i+1 ]
		char = compl( char )
		rtn_seq = rtn_seq + char
	return rtn_seq

######################################################
def multiply_nucleotide( char, num_times ):

	num_times = int(num_times)
	if (num_times > 100):
		num_times = 100
	rtn_seq = char * num_times
	return rtn_seq

######################################################
def get_minimum_2_idx( in1, in2, in3, in4 ):

	out1 = -1
	out2 = -1
	valid_idx = []
	if (in1 > -1):
		valid_idx.append( in1 )
	if (in2 > -1):
		valid_idx.append( in2 )
	if (in3 > -1):
		valid_idx.append( in3 )
	if (in4 > -1):
		valid_idx.append( in4 )
	if (len(valid_idx) > 0):
		out1 = min( valid_idx )
		valid_idx2 = []
		for this_idx in valid_idx:
			if (this_idx != out1):
				valid_idx2.append( this_idx )
		if (len(valid_idx2) > 0):
			out2 = min( valid_idx2 )
	return out1, out2

######################################################
def extract_pos( pos ):

	return_pos = pos.replace(" ", "")
	idx = return_pos.find('-')
	if (idx > -1):
		return_pos = return_pos[ 0 : idx ]
	return return_pos

######################################################
def extract_significance( significance ):

	return_significance = significance
	idx = return_significance.find('(')
	if (idx > -1):
		return_significance = return_significance[ 0 : idx ]
		return_significance = return_significance.replace( ' ', '_' )
	return return_significance

######################################################
def determine_if_filter_out_this_entry( gene, chrom, pos, assembly ):

	filter_out_this_entry = False

	# variant is a large deletion, gene = ALPL|RERE|C1QA|C1QB|C1QC|CA6|CAPZB|CASP9|TNFRSF8|CDA|CDK11B|CDC42|CLCN6|CLCNKA|CLCNKB|CNR2|CORT|DDOST|DFFA
	# idx = gene.find('|')
	# if (idx > -1):
	#	filter_out_this_entry = True

	if (assembly != 'GRCh37'):
		filter_out_this_entry = True

	# if we are not given chrom or pos, then can't build a variant record with CHROM and POS
	if ((chrom == '') or (pos == '')):
		filter_out_this_entry = True

	return filter_out_this_entry

######################################################
def determine_if_filter_out_this_entry_2( ref, alt ):

	filter_out_this_entry = False

	# Remove ref==alt so that it doesn't cause downstream problems in Hail.
	# ref==alt when we actually haven't properly recorded the ref or alt due to it being too large, eg.
	# NM_001001290.1(SLC2A9):c.727+8793C=
	# In that case, it is less likely we would match the SNP to a new sample anyway.
	# Clinvar has some ref==alt anyway, eg.
	# NM_001304561.1(BTNL2):c.1078A= (p.Ser360=)
	if (ref == alt):
		filter_out_this_entry = True

	# If we end up with rubbish, we should reject the record, eg.
	# NM_001918.3(DBT):c.1017_1018insNC_000001.11:g.100207187_100207312
	if (is_a_valid_ref_or_alt(ref) == False):
		filter_out_this_entry = True
	if (is_a_valid_ref_or_alt(alt) == False):
		filter_out_this_entry = True

	return filter_out_this_entry

######################################################
def build_info_field( significance, evidence, gene, pos, variant_name ):

	# significance:
	# Uncertain significance(Last reviewed: Sep 17, 2013)

	# gene:
	# SDHD

	# Hail will crash if it receives commas in the info field because it will try to create a tuple in a field that doesn't accept that?
	significance = significance.replace( ', ', '/' )
	evidence = evidence.replace( ', ', '/' )
	gene = gene.replace( ', ', '/' )
	variant_name = variant_name.replace( ', ', '/' )
	significance = significance.replace( ',', '/' )
	evidence = evidence.replace( ',', '/' )
	gene = gene.replace( ',', '/' )
	variant_name = variant_name.replace( ',', '/' )

	# Hail will crash if it receives spaces in the info field.
	significance = significance.replace( ' ', '_' )
	evidence = evidence.replace( ' ', '_' )
	gene = gene.replace( ' ', '_' )
	variant_name = variant_name.replace( ' ', '_' )

	return_info = 'CLINVAR_SIGNIFICANCE=' + significance + ';CLINVAR_EVIDENCE=' + evidence + ';CLINVAR_GENE=' + gene + ';CLINVAR_NAME=' + variant_name

	# GRCh37Location field:
	# 17354330
	# 17354337 - 17354339
	new_pos = pos.replace(" ", "")
	idx = new_pos.find('-')
	if (idx > -1):
		end_pos = str(new_pos[ idx+1 : ])
		return_info = return_info + ';END=' + end_pos

	return return_info

######################################################
def remove_trailing_protein_description_in_brackets( variant_name ):

	# 1946_1960delCATTTATTCCTAATG (p.Ala649_Asn653del)	==>	1946_1960delCATTTATTCCTAATG

	idx = variant_name.find(' (p.')
	if (idx > -1):
		variant_name = str(variant_name[ 0 : idx ])

	return variant_name

######################################################
def split_variant_name_into_positions_part_and_ins_del_desc_part( variant_name ):

	# 311delAinsGG			==>	311		and	delAinsGG	
	# 784_787dupGCTA		==>	784_787		and	dupGCTA
	# 1946_1960delCATTTATTCCTAATG	==>	1946_1960	and	delCATTTATTCCTAATG
	# *133T>C			==>	*133		and	T>C

	variant_positions = variant_name
	variant_trailing_desc = ''
	idx_for_del = variant_name.find('del')
	idx_for_ins = variant_name.find('ins')
	idx_for_dup = variant_name.find('dup')
	idx_for_inv = variant_name.find('inv')
	if ((idx_for_del > -1) or (idx_for_ins > -1) or (idx_for_dup > -1) or (idx_for_inv > -1)):
		list_of_idx = []
		if (idx_for_del > -1):
			list_of_idx.append( idx_for_del )
		if (idx_for_ins > -1):
			list_of_idx.append( idx_for_ins )
		if (idx_for_dup > -1):
			list_of_idx.append( idx_for_dup )
		if (idx_for_inv > -1):
			list_of_idx.append( idx_for_inv )
		idx = min( list_of_idx )
		variant_positions = variant_name[ 0 : idx ]
		variant_trailing_desc = variant_name[ idx : ]

	if (variant_trailing_desc == ''): # if we didn't find del, ins, dup, inv, etc.
		idx = variant_name.find('>')
		if (idx > -1):
			variant_positions = variant_name[ 0 : idx-1 ]
			variant_trailing_desc = variant_name[ idx-1 : ]

	if (variant_trailing_desc == ''): # if we didn't find G>T, del, ins, dup, inv, etc.
		idx = variant_name.find('=')
		if (idx > -1):
			variant_positions = variant_name[ 0 : idx-1 ]
			ref = variant_name[ idx-1 : idx ]
			variant_trailing_desc = ref + '>' + ref

	return variant_positions, variant_trailing_desc

######################################################
def extract_ref_and_alt_from_ref_alt_description( variant_trailing_desc, variant_name ):

	# G>A
	# delTTCTTC
	# delAinsGG
	# dupGGTC

	ref = 'N'
	alt = '.'

	idx = variant_trailing_desc.find('>')
	if (idx > -1): # is of type G>A
		ref = variant_trailing_desc[ 0 : idx ]
		alt = variant_trailing_desc[ idx+1 : ]

	else: # is not of type G>A
		idx_for_del = variant_trailing_desc.find('del')
		idx_for_ins = variant_trailing_desc.find('ins')
		idx_for_dup = variant_trailing_desc.find('dup')
		idx_for_inv = variant_trailing_desc.find('inv')
		if ((idx_for_del > -1) or (idx_for_ins > -1) or (idx_for_dup > -1) or (idx_for_inv > -1)): # is of type eg. delTTCTTC or delAinsGG

			if ((idx_for_ins > -1) and (idx_for_del == -1) and (idx_for_dup == -1) and (idx_for_inv == -1)): # insAA or insA
				ref = 'N' # we don't have the ref in this file
				alt = variant_trailing_desc[ idx_for_ins+3 : ]
				if (is_integer(alt)):
					alt = multiply_nucleotide( 'N', int(alt) )
				alt = 'N' + alt
			elif ((idx_for_del > -1) and (idx_for_ins == -1) and (idx_for_dup == -1) and (idx_for_inv == -1)): # delATT or del or delT
				ref = variant_trailing_desc[ idx_for_del+3 : ]
				alt = '.'
				if (ref == ''):
					ref = 'N' # we don't have the ref in this file
			elif ((idx_for_del > -1) and (idx_for_ins > -1) and (idx_for_dup == -1) and (idx_for_inv == -1)): # delATinsGC or delinsACCC
				if ((idx_for_ins < idx_for_del)):
					print 'Unexpected ins before del for:', variant_name
				else:
					ref = variant_trailing_desc[ idx_for_del+3 : idx_for_ins ]
					alt = variant_trailing_desc[ idx_for_ins+3 : ]
					if (ref == ''):
						ref = 'N'
			elif ((idx_for_dup > -1) and (idx_for_del == -1) and (idx_for_ins == -1) and (idx_for_inv == -1)): # dupGGTC
				ref = variant_trailing_desc[ idx_for_dup+3 : ]
				alt = ref + ref
			elif ((idx_for_ins > -1) and (idx_for_dup > -1) and (idx_for_del == -1) and (idx_for_inv == -1)): # dupinsGGTC
				if ((idx_for_ins < idx_for_dup)):
					print 'Unexpected ins before dup for:', variant_name
				else:
					ref = variant_trailing_desc[ idx_for_ins+3 : ]
					alt = ref + ref
			elif ((idx_for_inv > -1) and (idx_for_del == -1) and (idx_for_ins == -1) and (idx_for_dup == -1)): # invTTAA
				ref = variant_trailing_desc[ idx_for_inv+3 : ]
				alt = reverse_compl( ref )
			else:
				print 'Unexpected del,ins,dup,inv combination for:', variant_name

		else: # is not of type eg. is of type eg. delTTCTTC or delAinsGG
			print 'Did not find matching variant_trailing_desc for: ', variant_trailing_desc, 'in', variant_name

	if (is_integer(ref)):
		ref = multiply_nucleotide( 'N', int(ref) )
	alt = alt.upper()
	if (is_integer(alt)):
		alt = multiply_nucleotide( 'N', int(alt) )
	elif (alt == 'ALU'):
		alt = multiply_nucleotide( 'N', 100 )

	if (ref == ''):
		ref = 'N'
	if (alt == ''):
		alt = '.'

	return ref, alt

######################################################
def fix_typos( variant_name ):

	variant_name = variant_name.replace( 'deins', 'delins' )
	return variant_name

######################################################
def remove_dupl_chrom( chrom ):

	new_chrom = chrom
	idx = chrom.find('|')
	if (idx > -1):
		new_chrom = chrom[ 0 : idx ]
	return new_chrom

######################################################
def is_a_valid_ref_or_alt( ref_or_alt ):

	is_valid = True
	if (ref_or_alt.find('_') > -1):
		is_valid = False
	elif (ref_or_alt.find(':') > -1):
		is_valid = False
	elif (ref_or_alt.find('.') > -1):
		if (ref_or_alt != '.'):
			is_valid = False
	else:
		bit1 = ref_or_alt.replace( 'A', '' )
		bit2 = bit1.replace( 'C', '' )
		bit1 = bit2.replace( 'G', '' )
		bit2 = bit1.replace( 'T', '' )
		new_ref_or_alt = bit2
		if (len(new_ref_or_alt) > 0):
			is_valid = False

	# eg. of typing mistake which will cause Hail to crash with invalid nucleotide:
	# NM_139241.2(FGD4):c.1698G>H

	return is_valid

######################################################
def determine_whether_got_ref_and_alt_then_extract_them( variant_name ):

	ref = 'N'
	alt = '.'
	saved_variant_name = variant_name

	# variant_name field:
	# NM_003000.2(SDHB):c.*159_*184delinsGAACCTGTTCCTTTACTTGCCCCAA		= 3'-UTR, nucleotides 159 to 184
	#									for 5'-UTR, use -30 to -1
	# NM_000553.4(WRN):c.3234-130T>C
	# NM_003000.2(SDHB):c.*133T>C
	# NM_003001.3(SDHC):c.*1181_*1182insA
	# NM_003000.2(SDHB):c.814A>G (p.Thr272Ala)
	# NM_003000.2(SDHB):c.768T>C (p.Gly256=)
	# NM_003000.2(SDHB):c.17_42dup (p.Ala15Profs)
	# NM_003000.2(SDHB):c.210dupC (p.Met71Hisfs)
	# NM_003000.2(SDHB):c.784_787dupGCTA (p.Ile263Serfs)
	# NM_000251.2(MSH2):c.260_261insT (p.Val89Cysfs)
	# NM_003000.2(SDHB):c.761_762insC (p.Lys255Terfs)
	# NM_003000.2(SDHB):c.-151_*159del
	# NM_003000.2(SDHB):c.766-20T>C
	# NM_003000.2(SDHB):c.424-19_424-14delTTCTTC
	# NM_003000.2(SDHB):c.424-16_424-14dupTTC				= 16 nt in intron before exon nt 424, to 14 nt in intron before exon nt 424
	# NM_003000.2(SDHB):c.765+29G>A						NM_003000.2 = mRNA
	# NR_103459.1(SDHC):n.107+36T>A						NR_103459.1 = RNA
	# NM_003000.2(SDHB):c.713delT (p.Phe238Serfs)
	# NM_003000.2(SDHB):c.540+4_540+10delATTAGTC				= 4 nt in intron after exon nt 540, to 10 nt in intron after exon nt 540
	# NM_003000.2(SDHB):c.642_642+6delGGTGAGG
	# NM_003000.2(SDHB):c.311delAinsGG (p.Asn104Argfs)
	# NM_003000.2(SDHB):c.595_604delTACTGGTGGAinsGG (p.Tyr199Glyfs)
	# NM_003000.2(SDHB):c.*159_*184delinsGAACCTGTTCCTTTACTTGCCCCAA
	# NC_000001.10:g.161279434_161299373del19940
	# NC_000001.10:g.17375249_17390927del15679
	# NC_000001.11:g.(?_17044761)_(17044888_?)del				NC_000001.11:g = DNA
	# NC_000001.11:g.(?_17018722)_(17054170_?)dup
	# NC_000002.11:g.38121110_47669522inv
	# NM_000535.6(PMS2):c.1621A= (p.Lys541=)
	# NM_000059.3:c.9452_9466dupinsAlu
	# NM_000059.3:c.9342_9343insALU
	# NM_000267.3:c.7566_8314+3321_delinsCGCCACGGC
	# NM_000267.3:c.3975-1046_4110+902deins333
	# GRCh37/hg19 1p36.13(chr1:17366385-17440734)x3
	# GRCh37/hg19 1p36.33-q44(chr1:849467-249224684)
	# GRCh37/hg19 1p36.33-q44(chr1:849467-249224684)x3
	# GRCh37/hg19 1p36.32-36.12(chr1:2749920-22564787)x1
	# GRCh38/hg38 1q23.1-25.1(chr1:157747246-176021247)x3
	# GRCh37/hg19 9q22.32(97579146-99280739)x1
	# NM_000251.2(MSH2):c.1311_1334del24insNM_000251.1:c.1338_1361inv24 (p.Thr438_Ser445delinsPheSerLysPheGlnGluMetIle)

	found_pattern = False

	is_GRCh_pattern = False
	if (len(variant_name) > 4):
		char1to4 = variant_name[0:4]
		if (char1to4 == 'GRCh'):
			is_GRCh_pattern = True
	if (is_GRCh_pattern) : 			# GRCh37/hg19 1p36.32-36.12(chr1:2749920-22564787)x1
		#				# GRCh37/hg19 9q22.32(97579146-99280739)x1
		# simply use the genomic co-ordinates, they are extracted elsewhere, we don't have REF to extract
		found_pattern = True

	else: # is not like GRCh37/hg19 1p36.32-36.12(chr1:2749920-22564787)x1

		# see if it is like NM_003000.2(SDHB):c.765+29G>A
		idx = variant_name.find(':')
		if (idx > -1): # NM_003000.2(SDHB):c.*133T>C, NM_003000.2(SDHB):c.540+4_540+10delATTAGTC, NC_000001.10:g.17375249_17390927del15679, etc.
			variant_name = variant_name[ idx+1 : ]
			idx = variant_name.find('.')
			if (idx > -1):
				found_pattern = True
				variant_name = fix_typos( variant_name )
				variant_name = str(variant_name[ idx+1 : ])
				variant_name = remove_trailing_protein_description_in_brackets( variant_name )
				variant_positions, variant_trailing_desc = split_variant_name_into_positions_part_and_ins_del_desc_part( variant_name )
				ref, alt = extract_ref_and_alt_from_ref_alt_description( variant_trailing_desc, variant_name )

	if (found_pattern == False):
		print 'Did not match pattern for this variant_name:', saved_variant_name

	return ref, alt

######################################################
def read_gene_strands( genes_file_name ):

	gene_strands = {}
	genes_file = open(genes_file_name, 'r')
	for ingene in genes_file:
		infields = ingene.split("\t")
		gene_name = infields[3]
		gene_strand = infields[4]
		gene_strand = gene_strand.strip()
		gene_strands[ gene_name ] = gene_strand

	return gene_strands

######################################################
def reverse_compl_if_on_minus_strand( seq, gene, gene_strands ):

	print seq, gene, gene_strands
	if (gene in gene_strands):
		this_strand = gene_strands[ gene ]
		print this_strand
		if (this_strand == '-'):
			print 'here1'
			seq = reverse_compl( seq )
	# else, if we don't find the gene because it is multiple genes hit by this ClinVar variant,
	# then let's assume that the plus strand is what is reported and thus we don't need to reverse complement it.
	return seq

######################################################
def reverse_compl_if_ref_on_minus_strand( ref, alt, chrom, pos, reference_file ):

	if (ref != ""):
		ref1 = ref[0:1]
		pos2 = pos
		command_output_lines = []
		samtools_faidx_position = chrom + ':' + str(pos) + '-' + str(pos2)
		samtools_faidx_command = 'samtools faidx ' + reference_file + ' ' + samtools_faidx_position
		command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
		if (command_status != 0):
			raise ValueError("\n\nIn reverse_compl_if_ref_on_minus_strand: Was not able to get the reference sequence from reference genome for this region of the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more sequences.\n")
		else:
			command_output_lines = parse_out_any_warnings( command_output, samtools_faidx_command )
		# some positions may not be found in reference genome, will have a header but no sequence
		if (len(command_output_lines) > 1):
			new_ref = command_output_lines[1]
			new_ref1 = new_ref[0:1]
			compl_new_ref1 = compl(new_ref1)
			if (ref == compl_new_ref1):
				ref = reverse_compl(ref)
				alt = reverse_compl(alt)
		else:
			print 'No reference sequence found for', samtools_faidx_position

	return ref, alt

######################################################
def parse_out_any_warnings( command_output, command ):

	outlines = []
	all_outlines = command_output.split("\n")
	for outline in all_outlines:
		is_a_warning_line = False
		if (len(outline) >= 7):
			if (outline[0:7] == 'Warning'):
				is_a_warning_line = True
		if (is_a_warning_line == False):
			outlines.append( outline )

	return outlines

######################################################
def get_true_seq_instead_of_N( ref, alt, chrom, pos, reference_file ):

	idx = ref.find('N')
	if (idx > -1):
		pos2 = (int(pos) + len(ref) - 1)
		command_output_lines = []
		samtools_faidx_position = chrom + ':' + str(pos) + '-' + str(pos2)
		samtools_faidx_command = 'samtools faidx ' + reference_file + ' ' + samtools_faidx_position
		command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
		if (command_status != 0):
			raise ValueError("\n\nIn get_true_seq_instead_of_N: Was not able to get the reference sequence from reference genome for this region of the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more sequences.\n")
		else:
			command_output_lines = parse_out_any_warnings( command_output, samtools_faidx_command )
		# some positions may not be found in reference genome, will have a header but no sequence
		if (len(command_output_lines) > 1):
			ref = command_output_lines[1]
			ref1 = ref[0:1]
			if (alt != ''):
				if (alt[0:1] == 'N'):
					if (len(alt) == 1):
						alt = ref1
					else:
						alt = ref1 + alt[1:]
		else:
			print 'No reference sequence found for', samtools_faidx_position
	return ref, alt

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in a ClinVar tab-delimited file of variants, output the variants as a VCF file.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input ClinVar file')
	parser.add_argument('-r', action="store", dest="reference", required=True, help='Input reference fasta file')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output VCF file')
	parser.add_argument('-g', action="store", dest="genes", required=True, help='Input genes file, tab-delimited, format chrom pos1 pos2 gene strand')
	args = parser.parse_args()

	gene_strands = read_gene_strands( args.genes )

	# output the VCF file headers

	outfile = open(args.outfile, 'w')
	outfile.write('##fileformat=VCFv4.2' + "\n")
	today = str(datetime.date.today())
	today = today.replace('-','')
	outfile.write('##fileDate=' + today + "\n")
	outfile.write('##INFO=<ID=CLINVAR_SIGNIFICANCE,Number=1,Type=String,Description="ClinVar signficance">' + "\n")
	outfile.write('##INFO=<ID=CLINVAR_EVIDENCE,Number=1,Type=String,Description="ClinVar evidence">' + "\n")
	outfile.write('##INFO=<ID=CLINVAR_GENE,Number=1,Type=String,Description="ClinVar gene name">' + "\n")
	outfile.write('##INFO=<ID=CLINVAR_NAME,Number=1,Type=String,Description="ClinVar HGVSp- or HGVSc-like name">' + "\n")
	outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + "\n")

	# read input file and write out each record

	infile = open(args.infile, 'r')
	for inline in infile:

		inline = inline.strip()
		if (inline != ''):

			is_header = False
			if (len(inline) > 1):
				if (inline[0:1] == '#'):
					is_header = True

			if (is_header == False):

				infields = inline.split("\t")
				variant_name = str(infields[2])
				gene = str(infields[4])
				significance = str(infields[6])
				assembly = str(infields[16])
				evidence = str(infields[24])
				chrom = str(infields[18])
				location = str(infields[19])
				pos = extract_pos(location)
				significance = extract_significance(significance)
				filter_out_this_entry = determine_if_filter_out_this_entry( gene, chrom, pos, assembly )

				if ((filter_out_this_entry == False) and (pos != '')):

					chrom = remove_dupl_chrom( chrom )
					info = build_info_field( significance, evidence, gene, pos, variant_name )

					vcfid = '.'
					qual = '.'
					vcffilter = '.'

					ref, alt = determine_whether_got_ref_and_alt_then_extract_them( variant_name )

					ref, alt = get_true_seq_instead_of_N( ref, alt, chrom, pos, args.reference )

					#ref = reverse_compl_if_on_minus_strand( ref, gene, gene_strands )
					#alt = reverse_compl_if_on_minus_strand( alt, gene, gene_strands )
					ref, alt = reverse_compl_if_ref_on_minus_strand( ref, alt, chrom, pos, args.reference )

					filter_out_this_entry_2 = determine_if_filter_out_this_entry_2( ref, alt )
					if (filter_out_this_entry_2 == False):

						outline = chrom + "\t" + pos + "\t" + vcfid + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + vcffilter + "\t" + info + "\n"
						outfile.write( outline )

	outfile.close()

if __name__=='__main__':
    main()


