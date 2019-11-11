#!/usr/bin/python
# python create_vcf_from_CosmicMutantExport_tsv_file.py -i infile -o outfile -r reference_genome_fasta_file -g file_of_gene_names_with_strand_info
# python create_vcf_from_CosmicMutantExport_tsv_file.py -i CosmicMutantExport_genesOfInterest.tsv -r ../data/reference_genomes/hs37d5x/hs37d5x.fa -g ../data/all_genes_GRCh37_exons.txt -o CosmicMutantExport_genesOfInterest.vcf

# Gene name	Accession Number	Gene CDS length	HGNC ID	Sample name	ID_sample	ID_tumour	Primary site	Site subtype 1	Site subtype 2	Site subtype 3	Primary histology	Histology subtype 1	Histology subtype 2	Histology subtype 3	Genome-wide screen	Mutation ID	Mutation CDS	Mutation AA	Mutation Description	Mutation zygosity	LOH	GRCh	Mutation genome position	Mutation strand	SNP	Resistance Mutation	FATHMM prediction	FATHMM score	Mutation somatic status	Pubmed_PMID	ID_STUDY	Sample source	Tumour origin	Age
# MYH4	ENST00000255381	5820	7574	PD4098a	1317010	1227889	breast	NS	NS	NS	carcinoma	lobular_carcinoma	NS	NS	y	COSM162723	c.5718G>C	p.E1906D	Substitution - Missense		u	37	17:10346794-10346794	-	n	-	PATHOGENIC	.94661	Confirmed somatic variant	22722201	385	NS	primary	43
# PMS2	ENST00000265849	2589	9122	TCGA-46-3769-01	1781825	1685824	lung	NS	NS	NS	carcinoma	squamous_cell_carcinoma	NS	NS	y	COSM747344	c.8G>A	p.R3Q	Substitution - Missense		u	37	7:6048643-6048643	-	n	-	NEUTRAL	.12424	Confirmed somatic variant		418	fresh/frozen - NOS	primary	57
# NF1	ENST00000358273	8520	7765	TCGA-GC-A3I6-01	1898129	1786813	urinary_tract	bladder	NS	NS	carcinoma	NS	NS	NS	y	COSM1302612	c.1339C>T	p.L447F	Substitution - Missense		u	37	17:29533336-29533336	+	n	-	PATHOGENIC	.96832	Variant of unknown origin		413	fresh/frozen - NOS	primary	45
# CDKN2A	ENST00000304494	471	1787	TCGA-G2-A2EO-01	1898125	1786809	urinary_tract	bladder	NS	NS	carcinoma	NS	NS	NS	y	COSM12484	c.322G>A	p.D108N	Substitution - Missense		u	37	9:21971036-21971036	-	n	-	PATHOGENIC	.9705	Reported in another cancer sample as somatic		413	fresh/frozen - NOS	primary	69
# PTEN	ENST00000371953	1212	9588	TCGA-A5-A0GH-01	1783310	1687309	endometrium	NS	NS	NS	carcinoma	endometrioid_carcinoma	NS	NS	y	COSM4898	c.950_953delT	ACT	p.T319fs*1	Deletion - Frameshift	het	u	37	10:89720799-89720802	+		-			Reported in another cancer sample as somatic		419	fresh/frozen - NOS	primary	57
# MYH3	ENST00000226209	5823	7573	TCGA-AP-A056-01	1783334	1687333	endometrium	NS	NS	NS	carcinoma	endometrioid_carcinoma	NS	NS	y	COSM975402	c.1350G>A	p.T450T	Substitution - coding silent	het	u	37	17:10547728-10547728	-	n	-	NEUTRAL	.01148	Variant of unknown origin		419	fresh/frozen - NOS	primary	64
# MYH8	ENST00000403437	5814	7578	TCGA-B5-A0JY-01	1783388	1687387	endometrium	NS	NS	NS	carcinoma	endometrioid_carcinoma	NS	NS	y	COSM975209	c.2295+2T>C	p.?	Unknown	het	u	37	17:10309589-10309589	-	n	-	PATHOGENIC	.99744	Variant of unknown origin		419	fresh/frozen - NOS	primary	50
# APC	ENST00000457016	8532	583	TCGA-AP-A0LM-01	1783352	1687351	endometrium	NS	NS	NS	carcinoma	endometrioid_carcinoma	NS	NS	y	COSM1059637	c.6757C>A	p.L2253I	Substitution - Missense	het	u	37	5:112178048-112178048	+	n	-	NEUTRAL	.38627	Variant of unknown origin		419	fresh/frozen - NOS	primary	33


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
def extract_chrom_and_pos( Mutation_genome_position ):

	# 10:26993598-26993598
	chrom = ''
	pos = ''
	pos_end = ''
	if (Mutation_genome_position != ''):
		bits = Mutation_genome_position.split(":")
		chrom = str(bits[0])
		bits2 = str(bits[1])
		bits3 = bits2.split("-")
		pos = str(bits3[0])
		pos_end = str(bits3[1])
	return chrom, pos, pos_end

######################################################
def determine_if_filter_out_this_entry( gene, chrom, pos ):

	filter_out_this_entry = False

	# if we are not given chrom or pos, then can't build a variant record with CHROM and POS
	if ((chrom == '') or (pos == '')):
		filter_out_this_entry = True

	return filter_out_this_entry

######################################################
def build_info_field( Gene_name, HGNC_ID, Mutation_ID, Mutation_CDS, Mutation_AA, Mutation_Description, FATHMM_prediction, FATHMM_score, pos_end ):

	return_info = 'COSMICMUTANTEXPORT_GENENAME=' + Gene_name + ';COSMICMUTANTEXPORT_HGNCID=' + HGNC_ID + ';COSMICMUTANTEXPORT_MUTATIONID=' + Mutation_ID + ';COSMICMUTANTEXPORT_MUTATIONCDS=' + Mutation_CDS + ';COSMICMUTANTEXPORT_MUTATIONAA=' + Mutation_AA + ';COSMICMUTANTEXPORT_MUTATIONDESCRIPTION=' + Mutation_Description + ';COSMICMUTANTEXPORT_FATHMMPREDICTION=' + FATHMM_prediction + ';COSMICMUTANTEXPORT_FATHMMSCORE=' + FATHMM_score

	if (pos_end != ''):
		return_info = return_info + ';COSMICMUTANTEXPORT_END=' + pos_end

	return_info = return_info.replace( ' ', '_' )

	return return_info

######################################################
def determine_whether_got_ref_and_alt_then_extract_them( Mutation_CDS ):

	ref = 'N'
	alt = '.'

	# Mutation_CDS:
	# c.163-8delT
	# c.163-8T>A
	# c.1980A>G

	bits = Mutation_CDS[2:]
	bits2 = bits.replace( '-', '' )
	bits2 = bits2.replace( '_', '' )
	bits2 = bits2.replace( '+', '' )
	bits2 = bits2.replace( '(', '' )
	bits2 = bits2.replace( ')', '' )
	bits3 = bits2.lstrip('0123456789')
	bits4 = bits3
	if (bits4.find('ins') > -1):
		bits4 = bits4.replace( 'ins', '' )
		bits4 = 'N>' + bits4
	if (bits4.find('del') > -1):
		bits4 = bits4.replace( 'del', '' )
		bits4 = bits4 + '>'
	bits5 = bits4.strip('0123456789?')
	bits6 = bits5.split('>')
	if (len(bits6) == 2):
		ref = bits6[0]
		alt = bits6[1]
	if (ref == ''):
		ref = 'N'
	if (alt == ''):
		alt = '.'

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
def is_ref_on_minus_strand( ref, chrom, pos, strand, reference_file ):

	# In some cases, a minus strand gene's REF has the plus strand, in other cases it is the minus strand.
	# To know which it is (so as to know whether to change REF and ALT to plus strand), compare it to the reference.

	is_ref_on_minus_strand = False

	if (strand == '-'):

		pos2 = int(pos) + len(ref) - 1
		command_output_lines = []
		samtools_faidx_position = chrom + ':' + str(pos) + '-' + str(pos2)
		samtools_faidx_command = 'samtools faidx ' + reference_file + ' ' + samtools_faidx_position
		command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
		if (command_status != 0):
			if (str(command_status) == '256'):
				ignore_this_warning = True
			else:
				raise ValueError("\n\nWas not able to get the possible minus strand reference sequence from reference genome for this region of the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more sequences.\n")
		else:
			command_output_lines = parse_out_any_warnings( command_output, samtools_faidx_command )
		# some positions may not be found in reference genome, will have a header but no sequence
		if (len(command_output_lines) > 1):
			retrieved_ref = command_output_lines[1]
			revcompl_retrieved_ref = reverse_compl(retrieved_ref)
			if (revcompl_retrieved_ref == ref):
				is_ref_on_minus_strand = True

	return is_ref_on_minus_strand

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
def get_true_seq_instead_of_N( ref, alt, chrom, pos, Mutation_CDS, reference_file ):

	# Examples where this function is called:
	# c.1509_1510insGCCTAT    p.Y503_F504insAY        Insertion - In frame            y       37      4:55592185-55592186
	# c.2235_2249del15        p.E746_A750delELREA     Deletion - In frame     het     u       37      7:55242465-55242479

	orig_ref = ref

	idx = ref.find('N')
	if (idx > -1):
		pos2 = (int(pos) + len(ref) - 1)
		command_output_lines = []
		samtools_faidx_position = chrom + ':' + str(pos) + '-' + str(pos2)
		samtools_faidx_command = 'samtools faidx ' + reference_file + ' ' + samtools_faidx_position
		command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
		if (command_status != 0):
			if (str(command_status) == '256'):
				print 'Warning, no reference sequence found for', samtools_faidx_position
			else:
				raise ValueError("\n\nWas not able to get the reference sequence from reference genome for this region of the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more sequences.\n")
		else:
			command_output_lines = parse_out_any_warnings( command_output, samtools_faidx_command )
		# some positions may not be found in reference genome, will have a header but no sequence
		if (len(command_output_lines) > 1):
			ref = command_output_lines[1]
			ref1 = ref[0:1]
			if (alt != ''):
				if (alt[0:1] == 'N'):
					# An example where this function is called:
					# c.2235_2249del15        p.E746_A750delELREA     Deletion - In frame     het     u       37      7:55242465-55242479
					alt1 = alt.replace( 'N', '' )
					alt = ref1 + alt1

					idx2 = Mutation_CDS.find('del')
					if (idx2 > -1):
						bits = Mutation_CDS.split('del')
						deletion_length = bits[1]
						if (is_integer(deletion_length)):
							pos2 = (int(pos) + int(deletion_length) - 1)
							command_output_lines = []
							samtools_faidx_position_2 = chrom + ':' + str(pos) + '-' + str(pos2)
							samtools_faidx_command_2 = 'samtools faidx ' + reference_file + ' ' + samtools_faidx_position_2
							command_status_2, command_output_2 = commands.getstatusoutput( samtools_faidx_command_2 )
							if (command_status_2 != 0):
								raise ValueError("\n\nWas not able to get the reference sequence from reference genome for this deletion region of " + str(Mutation_CDS) + " the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more sequences.\n")
							else:
								command_output_lines_2 = parse_out_any_warnings( command_output_2, samtools_faidx_command_2 )
							# some positions may not be found in reference genome, will have a header but no sequence
							if (len(command_output_lines_2) > 1):
								ref = command_output_lines_2[1]

				else:
					# An example where this function is called:
					# c.1509_1510insGCCTAT    p.Y503_F504insAY        Insertion - In frame            y       37      4:55592185-55592186
					alt = ref1 + alt
		else:
			print 'No reference sequence found for', samtools_faidx_position

	# Note, this function returns the first nucleotide for reference sequence on the positive strand.
	# The final result will be correct if the ALT sequence also has the positive strand values or no values.

	return ref, alt

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in a Cosmic tab-delimited file of variants, output the variants as a VCF file.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input CosmicMutantExport.tsv file')
	parser.add_argument('-r', action="store", dest="reference", required=True, help='Input reference fasta file')
	parser.add_argument('-g', action="store", dest="genes", required=True, help='Input genes file containing strand in 5th column')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output VCF file')

	args = parser.parse_args()

	gene_strands = read_gene_strands( args.genes )

	# output the VCF file headers

	outfile = open(args.outfile, 'w')
	outfile.write('##fileformat=VCFv4.2' + "\n")
	today = str(datetime.date.today())
	today = today.replace('-','')
	outfile.write('##fileDate=' + today + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_GENENAME,Number=1,Type=String,Description="CosmicMutantExport Gene_name">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_HGNCID,Number=1,Type=String,Description="CosmicMutantExport HGNC_ID">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_MUTATIONID,Number=1,Type=String,Description="CosmicMutantExport Mutation_ID">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_MUTATIONCDS,Number=1,Type=String,Description="CosmicMutantExport Mutation_CDS">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_MUTATIONAA,Number=1,Type=String,Description="CosmicMutantExport Mutation_AA">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_MUTATIONDESCRIPTION,Number=1,Type=String,Description="CosmicMutantExport Mutation_Description">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_FATHMMPREDICTION,Number=1,Type=String,Description="CosmicMutantExport FATHMM_prediction">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_FATHMMSCORE,Number=1,Type=String,Description="CosmicMutantExport FATHMM_score">' + "\n")
	outfile.write('##INFO=<ID=COSMICMUTANTEXPORT_END,Number=1,Type=String,Description="CosmicMutantExport genomic position end from Mutation_genome_position">' + "\n")
	outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + "\n")

	# read input file and write out each record

	infile = open(args.infile, 'r')
	is_header = True
	for inline in infile:

		inline = inline.strip()
		if (inline != ''):

			infields = inline.split("\t")
			Gene_name = str(infields[0])

			if (Gene_name == 'Gene name'):
				is_header = True
			else:
				is_header = False

			if (is_header == False):

				infields = inline.split("\t")

				Gene_name = str(infields[0])
				Accession_Number = str(infields[1])
				Gene_CDS_length = str(infields[2])
				HGNC_ID = str(infields[3])
				Sample_name = str(infields[4])
				ID_sample = str(infields[5])
				ID_tumour = str(infields[6])
				Primary_site = str(infields[7])
				Site_subtype_1 = str(infields[8])
				Site_subtype_2 = str(infields[9])
				Site_subtype_3 = str(infields[10])
				Primary_histology = str(infields[11])
				Histology_subtype_1 = str(infields[12])
				Histology_subtype_2 = str(infields[13])
				Histology_subtype_3 = str(infields[14])
				Genome_wide_screen = str(infields[15])
				Mutation_ID = str(infields[16])
				Mutation_CDS = str(infields[17])
				Mutation_AA = str(infields[18])
				Mutation_Description = str(infields[19])
				Mutation_zygosity = str(infields[20])
				LOH = str(infields[21])
				GRCh = str(infields[22])
				Mutation_genome_position = str(infields[23])
				Mutation_strand = str(infields[24])
				SNP = str(infields[25])
				Resistance_Mutation = str(infields[26])
				FATHMM_prediction = str(infields[27])
				FATHMM_score = str(infields[28])
				# Mutation_somatic_status = str(infields[29])
				# Pubmed_PMID = str(infields[30])
				# ID_STUDY = str(infields[31])
				# Sample_source = str(infields[32])
				# Tumour_origin = str(infields[33])
				# Age = str(infields[34])

				chrom, pos, pos_end = extract_chrom_and_pos( Mutation_genome_position )
				filter_out_this_entry = determine_if_filter_out_this_entry( Gene_name, chrom, pos )

				if ((filter_out_this_entry == False) and (pos != '')):

					info = build_info_field( Gene_name, HGNC_ID, Mutation_ID, Mutation_CDS, Mutation_AA, Mutation_Description, FATHMM_prediction, FATHMM_score, pos_end )

					vcfid = '.'
					qual = '.'
					vcffilter = '.'

					ref, alt = determine_whether_got_ref_and_alt_then_extract_them( Mutation_CDS )

					ref, alt = get_true_seq_instead_of_N( ref, alt, chrom, pos, Mutation_CDS, args.reference )

					ref_is_on_minus_strand = is_ref_on_minus_strand( ref, chrom, pos, Mutation_strand, args.reference )
					if (ref_is_on_minus_strand):
						ref = reverse_compl( ref )
						alt = reverse_compl( alt )

					outline = chrom + "\t" + pos + "\t" + vcfid + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + vcffilter + "\t" + info + "\n"
					outfile.write( outline )

	outfile.close()

if __name__=='__main__':
    main()


