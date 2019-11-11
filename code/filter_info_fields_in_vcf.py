#!/usr/bin/python
# python filter_info_fields_in_vcf.py -i infile -o outfile -f list_of_info_fields_to_keep
# python filter_info_fields_in_vcf.py -i CosmicCodingMuts_geneOfInterest_tooManyInfoFields.vcf -f 'COSMICCODINGMUTS_GENE,COSMICCODINGMUTS_STRAND,COSMICCODINGMUTS_CDS,COSMICCODINGMUTS_AA,COSMICCODINGMUTS_CNT' -o CosmicCodingMuts_geneOfInterest.vcf

# Input file is a VCF file. For the INFO field in each record, keep only the fields specified on input.


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
def extract_info_field_name_from_header_record( inline ):

	this_info_field_name = ''

	bits = inline.split('##INFO=<ID=')
	bit2 = bits[1]
	bits2 = bit2.split(',')
	this_info_field_name = bits2[0]

	return this_info_field_name

######################################################
def is_this_info_field_in_list_of_fields_to_keep( this_info_field_name, info_fields ):

	is_in_list = False

	if this_info_field_name in info_fields:
		is_in_list = True

	return is_in_list

######################################################
def keep_only_info_fields_to_keep( info, info_fields ):

	new_info = ''

	bits = info.split(';')
	for key_and_value in bits:
		bits2 = key_and_value.split('=')
		this_key = bits2[0] # to get value, don't split on '=' because INFO field value can contain '='. it's ok to split on '=' to get INFO key though.
		this_key_idx = key_and_value.find('=')
		if (this_key_idx > -1):
			this_value = key_and_value[ (this_key_idx+1): ]


		if this_key in info_fields:
			if (new_info == ''):
				new_info = str(key_and_value)
			else:
				new_info = new_info + ';' + str(key_and_value)

	if (new_info == ''):
		new_info = '.'

	return new_info

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Input file is a VCF file. For the INFO field in each record, keep only the fields specified on input.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output VCF file')
	parser.add_argument('-f', action="store", dest="fields", required=True, help='INFO field names, separated by commas')
	args = parser.parse_args()

	info_fields = {}
	info_fields_string = args.fields
	bits = info_fields_string.split(',')
	for this_field in bits:
		this_field = this_field.strip()
		info_fields[this_field] = this_field

	# output the VCF file headers

	outfile = open(args.outfile, 'w')

	# read input file and write out each record

	infile = open(args.infile, 'r')
	in_header = True
	for inline in infile:

		inline = inline.strip()
		if (inline != ''):

			if (in_header == True):
				output_this_line = True

				if (len(inline) > 11):
					if (inline[0:11] == '##INFO=<ID='):

						output_this_line = False
						this_info_field_name = extract_info_field_name_from_header_record( inline )
						is_in_list = is_this_info_field_in_list_of_fields_to_keep( this_info_field_name, info_fields )
						if (is_in_list):
							output_this_line = True

				if (output_this_line):
					outfile.write(inline + "\n")

				if (len(inline) >= 1):
					if (inline[0:1] == '#'):
						in_header = True
					else:
						in_header = False

			if (in_header == False):

				infields = inline.split("\t")
				chrom = str(infields[0])
				pos = str(infields[1])
				vcfid = str(infields[2])
				ref = str(infields[3])
				alt = str(infields[4])
				qual = str(infields[5])
				vcffilter = str(infields[6])
				info = str(infields[7])

				fields_after_info = ''
				if (len(infields) > 8):
					for i in range( 0, len(infields) ):
						fields_after_info = fields_after_info + "\t" + str(infields[i])

				new_info = keep_only_info_fields_to_keep( info, info_fields )

				outline = chrom + "\t" + pos + "\t" + vcfid + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + vcffilter + "\t" + new_info + fields_after_info + "\n"
				outfile.write( outline )

	outfile.close()

if __name__=='__main__':
    main()


