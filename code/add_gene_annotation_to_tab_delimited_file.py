#!/usr/bin/python
# python add_gene_annotation_to_tab_delimited_file.py -i infile -o outfile -num_aa gene_length_file
# python add_gene_annotation_to_tab_delimited_file.py -i ISKS.WholeGenome.shard0048.chrom7.pos257486-50390632.renamed_samples.vqsr.whole_exome.split.sorted.vep_reformatted_hail.tsv -o temp_shard0048_out.tsv -num_aa UCSC_tables_GRCh37_RefSeq_genes_20190716_canonical_transcripts_gene_length.txt

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
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Add a column to tab-delimited file. Match record\'s gene name to reference file to get new column\'s value.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input TSV')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output TSV')
	parser.add_argument('-num_aa', action="store", dest="num_aa", required=True, help='Reference file')
	args = parser.parse_args()

	outfile = open(args.outfile, 'w')

	num_canonical_aa = {}
	num_longest_aa = {}
	if (args.num_aa is not None):
		infile = open(args.num_aa, 'r')
		is_header = True
		for inline in infile:
			if (is_header):
				is_header = False # first record is header, so don't process
			else:
				inline = inline.strip()
				if (inline != ''):
					infields = inline.split("\t")
					this_gene_name = str(infields[0])
					this_num_canonical_aa = str(infields[1])
					this_num_longest_aa = str(infields[2])
					if (this_num_canonical_aa != 'NA'):
						this_num_canonical_aa = int(this_num_canonical_aa)
						num_canonical_aa[this_gene_name] = this_num_canonical_aa
					if (this_num_longest_aa != 'NA'):
						this_num_longest_aa = int(this_num_longest_aa)
						num_longest_aa[this_gene_name] = this_num_longest_aa

	# Open the input file, first record is header

	infile = open(args.infile, 'r')
	is_header = True
	for inline in infile:

		inline = inline.strip()
		if (len(inline) >= 1):

			if (is_header):

				new_col1 = "CANONICAL_GENE_LENGTH_IN_AMIMO_ACIDS"
				new_col2 = "LONGEST_GENE_LENGTH_IN_AMINO_ACIDS"
				inline = inline + "\t" + new_col1 + "\t" + new_col2
				outfile.write( inline + "\n" )
				is_header = False

			else:

				infields = inline.split("\t")
				this_gene_name = str(infields[8])
				this_num_canonical_aa = 'NOT_FOUND'
				this_num_longest_aa = 'NOT_FOUND'
				if this_gene_name in num_canonical_aa:
					this_num_canonical_aa = str(num_canonical_aa[this_gene_name])
				if this_gene_name in num_longest_aa:
					this_num_longest_aa = str(num_longest_aa[this_gene_name])
				inline = inline + "\t" + this_num_canonical_aa + "\t" + this_num_longest_aa
				outfile.write( inline + "\n" )

	outfile.close()

if __name__=='__main__':
    main()


