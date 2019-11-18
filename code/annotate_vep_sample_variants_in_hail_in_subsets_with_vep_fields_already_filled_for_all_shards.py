#!/usr/bin/python
# python annotate_vep_sample_variants_in_hail_in_subsets_with_vep_fields_already_filled_for_all_shards.py -invcf_list invcf [-inregions inregions_file] [-inregion inregion] -cadd caddfile -gerp gerpfile -gnomad gnomadfile [-size subset_size_in_lines_of_inregions_file] [-log log_file] [-start start_from_line_in_inregions_file]

# python annotate_vep_sample_variants_in_hail_in_subsets_with_vep_fields_already_filled_for_all_shards.py -invcf_list list_shards_for_hail.txt -cadd_list list_cadd_shards.txt -gerp_list list_gerp_shards.txt -gnomad_list list_gnomad_shards.txt

# Read in a vcf file of a cohort containing certain genomic regions and sample, already annotated by vep . Read in a list of regions. Load them into hail.
# For those regions, annotate the cohort vds with other annotations (vep already there), and the provided CADD file.
# Output tab-delimited file, one line per variant per cohort sample, with annotation columns.
# If specified, break file into subsets to process separately in hail.
# If specified, start processing from specified middle of regions file (for when previous subsets are already processed).

# The input is already annotated with VEP:
# /nvme/emmrat/software_various/vep_nov2017/ensembl-vep/vep --vcf --offline --cache -dir_cache /nvme/emmrat/software_various/vep_nov2017/vep_data -i ../whole_exome_extract/02_ISKS_MGRB_HC.shard0062.chrom9.pos1-39688687.whole_exome.vcf.gz -o ../whole_exome_extract_vep/02_ISKS_MGRB_HC.shard0062.chrom9.pos1-39688687.whole_exome.vep.vcf --everything --force_overwrite
# python convert_vcf_vep_info_into_standard_info_fields.py -i ../whole_exome_extract_vep/02_ISKS_MGRB_HC.shard0062.chrom9.pos1-39688687.whole_exome.vep.vcf -o ../whole_exome_extract_vep_reformatted/02_ISKS_MGRB_HC.shard0062.chrom9.pos1-39688687.whole_exome.vep_reformatted.vcf

# Example of inregions file:
# 1:2161175-2161226
# 1:2234365-2234594
# 1:2234672-2234891
# Y:1534369-1535454
# Y:1535476-1535527
# Y:1605762-1605813

# head list_gerp_shards.txt 
# shard0001	/nvme/emmrat/ISKS_2018jan/hs37d5x/gerp_db/gerp.shard0001.chrom1.pos1-13077999.vcf.bgz
# shard0002	/nvme/emmrat/ISKS_2018jan/hs37d5x/gerp_db/gerp.shard0002.chrom1.pos13078000-17150659.vcf.bgz
# shard0003	/nvme/emmrat/ISKS_2018jan/hs37d5x/gerp_db/gerp.shard0003.chrom1.pos17150660-29953083.vcf.bgz


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import commands
import subprocess
import argparse
import re
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pandas.plotting
import numpy as np
import seaborn
from math import log, isnan
from pprint import pprint
from hail import *
from hail.expr import TString, TBoolean

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def read_in_whole_genome_annotation_files( hc ):

	### annotations from GNOMAD
	print 'read in annotations from GNOMAD'
	#vds_gnomad = hc.read(['/nvme/emmrat/ISKS_2018jan/hs37d5x/gnomad_db/gnomad.genomes.r2.0.1.sites.combined.split.minrep.vds']) # gnomad v2.0.1
	vds_gnomad = hc.import_vcf( '/nvme/emmrat/ISKS_2018jan/hs37d5x/gnomad_db_r2_1_1_20190823/gnomad.genomes.r2.1.1.sites.vcf.bgz' ) # gnomad v2.1.1

	### annotations from Clinvar
	print 'read in annotations from Clinvar'
	vds_clinvar = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/clinvar_20190718/clinvar_20190718.vds')
	vds_clinvarotherallele = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/clinvar_20190718/clinvar_other_allele_is_pathogenic_20190718.vds')

	### annotations from COSMIC
	print 'annotations from COSMIC'
	CosmicCodingMuts_vds = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/cosmic_db/CosmicCodingMuts.vds')
	CosmicMutantExport_vds = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/cosmic_db/CosmicMutantExport.vds')

	### annotations from Eigen
	print 'read in annotations from Eigen'
	vds_eigen = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/Eigen_hail_db/Eigen_coding_04092016.vds') # .repartition(10000)
	vds_eigen = vds_eigen.annotate_variants_expr('''
	    va.predictions = {Eigen: {
		EigenRaw: va.info.EigenRaw,
		EigenPhred: va.info.EigenPhred,
		EigenPCRaw: va.info.EigenPCRaw,
		EigenPCPhred: va.info.EigenPCPhred
	    }}''')

	### annotations from Condel
	print 'read in annotations from Condel'
	#vds_condel = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/condel_db/condel_scores.vds')
	vds_condel = hc.import_vcf( '/nvme/emmrat/ISKS_2018jan/hs37d5x/condel_db/condel_scores.vcf' )

	### annotations from Cato
	print 'read in annotations from Cato'
	vds_cato = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/Eigen_hail_db/CATO_1.1.vds')

	### annotations from SweGen
	print 'read in annotations from SweGen'
	vds_swegen = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/swegen_hail_db/swegen_autosomes_allelefreqs.vds')

	### annotations from Revel
	print 'read in annotations from Revel'
	vds_revel = hc.read('/nvme/emmrat/ISKS_2018jan/hs37d5x/revel_20190711/revel_20190711.vcf.vds')

	return vds_clinvar, vds_clinvarotherallele, CosmicCodingMuts_vds, CosmicMutantExport_vds, vds_eigen, vds_condel, vds_cato, vds_swegen, vds_revel, vds_gnomad

######################################################
def read_in_shard_annotation_files( hc, caddfile, gerpfile ): # gnomadfile ):

	### annotations from CADD
	print 'read in annotations from CADD', caddfile 
	#vds_cadd = hc.read(caddfile)
	vds_cadd = hc.import_vcf( caddfile )

	### annotations from GERP
	print 'read in annotations from GERP', gerpfile
	#vds_gerp = hc.read(gerpfile)
	vds_gerp = hc.import_vcf( gerpfile )

	### annotations from GNOMAD
	# print 'read in annotations from GNOMAD', gnomadfile
	# vds_gnomad = hc.import_vcf( gnomadfile )

	return vds_cadd, vds_gerp # vds_gnomad

######################################################
def extract_intervals_of_interests_in_cohort_vds( hc, cohort_vds, orig_inlines ):

	print 'extract_intervals_of_interests_in_cohort_vds'

	cohort_geneset_vds = cohort_vds

	if (len(orig_inlines) > 0):

		inlines = []
		for i in range( 0, len(orig_inlines) ):
			this_region = orig_inlines[i]
			bits1 = this_region.split(':')
			chrom = str(bits1[0])
			pos1_pos2 = str(bits1[1])
			bits2 = pos1_pos2.split('-')
			pos1 = int(bits2[0])
			pos2 = int(bits2[1])
			output_this_region = True
			inlines.append( this_region )

		inlines = orig_inlines

		intervals_all_regions_of_interest = map(Interval.parse, inlines)
		cohort_geneset_vds = cohort_vds.filter_intervals(intervals_all_regions_of_interest)

	return cohort_geneset_vds

######################################################
def annotate_and_write_output( this_subset, hc, cohort_geneset_vds, outtsv, outvds, vds_gnomad, vds_clinvar, vds_clinvarotherallele, CosmicCodingMuts_vds, CosmicMutantExport_vds, vds_eigen, vds_condel, vds_cadd, vds_gerp, vds_cato, vds_swegen, vds_revel, gnomad_af_nfe, gnomad_af ):

	print 'annotate_and_write_output for subset', this_subset

	# pprint(cohort_geneset_vds.variant_schema)

	### annotate and filter with GNOMAD
	print 'annotate and filter with GNOMAD for subset', this_subset

	# There are 3 places in this program that must be changed according to whether gnomad 2.0.1 or 2.1.1 is being used,
	# because our two dataset do not have the same name fields.

	# These are the fields for gnomad r2.0.1
	# cohort_geneset_gnomad_vds = cohort_geneset_vds.annotate_variants_vds(vds_gnomad, expr='va.gnomad = vds.gnomad') # gnomad v2.0.1
	# cohort_geneset_gnomad_NFElt0p05_vds = cohort_geneset_gnomad_vds.filter_variants_expr( '(va.gnomad.AF_NFE < 0.001) || isMissing(va.gnomad.AF_NFE)' ) # gnomad v2.0.1

	# These are the fields for gnomad r2.1.1
	cohort_geneset_gnomad_vds = cohort_geneset_vds.annotate_variants_vds(vds_gnomad, expr='va.gnomad = vds.info') # gnomad v2.1.1
	# cohort_geneset_gnomad_NFElt0p05_vds = cohort_geneset_gnomad_vds.filter_variants_expr( '(va.gnomad.AF_nfe[0] < 0.001) || isMissing(va.gnomad.AF_nfe[0])' ) # gnomad v2.1.1
	cohort_geneset_gnomad_NFElt0p05_vds = cohort_geneset_gnomad_vds.filter_variants_expr( '( ((va.gnomad.AF_nfe[0] < ' + str(gnomad_af_nfe) + ') || isMissing(va.gnomad.AF_nfe[0])) && ((va.gnomad.AF[0] < ' + str(gnomad_af) + ') || isMissing(va.gnomad.AF[0])) )' ) # gnomad v2.1.1
	# vds_gnomad_df = cohort_geneset_gnomad_NFElt0p05_vds.variants_table().to_pandas()
	# vds_gnomad_df[1:5]	

	### the input is already annotated with VEP
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_vds = cohort_geneset_gnomad_NFElt0p05_vds
	# For these genes, take only vep nonsense and missense variants
	# print 'For these genes, take only vep nonsense and missense variants for subset', this_subset
	# cohort_geneset_gnomad_NFElt0p05_vep_filtered_vds = cohort_geneset_gnomad_NFElt0p05_vds.filter_variants_expr( '(va.info.VEP_Consequence == "splice_acceptor_variant") || (va.info.VEP_Consequence == "splice_donor_variant") || (va.info.VEP_Consequence == "stop_gained") || (va.info.VEP_Consequence == "frameshift_variant") || (va.info.VEP_Consequence == "stop_lost") || (va.info.VEP_Consequence == "start_lost") || (va.info.VEP_Consequence == "inframe_insertion") || (va.info.VEP_Consequence == "inframe_deletion") || (va.info.VEP_Consequence == "missense_variant") || (va.info.VEP_Consequence == "incomplete_terminal_codon_variant") || (va.info.VEP_Consequence == "coding_sequence_variant")' )

	### annotate with Clinvar
	print 'annotate with Clinvar for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_vds.annotate_variants_vds(vds_clinvar, expr='va.clinvar = vds.info')
	print 'annotate with clinvar_other_allele_is_pathogenic for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_clinvarotherallele_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_vds.annotate_variants_vds(vds_clinvarotherallele, expr='va.clinvarotherallele = vds.info')

	### annotate with COSMIC
	print 'annotate with COSMIC for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_clinvarotherallele_vds.annotate_variants_vds(CosmicCodingMuts_vds, expr='va.cosmic.CosmicCodingMuts = vds.info')
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_vds.annotate_variants_vds(CosmicMutantExport_vds, expr='va.cosmic.CosmicMutantExport = vds.info')

	### annotate with Eigen
	print 'annotate with Eigen for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_vds.annotate_variants_vds(vds_eigen, expr='va.eigen = vds.info')

	### annotate with Condel
	print 'annotate with Condel for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_vds.annotate_variants_vds(vds_condel, expr='va.condel = vds.info')

	### annotate with CADD
	print 'annotate with CADD for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_vds.annotate_variants_vds(vds_cadd, expr='va.cadd = vds.info')

	### annotate with GERP
	print 'annotate with GERP for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_vds.annotate_variants_vds(vds_gerp, expr='va.gerp = vds.info')

	### annotate with Cato
	print 'annotate with Cato for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_vds.annotate_variants_vds(vds_cato, expr='va.cato = vds.info')

	### annotate with Swegen
	print 'annotate with Swegen for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_vds.annotate_variants_vds(vds_swegen, expr='va.swegen = vds.info')

	### annotate with Revel
	print 'annotate with Revel for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_revel_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_vds.annotate_variants_vds(vds_revel, expr='va.revel = vds.info')

	### annotate with cohort variant statistics
	print 'annotate with cohort variant statistics for subset', this_subset
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_revel_varqc_vds = cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_revel_vds.variant_qc()

	### export cohort variants to tab-delimited file
	print 'export cohort variants to tab-delimited file for subset', this_subset
	# when GT=1, visual_GT=0/1, gtj(g.gt)=0, gtk(g.gt)=1, g.ad[0]=depth_of_ref_allele, g.ad[1]=depth_of_alt_allele
	# when GT=2, visual_GT=1/1, gtj(g.gt)=1, gtk(g.gt)=1, g.ad[0]=depth_of_ref_allele, g.ad[1]=depth_of_alt_allele

	# These are the fields for gnomad r2.0.1
	# cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_revel_varqc_vds.export_genotypes(outtsv, 'SAMPLE=s, VARIANT=v, GT=g.gt, VAF=(g.ad[1] / g.dp), depth_of_alt=g.ad[1], depth_of_ref=g.ad[0], GQ=g.gq, DP=g.dp, gene_symbol=va.info.VEP_SYMBOL, gene_symbol_2=va.info.VEP_2_SYMBOL, vep_consequence=va.info.VEP_Consequence, vep_consequence_2=va.info.VEP_2_Consequence, vep_impact=va.info.VEP_IMPACT, vep_impact_2=va.info.VEP_2_IMPACT, vep_canonical=va.info.VEP_CANONICAL, vep_canonical_2=va.info.VEP_2_CANONICAL, vep_hgvsp=va.info.VEP_HGVSp, vep_transcript_id=va.info.VEP_Feature, vep_transcript_id_2=va.info.VEP_2_Feature, vep_polyphen=va.info.VEP_PolyPhen, vep_sift=va.info.VEP_SIFT, cohort_MAF=va.qc.AF, gnomad_AF=va.gnomad.AF, gnomad_AF_NFE=va.gnomad.AF_NFE, swegen_AF=va.swegen.swegen_AF, genotype_alt_allele_depth=g.ad[1], genotype_total_depth=g.dp, genotype_ad_div_sum_ad=g.ad[1]/sum(g.ad), genotype_ad_div_dp=g.ad[1]/g.dp, clinvar_Clinical_Significance=va.clinvar.CLINVAR_SIGNIFICANCE, clinvar_Clinical_Evidence=va.clinvar.CLINVAR_EVIDENCE, clinvar_Gene=va.clinvar.CLINVAR_GENE, clinvar_Name=va.clinvar.CLINVAR_NAME, clinvarOtherAllele_Clinical_Significance=va.clinvarotherallele.CLINVAROTHERALLELE_SIGNIFICANCE, clinvarOtherAllele_Clinical_Evidence=va.clinvarotherallele.CLINVAROTHERALLELE_EVIDENCE, clinvarOtherAllele_Gene=va.clinvarotherallele.CLINVAROTHERALLELE_GENE, clinvarOtherAllele_Name=va.clinvarotherallele.CLINVAROTHERALLELE_NAME, CosmicCodingMuts_Gene=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_GENE, CosmicCodingMuts_Strand=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_STRAND, CosmicCodingMuts_CDS=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_CDS, CosmicCodingMuts_AA=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_AA, CosmicCodingMuts_CNT=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_CNT, CosmicMutantExport_GeneName=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_GENENAME, CosmicMutantExport_HGNCId=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_HGNCID, CosmicMutantExport_MutationId=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONID, CosmicMutantExport_MutationsCDS=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONCDS, CosmicMutantExport_MutationAA=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONAA, CosmicMutantExport_MutationDescription=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONDESCRIPTION, CosmicMutantExport_FATHMMPrediction=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_FATHMMPREDICTION, CosmicMutantExport_FATHMMScore=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_FATHMMSCORE, CosmicMutantExport_End=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_END, EigenRaw=va.eigen.EigenRaw, EigenPhred=va.eigen.EigenPhred, EigenPCRaw=va.eigen.EigenPCRaw, EigenPCPhred=va.eigen.EigenPCPhred, CADD_RAW_SCORE=va.cadd.CADD_RAW_SCORE, CADD_PHRED=va.cadd.CADD_PHRED, GERP_SCORE=va.gerp.GERP_SCORE, GERP_RATE=va.gerp.GERP_RATE, Cato_score=va.cato.score, CONDEL_Score=va.condel.CONDEL_SCORE, CONDEL_FatHMM_Likely_Deleteriousness=va.condel.CONDEL_FATHMM, CONDEL_Mutation_Assessor=va.condel.CONDEL_MA, CONDEL_PolyPhen2=va.condel.CONDEL_PPH2, CONDEL_SIFT=va.condel.CONDEL_SIFT, CONDEL_Protein=va.condel.CONDEL_PROTEIN, CONDEL_Transcript=va.condel.CONDEL_TRANSCRIPT, CONDEL_Amino_Acid_Position=va.condel.CONDEL_AA_POS, CONDEL_Amino_Acid_REF=va.condel.CONDEL_AA_REF, CONDEL_Amino_Acid_ALT=va.condel.CONDEL_AA_ALT, CONDEL_Strand=va.condel.CONDEL_STRAND, vep_gene=va.info.VEP_Gene, vep_biotype=va.info.VEP_BIOTYPE, vep_exon=va.info.VEP_EXON, vep_intron=va.info.VEP_INTRON, vep_HGVSc=va.info.VEP_HGVSp, vep_HGVSc=va.info.VEP_HGVSp, vep_cDNA_position=va.info.VEP_cDNA_position, vep_CDS_position=va.info.VEP_CDS_position, vep_Protein_position=va.info.VEP_Protein_position, vep_Amino_acids=va.info.VEP_Amino_acids, vep_Codons=va.info.VEP_Codons, vep_Existing_variation=va.info.VEP_Existing_variation, vep_DISTANCE=va.info.VEP_DISTANCE, vep_STRAND=va.info.VEP_STRAND, vep_FLAGS=va.info.VEP_FLAGS, vep_VARIANT_CLASS=va.info.VEP_VARIANT_CLASS, vep_SYMBOL_SOURCE=va.info.VEP_SYMBOL_SOURCE, vep_HGNC_ID=va.info.VEP_HGNC_ID, vep_MANE=va.info.VEP_MANE, vep_TSL=va.info.VEP_TSL, vep_APPRIS=va.info.VEP_APPRIS, vep_CCDS=va.info.VEP_CCDS, vep_ENSP=va.info.VEP_ENSP, vep_SWISSPROT=va.info.VEP_SWISSPROT, vep_TREMBL=va.info.VEP_TREMBL, vep_UNIPARC=va.info.VEP_UNIPARC, vep_GENE_PHENO=va.info.VEP_GENE_PHENO, vep_DOMAINS=va.info.VEP_DOMAINS, vep_miRNA=va.info.VEP_miRNA, vep_HGVS_OFFSET=va.info.VEP_HGVS_OFFSET, vep_CLIN_SIG=va.info.VEP_CLIN_SIG, vep_SOMATIC=va.info.VEP_SOMATIC, vep_PHENO=va.info.VEP_PHENO, vep_PUBMED=va.info.VEP_PUBMED, vep_MOTIF_NAME=va.info.VEP_MOTIF_NAME, vep_MOTIF_POS=va.info.VEP_MOTIF_POS, vep_HIGH_INF_POS=va.info.VEP_HIGH_INF_POS, vep_MOTIF_SCORE_CHANGE=va.info.VEP_MOTIF_SCORE_CHANGE, vep_LoF=va.info.VEP_LoF, vep_LoF_filter=va.info.VEP_LoF_filter, vep_LoF_flags=va.info.VEP_LoF_flags, vep_LoF_info=va.info.VEP_LoF_info, vep_LoFtool=va.info.VEP_LoFtool, vep_REVEL=va.info.VEP_REVEL, vep_LOVD=va.info.VEP_LOVD, vep_CSN=va.info.VEP_CSN, vep_SpliceRegion=va.info.VEP_SpliceRegion, vep_TSSDistance=va.info.VEP_TSSDistance, revel_AAREF=va.revel.REVEL_AAREF, revel_AAALT=va.revel.REVEL_AAALT, revel_VALUE=va.revel.REVEL_VALUE') # gnomad v2.0.1

	# These are the fields for gnomad r2.1.1
	cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_revel_varqc_vds.export_genotypes(outtsv, 'SAMPLE=s, VARIANT=v, GT=g.gt, VAF=(g.ad[1] / g.dp), depth_of_alt=g.ad[1], depth_of_ref=g.ad[0], GQ=g.gq, DP=g.dp, gene_symbol=va.info.VEP_SYMBOL, gene_symbol_2=va.info.VEP_2_SYMBOL, vep_consequence=va.info.VEP_Consequence, vep_consequence_2=va.info.VEP_2_Consequence, vep_impact=va.info.VEP_IMPACT, vep_impact_2=va.info.VEP_2_IMPACT, vep_canonical=va.info.VEP_CANONICAL, vep_canonical_2=va.info.VEP_2_CANONICAL, vep_hgvsp=va.info.VEP_HGVSp, vep_transcript_id=va.info.VEP_Feature, vep_transcript_id_2=va.info.VEP_2_Feature, vep_polyphen=va.info.VEP_PolyPhen, vep_sift=va.info.VEP_SIFT, cohort_MAF=va.qc.AF, gnomad_AF=va.gnomad.AF, gnomad_AF_NFE=va.gnomad.AF_nfe, swegen_AF=va.swegen.swegen_AF, genotype_alt_allele_depth=g.ad[1], genotype_total_depth=g.dp, genotype_ad_div_sum_ad=g.ad[1]/sum(g.ad), genotype_ad_div_dp=g.ad[1]/g.dp, clinvar_Clinical_Significance=va.clinvar.CLINVAR_SIGNIFICANCE, clinvar_Clinical_Evidence=va.clinvar.CLINVAR_EVIDENCE, clinvar_Gene=va.clinvar.CLINVAR_GENE, clinvar_Name=va.clinvar.CLINVAR_NAME, clinvarOtherAllele_Clinical_Significance=va.clinvarotherallele.CLINVAROTHERALLELE_SIGNIFICANCE, clinvarOtherAllele_Clinical_Evidence=va.clinvarotherallele.CLINVAROTHERALLELE_EVIDENCE, clinvarOtherAllele_Gene=va.clinvarotherallele.CLINVAROTHERALLELE_GENE, clinvarOtherAllele_Name=va.clinvarotherallele.CLINVAROTHERALLELE_NAME, CosmicCodingMuts_Gene=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_GENE, CosmicCodingMuts_Strand=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_STRAND, CosmicCodingMuts_CDS=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_CDS, CosmicCodingMuts_AA=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_AA, CosmicCodingMuts_CNT=va.cosmic.CosmicCodingMuts.COSMICCODINGMUTS_CNT, CosmicMutantExport_GeneName=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_GENENAME, CosmicMutantExport_HGNCId=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_HGNCID, CosmicMutantExport_MutationId=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONID, CosmicMutantExport_MutationsCDS=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONCDS, CosmicMutantExport_MutationAA=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONAA, CosmicMutantExport_MutationDescription=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_MUTATIONDESCRIPTION, CosmicMutantExport_FATHMMPrediction=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_FATHMMPREDICTION, CosmicMutantExport_FATHMMScore=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_FATHMMSCORE, CosmicMutantExport_End=va.cosmic.CosmicMutantExport.COSMICMUTANTEXPORT_END, EigenRaw=va.eigen.EigenRaw, EigenPhred=va.eigen.EigenPhred, EigenPCRaw=va.eigen.EigenPCRaw, EigenPCPhred=va.eigen.EigenPCPhred, CADD_RAW_SCORE=va.cadd.CADD_RAW_SCORE, CADD_PHRED=va.cadd.CADD_PHRED, GERP_SCORE=va.gerp.GERP_SCORE, GERP_RATE=va.gerp.GERP_RATE, Cato_score=va.cato.score, CONDEL_Score=va.condel.CONDEL_SCORE, CONDEL_FatHMM_Likely_Deleteriousness=va.condel.CONDEL_FATHMM, CONDEL_Mutation_Assessor=va.condel.CONDEL_MA, CONDEL_PolyPhen2=va.condel.CONDEL_PPH2, CONDEL_SIFT=va.condel.CONDEL_SIFT, CONDEL_Protein=va.condel.CONDEL_PROTEIN, CONDEL_Transcript=va.condel.CONDEL_TRANSCRIPT, CONDEL_Amino_Acid_Position=va.condel.CONDEL_AA_POS, CONDEL_Amino_Acid_REF=va.condel.CONDEL_AA_REF, CONDEL_Amino_Acid_ALT=va.condel.CONDEL_AA_ALT, CONDEL_Strand=va.condel.CONDEL_STRAND, vep_gene=va.info.VEP_Gene, vep_biotype=va.info.VEP_BIOTYPE, vep_exon=va.info.VEP_EXON, vep_intron=va.info.VEP_INTRON, vep_HGVSc=va.info.VEP_HGVSp, vep_HGVSc=va.info.VEP_HGVSp, vep_cDNA_position=va.info.VEP_cDNA_position, vep_CDS_position=va.info.VEP_CDS_position, vep_Protein_position=va.info.VEP_Protein_position, vep_Amino_acids=va.info.VEP_Amino_acids, vep_Codons=va.info.VEP_Codons, vep_Existing_variation=va.info.VEP_Existing_variation, vep_DISTANCE=va.info.VEP_DISTANCE, vep_STRAND=va.info.VEP_STRAND, vep_FLAGS=va.info.VEP_FLAGS, vep_VARIANT_CLASS=va.info.VEP_VARIANT_CLASS, vep_SYMBOL_SOURCE=va.info.VEP_SYMBOL_SOURCE, vep_HGNC_ID=va.info.VEP_HGNC_ID, vep_MANE=va.info.VEP_MANE, vep_TSL=va.info.VEP_TSL, vep_APPRIS=va.info.VEP_APPRIS, vep_CCDS=va.info.VEP_CCDS, vep_ENSP=va.info.VEP_ENSP, vep_SWISSPROT=va.info.VEP_SWISSPROT, vep_TREMBL=va.info.VEP_TREMBL, vep_UNIPARC=va.info.VEP_UNIPARC, vep_GENE_PHENO=va.info.VEP_GENE_PHENO, vep_DOMAINS=va.info.VEP_DOMAINS, vep_miRNA=va.info.VEP_miRNA, vep_HGVS_OFFSET=va.info.VEP_HGVS_OFFSET, vep_CLIN_SIG=va.info.VEP_CLIN_SIG, vep_SOMATIC=va.info.VEP_SOMATIC, vep_PHENO=va.info.VEP_PHENO, vep_PUBMED=va.info.VEP_PUBMED, vep_MOTIF_NAME=va.info.VEP_MOTIF_NAME, vep_MOTIF_POS=va.info.VEP_MOTIF_POS, vep_HIGH_INF_POS=va.info.VEP_HIGH_INF_POS, vep_MOTIF_SCORE_CHANGE=va.info.VEP_MOTIF_SCORE_CHANGE, vep_LoF=va.info.VEP_LoF, vep_LoF_filter=va.info.VEP_LoF_filter, vep_LoF_flags=va.info.VEP_LoF_flags, vep_LoF_info=va.info.VEP_LoF_info, vep_LoFtool=va.info.VEP_LoFtool, vep_REVEL=va.info.VEP_REVEL, vep_LOVD=va.info.VEP_LOVD, vep_CSN=va.info.VEP_CSN, vep_SpliceRegion=va.info.VEP_SpliceRegion, vep_TSSDistance=va.info.VEP_TSSDistance, revel_AAREF=va.revel.REVEL_AAREF, revel_AAALT=va.revel.REVEL_AAALT, revel_VALUE=va.revel.REVEL_VALUE') # gnomad v2.1.1

	outtsv_done = outtsv + '.done'
	subprocess.call(['touch', outtsv_done])

	### export this cohort variants vds file in case we need to rerun use this in future
	# print 'export this cohort variants vds file in case we need to rerun use this in future for subset', this_subset
	# cohort_geneset_gnomad_NFElt0p05_vep_filtered_clinvar_CosmicCodingMuts_CosmicMutantExport_eigen_condel_cadd_gerp_cato_swegen_revel_varqc_vds.write( outvds, overwrite=True)

	return

######################################################
def log_this_subset_done( logfile, this_subset, num_subsets, subset_size, start_idx, end_idx, outtsv, outvds, inregions, inregion, invcf ):

	print 'log_this_subset_done'
	print_start_idx = start_idx + 1 # change from 0-based to 1-based line numbers
	print_end_idx = end_idx + 1 # change from 0-based to 1-based line numbers
	outline = 'Processed subset ' + str(this_subset) + ' of ' + str(num_subsets) + ' subsets. '
	outline = outline + 'Lines ' + str(print_start_idx) + ' to ' + str(print_end_idx) + ' of regions file ' + str(inregions) + ' '
	outline = outline + 'Or else inregion is ' + str(inregion) + ' '
	outline = outline + 'outtsv=' + str(outtsv) + ' outvds=' + str(outvds) + ' invcf=' + str(invcf)
	logfile_handle = open(logfile, 'a')
	logfile_handle.write( outline + "\n" )
	logfile_handle.close()

	return

######################################################
def main():

	# what input arguments have been supplied to the program

	print 'Process the arguments to this program.'
	parser = argparse.ArgumentParser(description='Read in vcf file. Load into hail and annotate. Output as tab-delimited, one line per variant per sample.')
	parser.add_argument('-invcf_list', action="store", dest="invcf_list", required=True, help='Tab-delimited input list, each line contains shard_id and cohort vds shard file, and outuput shard prefix. A tab-delimited file named outprefix.tsv will be produced by this program.')
	parser.add_argument('-inregions', action="store", dest="inregions", required=False, help='Input list of genomic regions')
	parser.add_argument('-inregion', action="store", dest="inregion", required=False, help='Input genomic region') # X:1-1000000
	parser.add_argument('-cadd_list', action="store", dest="caddfile_list", required=True, help='Tab-delimited input list, each line contains shard_id and CADD reference shard file.')
	parser.add_argument('-gerp_list', action="store", dest="gerpfile_list", required=True, help='Tab-delimited input list, each line contains shard_id and GERP reference shard file.')
	# parser.add_argument('-gnomad_list', action="store", dest="gnomadfile_list", required=True, help='Tab-delimited input list, each line contains shard_id and GNOMAD reference shard file.')
	parser.add_argument('-log', action="store", dest="logfile", required=False, help='Name of program log file to record subsets processed (default is temp_python_hail_logfile.txt')
	parser.add_argument('-size', action="store", dest="subset_size", required=False, help='Process in subsets of this many lines in the inregions file')
	parser.add_argument('-start', action="store", dest="startrec", required=False, help='Start processing from this variant record in the inregions files')
	parser.add_argument('-gnomad_af_nfe', action="store", dest="gnomad_af_nfe", required=False, help='Filter variants having gnomad_af_nfe less than this (defaults to 1 which means not filtered by this field)')
	parser.add_argument('-gnomad_af', action="store", dest="gnomad_af", required=False, help='Filter variants having gnomad_af less than this (defaults to 1 which means not filtered by this field)')
	args = parser.parse_args()

	inregions = ''
	inregion = ''
	invcf_list = str(args.invcf_list)
	caddfile_list = str(args.caddfile_list)
	gerpfile_list = str(args.gerpfile_list)
	# gnomadfile_list = str(args.gnomadfile_list)
	gnomad_af_nfe = float(1.1)
	if (args.gnomad_af_nfe is not None):
		gnomad_af_nfe = float(args.gnomad_af_nfe)
	gnomad_af = float(1.1)
	if (args.gnomad_af is not None):
		gnomad_af = float(args.gnomad_af)
	subset_size = 0
	if (args.subset_size is not None):
		subset_size = int(args.subset_size)
	startrec = 0
	if (args.startrec is not None):
		startrec = int(args.startrec) - 1 # change from 1-based to 0-based line numbers
	logfile = 'temp_python_hail_logfile.txt'
	if (args.logfile is not None):
		logfile = str(args.logfile)
	open(logfile, 'w').close() # initialise log file to be empty

	# initialise hail context

	hc = HailContext( tmp_dir='/nvme/emmrat/tmp_hail', log='/nvme/emmrat/tmp_hail_log/hail.log' )
	# hc.stop()
	# hc = HailContext(tmp_dir='/nvme/emmrat/tmp_hail')

	# read in the annotation files

	print 'Read in the whole genome annotation files.'
	vds_clinvar, vds_clinvarotherallele, CosmicCodingMuts_vds, CosmicMutantExport_vds, vds_eigen, vds_condel, vds_cato, vds_swegen, vds_revel, vds_gnomad = read_in_whole_genome_annotation_files( hc )
	#vds_clinvar, vds_clinvarotherallele, CosmicCodingMuts_vds, CosmicMutantExport_vds, vds_eigen, vds_condel, vds_cato, vds_swegen, vds_revel = read_in_whole_genome_annotation_files( hc )

	# Get the list of vcf shards to process

	in_shards = []
	in_vcfs = {}
	out_prefixes = {}
	invcf_list_file = open(invcf_list, 'r')
	for inline in invcf_list_file:
		inline = inline.strip()
		infields = inline.split("\t")
		this_shard = str(infields[0])
		this_in_vcf = str(infields[1])
		this_out_prefix = str(infields[2])
		in_shards.append( this_shard )
		in_vcfs[this_shard] = this_in_vcf
		out_prefixes[this_shard] = this_out_prefix

	caddfiles = {}
	caddfile_list_file = open(caddfile_list, 'r')
	for inline in caddfile_list_file:
		inline = inline.strip()
		infields = inline.split("\t")
		this_shard = str(infields[0])
		this_caddfile = str(infields[1])
		caddfiles[this_shard] = this_caddfile

	gerpfiles = {}
	gerpfile_list_file = open(gerpfile_list, 'r')
	for inline in gerpfile_list_file:
		inline = inline.strip()
		infields = inline.split("\t")
		this_shard = str(infields[0])
		this_gerpfile = str(infields[1])
		gerpfiles[this_shard] = this_gerpfile

	#gnomadfiles = {}
	#gnomadfile_list_file = open(gnomadfile_list, 'r')
	#for inline in gnomadfile_list_file:
	#	inline = inline.strip()
	#	infields = inline.split("\t")
	#	this_shard = str(infields[0])
	#	this_gnomadfile = str(infields[1])
	#	gnomadfiles[this_shard] = this_gnomadfile

	for i in range( 0, len(in_shards) ):

		this_shard = in_shards[i]

		invcf = in_vcfs[this_shard]
		outprefix = out_prefixes[this_shard]
		caddfile = caddfiles[this_shard]
		gerpfile = gerpfiles[this_shard]
		#gnomadfile = gnomadfiles[this_shard]

		outprefix_start = outprefix + '.start'
		subprocess.call(['touch', outprefix_start])

		print 'Processing shard', invcf
		print 'to produce files', outprefix

		print 'Read in the sharded annotation files.'
		vds_cadd, vds_gerp = read_in_shard_annotation_files( hc, caddfile, gerpfile )
		#vds_cadd, vds_gerp, vds_gnomad = read_in_shard_annotation_files( hc, caddfile, gerpfile, gnomadfile )

		# read in the cohort vds file of variants to be annotated

		print 'Read in the cohort vds file of variants to be annotated.', invcf
		cohort_vds = hc.import_vcf( invcf )
		# cohort_vds.summarize().report()
		# pprint(cohort_vds.variant_schema)
		# cohort_vds.count()
		cohort_split_vds = cohort_vds.split_multi()
		cohort_split_minrep_vds = cohort_split_vds.min_rep()
		cohort_vds = cohort_split_minrep_vds

		# read in the inregions file

		print 'Read in the file of regions to be annotated.'
		all_inlines = []
		if (args.inregions is not None):
			inregions = str(args.inregions)
			all_inlines = []
			infile = open(inregions, 'r')
			inlines_of_this_infile = infile.readlines()
			for i in range( 0, len(inlines_of_this_infile) ):
				this_inline = inlines_of_this_infile[i]
				this_inline = this_inline.strip()
				all_inlines.append(this_inline)
		else:
			if (args.inregion is not None):
				all_inlines.append(args.inregion)

		# create subsets of input and carry out processing for each subset
		if (subset_size > 0):
			print 'Create subsets of input and carry out processing for each subset.'
			num_subsets = float(len(all_inlines)) / float(subset_size)
			num_subsets = int(math.ceil(num_subsets))
			for i in range( 0, num_subsets ):
				this_subset = i + 1
				inlines = []
				start_idx = i * subset_size
				end_idx = start_idx + subset_size
				if (end_idx > len(all_inlines)):
					end_idx = len(all_inlines)
				if (startrec < end_idx):
					if (startrec > start_idx):
						start_idx = startrec
					if (startrec < end_idx):
						for j in range( start_idx, end_idx ):
							inlines.append( all_inlines[j] )
						outtsv = outprefix + '_SUBSET_' + str(this_subset) + '.tsv'
						outvds = outprefix + '_SUBSET_' + str(this_subset) + '.vds'
						print 'Processing subset', this_subset
						cohort_geneset_vds = extract_intervals_of_interests_in_cohort_vds( hc, cohort_vds, inlines )
						annotate_and_write_output( this_subset, hc, cohort_geneset_vds, outtsv, outvds, vds_gnomad, vds_clinvar, vds_clinvarotherallele, CosmicCodingMuts_vds, CosmicMutantExport_vds, vds_eigen, vds_condel, vds_cadd, vds_gerp, vds_cato, vds_swegen, vds_revel, gnomad_af_nfe, gnomad_af )
						log_this_subset_done( logfile, this_subset, num_subsets, subset_size, start_idx, (end_idx-1), outtsv, outvds, inregions, inregion, invcf )
			
		else: # process the input all in one go
			print 'Process the input all in one go.'
			inlines = all_inlines
			outtsv = outprefix + '.tsv'
			outvds = outprefix + '.vds'
			this_subset = 1
			cohort_geneset_vds = extract_intervals_of_interests_in_cohort_vds( hc, cohort_vds, inlines )
			# print 'cohort_vds'
			# cohort_vds.summarize().report()
			# print 'cohort_geneset_vds'
			# cohort_geneset_vds.summarize().report()
			annotate_and_write_output( this_subset, hc, cohort_geneset_vds, outtsv, outvds, vds_gnomad, vds_clinvar, vds_clinvarotherallele, CosmicCodingMuts_vds, CosmicMutantExport_vds, vds_eigen, vds_condel, vds_cadd, vds_gerp, vds_cato, vds_swegen, vds_revel, gnomad_af_nfe, gnomad_af )
			log_this_subset_done( logfile, 1, 1, "no subsetting", 0, len(inlines)-1, outtsv, outvds, inregions, inregion, invcf )

		print 'Finished processing!'

if __name__=='__main__':
    main()

