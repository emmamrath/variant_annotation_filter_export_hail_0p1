Variant annotation filter export hail 0p1
=========================================

## Overview
This is the variant pipeline used to annotate and filter variants for downstream analyses, 
using [VEP](https://asia.ensembl.org/info/docs/tools/vep/script/index.html) and [Hail 0.1](https://hail.is/docs/0.1/index.html), 
for the project to compare the MGRB and ISKS cohorts.  

This pipeline outputs a tab-delimited file of sample variants, one line per sample per variant, with the various annotations of each variant. 
The tab-delimited output is suitable for a clinician to use as a spreadsheet. 
The annotations are carried out for all variants in the input. 
However, the provided annotation reference files mainly contain coding exon regions and thus this pipeline annotates only coding exon regions.
Included code for annotations and data for annotations are for Clinvar, Cosmic, Cadd, Cato, Eigen, and Revel. 
Included code for annotations for which data must be obtained separately, due to their data being too large for github, are Gnomad, Condel, and Swegen.
The pipeline filters output to produce only variants whose Gnomad NFE_AF < 0.001 (or have no Gnomad 2.1.1 entry).
This filter is hard-coded in annotate_vep_sample_variants_in_hail_in_subsets_with_vep_fields_already_filled_for_all_shards.py and can be changed.  

This pipeline consists of bash scripts and python2 programs. 
It assumes that [VEP](https://asia.ensembl.org/info/docs/tools/vep/script/index.html) and some of its plugins, 
and [HAIL 0.1](https://hail.is/docs/0.1/index.html) are already installed. 
Input parameters to the scripts include the locations of input files, 
and locations of the [VEP](https://asia.ensembl.org/info/docs/tools/vep/script/index.html) and [HAIL](https://hail.is/docs/0.1/index.html) installations.
Although Hail 0.1 has a VEP annotation function, it crashes for regions close to the NANS gene, including TRIM14 gene.
Thus this pipeline runs VEP prior to running Hail. 
This architecture allows the latest version to VEP to be used, instead of the older version in Hail 0.1. 
This pipeline also uses [bcftools](http://samtools.github.io/bcftools/bcftools.html).  

Variants in this project are the germline genetic mutation SNVs recorded in a multi-sample [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file.  
* MGRB [(Medical Genome Reference Bank)](https://sgc.garvan.org.au/initiatives) data 
is a whole-genome data resource of 4000 healthy elderly individuals 
([manuscript 1](https://www.biorxiv.org/content/10.1101/473348v1), [manuscript 2](https://www.nature.com/articles/s41431-018-0279-z)).  
* ISKS [(International Sarcoma Kindred Study)](http://sarcomahelp.org/articles/sarcoma-kindred-study.html) data 
is a whole-genome data resource (manuscript in preparation).  

## Input data
The input vcf files of this project were broken into 145 subsets of genomic regions, called "shard" files. 
Thus this project assumes that it must process multiple input sharded vcfs.  

## Reference data
For this project, [VEP](https://asia.ensembl.org/info/docs/tools/vep/script/index.html) with various plugins is used to carry out many annotations. 
Thus this project assumes that a [VEP](https://asia.ensembl.org/info/docs/tools/vep/script/index.html) installation is available.  
This project then uses the [Hail 0.1](https://hail.is/docs/0.1/index.html) platform to load in additional annotation datasets for further annotations.  

## Temporary intermediate data
This pipeline creates intermediate files and does not remove them, 
so that should a subsequent processing step fail, processing can continue using the previous output's data files.  
After these files are used in the next step after the step where they were created, they are not used again and can be manually deleted.  

## Flowchart of processing
/my/cohort/before_shards/\*.vcf &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ==>  
```
annotate_vep_sample_variants_in_hail_00_01_extract_shards.sh
```
==> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/shard_data/\*.shard\*.chrom\*.pos\*.vcf &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ==>  
```
annotate_vep_sample_variants_in_hail_00_02_extract_exons_from_shards.sh  
```
==> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/shard_data_extract/\*.shard\*.chrom\*.pos\*.vcf &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ==>  
```
annotate_vep_sample_variants_in_hail_00_03_extract_shards.sh
```
==> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split/my_cohort.shard0001.chrom1.pos1-13077999.split.vcf  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split/my_cohort.shard0002.chrom1.pos13078000-17150659.split.vcf  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split/my_cohort.shard0003.chrom1.pos17150660-29953083.split.vcf &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ==>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ...  
```
annotate_vep_sample_variants_in_hail_01_run_vep.sh
```
==> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split_sorted/\*.shard\*.chrom\*.pos\*.split.sorted.vcf  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split_sorted_vep/\*.shard\*.chrom\*.pos\*.split.sorted.vep.vcf  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split_sorted_vep_reformatted/\*.shard\*.chrom\*.pos\*.split.sorted.vep_reformatted.vcf &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ==>  
```
annotate_vep_sample_variants_in_hail_02_run_hail.sh
```
==> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split_sorted_vep_reformatted_hail/\*.shard\*.chrom\*.pos\*.split.sorted.vep_reformatted_hail.tsv &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ==>  
```
annotate_vep_sample_variants_in_hail_02_run_hail.sh
```
==> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /my/cohort/data_extract_split_sorted_vep_reformatted_hail_extra/\*.shard\*.chrom\*.pos\*.split.sorted.vep_reformatted_hail.extra.tsv

## Input parameter for pipeline steps
The input parameters for calling the scripts in this pipeline can be seen by looking at the beginning of each script.

### annotate_vep_sample_variants_in_hail_00_01_extract_shards.sh
indir=$1 # /my/cohort/before_shards  
outdir=$2 # /my/cohort/shard_data  
tmpdir_base=$3 # ./tmp  

### annotate_vep_sample_variants_in_hail_00_02_extract_exons_from_shards.sh
indir=$1 # "/my/cohort/shard_data"  
outdir=$2 # "/my/cohort/shard_data_extract"  
tmpdir_base=$3 # "./tmp"  

### annotate_vep_sample_variants_in_hail_00_03_extract_shards.sh
indir=$1 # "/my/cohort/shard_data_extract"  
outdir=$2 # "/my/cohort/data_extract_split"  
tmpdir_base=$3 # "./tmp"  

### annotate_vep_sample_variants_in_hail_01_run_vep.sh
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

### annotate_vep_sample_variants_in_hail_02_run_hail.sh
inlist=$1 # "list_shards_for_hail.txt"  
indir=$2 # "/my/cohort/data_extract_split_sorted_vep_reformatted" # this indir is specified in inlist file too  
outdir=$3 # "/my/cohort/data_extract_split_sorted_vep_reformatted_hail" # this outdir is specified in inlist file too  
caddlist=$4 # "list_cadd_shards.txt"  
gerplist=$5 # "list_gerp_shards.txt"  
gnomadlist=$6 # "list_gnomad_shards.txt"  

### annotate_vep_sample_variants_in_hail_03_extra_annotation.sh
shard1=$1 # "0001"  
shard2=$2 # "0145"  
indir=$3 # "/my/cohort/data_extract_split_sorted_vep_reformatted_hail"  
outdir=$4 # "/my/cohort/data_extract_split_sorted_vep_reformatted_hail_extra"  

## Reference data used by this pipeline that is provided with this pipeline
Clinvar: data came from [https://www.ncbi.nlm.nih.gov/clinvar/](https://www.ncbi.nlm.nih.gov/clinvar/) on 20190718  
Cosmic: data came from [https://cancer.sanger.ac.uk/cosmic/download](https://cancer.sanger.ac.uk/cosmic/download) on 2019 November  
Eigen: data came from [http://www.columbia.edu/~ii2135/download.html](http://www.columbia.edu/~ii2135/download.html) . Only coding regions data is included.  
Cato: data came from [http://www.mauranolab.org/CATO/](http://www.mauranolab.org/CATO/)  
Revel: data came from [https://sites.google.com/site/revelgenomics/downloads](https://sites.google.com/site/revelgenomics/downloads)  
Cadd: data came from [https://cadd.gs.washington.edu/](https://cadd.gs.washington.edu/)  

## Reference data used by this pipeline that needs to be download from public sites because it is not included
The following large reference data are not provided with this pipeline, and are used by the pipeline to annotate variants.  
They need to be downloaded as specific below.  

### Gnomad 
Download [gnomad.genomes.r2.1.1.sites.vcf.bgz](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz) and [gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi)  
from [https://gnomad.broadinstitute.org/downloads](https://gnomad.broadinstitute.org/downloads)  
into data/gnomad directory. File size is 461G.  

### Condel 
Download [http://bbglab.irbbarcelona.org/fannsdb/downloads/fannsdb.tsv.gz](http://bbglab.irbbarcelona.org/fannsdb/downloads/fannsdb.tsv.gz) into data/condel directory  
gunzip fannsdb.tsv.gz # File size is 22G.  
cd data/condel  
./create_condel_vcf.sh # creates data/condel/condel_scores.vcf from fannsdb.tsv  

### Swegen
Obtain the Swegen vcf from [https://swefreq.nbis.se/](https://swefreq.nbis.se/).  
Convert to hail 01. vds format by running the following in [hail 0.1](https://hail.is/docs/0.1/index.html) python:  
hc.import_vcf('../data/swegen/swegen_autosomes_allelefreqs.vcf').min_rep().write('../data/swegen/swegen_autosomes_allelefreqs.vds', overwrite=True)  

### Gerp
cd data/gerp  
wget http://mendel.stanford.edu/sidowlab/downloads/gerp/hg19.GERP_scores.tar.gz  
tar xzvf hg19.GERP_scores.tar.gz  
./create_gerp_shard_vcfs_01_create_chrom_positions.sh  
./create_gerp_shard_vcfs_02_create_shard_vcfs.sh  

### GRCh37 fasta file
cd data/reference_genome_hs37d5x  
Obtain the GRCh37 hs375dx fasta file and its indexes for GRCh37.  

## Creating reference data from a new version of the source reference data
The following reference data are provided with this pipeline, and are used by the pipeline to annotate variants.  
They come from publicly available databases that are regularly updated.  
If you want to use the latest data from these publicly available databases instead of the older versions provided with this pipeline,  then follow instructions below to obtain the latest version and prepare it for use by this pipeline.  

### Clinvar
cd data/clinvar  
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz  
gunzip variant_summary.txt.gz  
mv variant_summary.txt clinvar_variant_summary.txt  
cd ../../code  
./create_vcf_from_clinvar_file.sh  

### Cosmic
cd data/cosmic  
sftp -o User=<my_email@my_email.org> sftp-cancer.sanger.ac.uk  
get /files/grch37/cosmic/v84/CosmicMutantExport.tsv.gz  
get /files/grch37/cosmic/v84/CosmicResistanceMutations.tsv.gz  
get /files/grch37/cosmic/v84/VCF/CosmicCodingMuts.vcf.gz  
get /files/grch37/cosmic/v84/VCF/CosmicNonCodingVariants.vcf.gz  
get /files/grch37/cosmic/v84/CosmicHGNC.tsv.gz  
cd ../../code  
./create_vcf_files_from_cosmic_files.sh  
Note that the above script uses [tabix](http://www.htslib.org/doc/tabix.html), [bgzip](http://www.htslib.org/doc/bgzip.html), and [vt normalize](https://genome.sph.umich.edu/wiki/Vt).  


