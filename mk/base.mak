######################################################################################################
# input files
######################################################################################################

ASSEMBLY_ID?=test
INPUT_DIR?=input/$(ASSEMBLY_ID)

# contig and length
CONTIG_TABLE?=$(INPUT_DIR)/contig_table

# contig fasta
CONTIG_FASTA?=$(INPUT_DIR)/contigs.fa

SAMPLES?=t1 t2 t3 t4
SAMPLE?=t1
INPUT_LIB_DIR?=$(INPUT_DIR)/$(SAMPLE)

# read sides mapped separarately using bwa into the SAM format
SAM_R1?=$(INPUT_LIB_DIR)/R1.sam
SAM_R2?=$(INPUT_LIB_DIR)/R2.sam

######################################################################################################
# basic output directories
######################################################################################################

# base output directory
OUTPUT_DIR?=output/$(ASSEMBLY_ID)

# library output dir
OUTPUT_LIB_DIR?=$(OUTPUT_DIR)/$(SAMPLE)

# table with paths to pop files
POP_TABLE?=$(OUTPUT_DIR)/libs.fns

######################################################################################################
# consrtuct.mak
######################################################################################################

# parse sam reads
READS_R1?=$(OUTPUT_LIB_DIR)/R1.tab
STATS_R1?=$(OUTPUT_LIB_DIR)/R1.stats
READS_R2?=$(OUTPUT_LIB_DIR)/R2.tab
STATS_R2?=$(OUTPUT_LIB_DIR)/R2.stats

# pair reads
READS_PAIRED?=$(OUTPUT_LIB_DIR)/paired
READS_NON_PAIRED_R1?=$(OUTPUT_LIB_DIR)/non_paired_R1
READS_NON_PAIRED_R2?=$(OUTPUT_LIB_DIR)/non_paired_R2
PAIR_STATS?=$(OUTPUT_LIB_DIR)/pair.stats

# PP file
LIB_POP?=$(OUTPUT_LIB_DIR)/lib.pop.gz

######################################################################################################
# seq_errors.mak
######################################################################################################

# merge all libraries
MERGE_POP?=$(OUTPUT_DIR)/merge.pop.gz

# all pop files
ALL_POPS=$(addsuffix /lib.pop,$(addprefix $(OUTPUT_DIR)/,$(SAMPLES)))

# segregating sites
BASE_SITES?=$(OUTPUT_DIR)/base_sites

# inferred error parameters
ERROR_PARAMS?=$(OUTPUT_DIR)/error_params

# pop file, resrticted to true segregating sites
RESTRICT_POP?=$(OUTPUT_LIB_DIR)/restricted.pop.gz

######################################################################################################
# refine.mak
######################################################################################################

REFINE_CONTIGS?=$(OUTPUT_DIR)/contigs
REFINE_SEGMENTS?=$(OUTPUT_DIR)/segments
REFINE_BREAKPOINTS?=$(OUTPUT_DIR)/breakpoints

######################################################################################################
# genomes.mak
######################################################################################################

# fasta for segments
SEGMENT_FASTA?=$(OUTPUT_DIR)/segments.fa

# coverage matrix for segments 
SEGMENT_COV_MATRIX?=$(OUTPUT_DIR)/segments.cov

# metaBAT output
SEGMENT_BIN_BASE?=$(OUTPUT_DIR)/seg_bin_base

# base bins
BASE_BINS?=$(OUTPUT_DIR)/base_bins

# associate segment with bin
SEG_BIN_FINAL?=$(OUTPUT_DIR)/seg_bin_final

######################################################################################################
# sites.mak
######################################################################################################

SITES_ANNOTATED?=$(OUTPUT_DIR)/sites_annotated

######################################################################################################
# trajectories.mak
######################################################################################################

# summed read counts of bins over samples
BIN_TRAJECTORY?=$(OUTPUT_DIR)/bin_trajectory

# variant and total reads counts for samples
SITE_TRAJECTORY_COUNT?=$(OUTPUT_DIR)/site_count
SITE_TRAJECTORY_TOTAL?=$(OUTPUT_DIR)/site_total

