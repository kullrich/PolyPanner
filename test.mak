######################################################################################################
# input files
######################################################################################################

INPUT_DIR?=input

# contig and length
CONTIG_TABLE?=$(INPUT_DIR)/contig_table

# contig fasta
CONTIG_FASTA?=$(INPUT_DIR)/contigs.fa

# table with paths to pop files
POP_TABLE?=$(INPUT_DIR)/libs.fns

SAMPLES?=t1 t2 t3 t4
SAMPLE?=t1
INPUT_LIB_DIR?=$(INPUT_DIR)/$(SAMPLE)

# read sides mapped separarately using bwa into the SAM format
SAM_R1?=$(INPUT_LIB_DIR)/R1.sam
SAM_R2?=$(INPUT_LIB_DIR)/R2.sam

######################################################################################################
# output directories
######################################################################################################

# base output directory
OUTPUT_DIR?=output

# library output dir
OUTPUT_LIB_DIR?=output/$(SAMPLE)



