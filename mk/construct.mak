######################################################################################################
# constructing PP files
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
LIB_POP?=$(OUTPUT_LIB_DIR)/lib.pop

# parse and pair SAM files 
parse_sam:
	mkdir -p $(OUTPUT_LIB_DIR)
	perl utils/parse_bwa_sam.pl $(SAM_R1) $(READS_R1) $(STATS_R1)
	perl utils/parse_bwa_sam.pl $(SAM_R2) $(READS_R2) $(STATS_R2)

pair_reads:
	perl utils/pair_reads.pl \
		$(READS_R1) \
		$(READS_R2) \
		match_length \
		$(READS_PAIRED) \
		$(READS_NON_PAIRED_R1) \
		$(READS_NON_PAIRED_R2) \
		$(PAIR_STATS)

construct:
	bin/polypanner construct \
		-ifn_paired $(READS_PAIRED) \
		-ifn_R1 $(READS_NON_PAIRED_R1) \
		-ifn_R2 $(READS_NON_PAIRED_R2) \
		-contig_table $(CONTIG_TABLE) \
	        -discard_clipped discard_both \
	        -min_score 30 \
	        -min_length 50 \
	        -max_edit 20 \
	        -trim 20 \
	        -unique_only T \
	        -ofn $(LIB_POP)

# construct polypanner files
construct_all:
	$(foreach S,$(SAMPLES),$(MAKE) parse_sam pair_reads construct SAMPLE=$S && ) true
