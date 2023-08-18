######################################################################################################
# filter seq errors
######################################################################################################

# merge all libraries
MERGE_POP?=$(OUTPUT_DIR)/merge.pop

# all pop files
ALL_POPS=$(addsuffix /lib.pop,$(addprefix $(OUTPUT_DIR)/,$(SAMPLES)))

# segregating sites
BASE_SITES?=$(OUTPUT_DIR)/base_sites

# inferred error parameters
ERROR_PARAMS?=$(OUTPUT_DIR)/error_params

# pop file, resrticted to true segregating sites
RESTRICT_POP?=$(OUTPUT_LIB_DIR)/restricted.pop

merge:
	bin/polypanner merge \
		$(MERGE_POP) \
		$(ALL_POPS)

filter:
	bin/polypanner filter \
	        -ifn $(MERGE_POP) \
	        -threads 16 \
	        -alpha 0.001 \
	        -min_variant_count 3 \
	        -min_contig_length 1000 \
	        -min_error_pct 0.001 \
	        -max_error_pct 5 \
	        -seed_error_pct 1 \
	        -ofn_sites $(BASE_SITES) \
	        -ofn_error_params $(ERROR_PARAMS)

restrict:
	bin/polypanner restrict \
	        -ifn_cav $(LIB_POP) \
	        -ifn_sites $(BASE_SITES) \
	        -ofn $(RESTRICT_POP)

restrict_all:
	$(foreach S,$(SAMPLES),$(MAKE) restrict SAMPLE=$S && ) true
