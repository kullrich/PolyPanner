######################################################################################################
# refine assembly
######################################################################################################

REFINE_CONTIGS?=$(OUTPUT_DIR)/contigs
REFINE_SEGMENTS?=$(OUTPUT_DIR)/segments
REFINE_BREAKPOINTS?=$(OUTPUT_DIR)/breakpoints

refine:
	bin//polypanner refine \
	        -ifn $(POP_TABLE) \
	        -fdr 0.25 \
	        -threads 80 \
	        -p_value 0.01 \
	        -pseudo_count 0.1 \
	        -margin 10 \
	        -min_contig_length 500 \
	        -min_length 200 \
	        -max_length 400 \
	        -max_lib_count 0 \
	        -cand_step 50 \
	        -only_dangles F \
	        -contig NA \
	        -ofn_contigs $(REFINE_CONTIGS) \
	        -ofn_segments $(REFINE_SEGMENTS) \
	        -ofn_breakpoints $(REFINE_BREAKPOINTS)

