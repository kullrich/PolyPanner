######################################################################################################
# identify polymorphic sites
######################################################################################################

SITES_ANNOTATED?=$(OUTPUT_DIR)/sites_annotated

sites:
	bin/polypanner sites \
	        -ifn_libs $(POP_TABLE) \
	        -ifn_sites $(BASE_SITES) \
	        -ifn_contigs $(CONTIG_TABLE) \
	        -ifn_segments $(REFINE_SEGMENTS) \
	        -ifn_segment_bins $(SEG_BIN_FINAL) \
	        -pseudo_count 0.1 \
	        -regional_window 2000 \
	        -margin 200 \
	        -ofn $(SITES_ANNOTATED)
