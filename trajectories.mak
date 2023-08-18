######################################################################################################
# output trajectories
######################################################################################################


# summed read counts of bins over samples
BIN_TRAJECTORY?=$(OUTPUT_DIR)/bin_trajectory

# variant and total reads counts for samples
SITE_TRAJECTORY_COUNT?=$(OUTPUT_DIR)/site_count
SITE_TRAJECTORY_TOTAL?=$(OUTPUT_DIR)/site_total

bin_trajectory:
	bin/polypanner bin_trajectory \
	        -ifn_libs $(POP_TABLE) \
	        -ifn_segments $(SEG_BIN_FINAL) \
	        -ofn $(BIN_TRAJECTORY)

site_trajectory:
	bin/polypanner site_trajectory \
	        -ifn_libs $(POP_TABLE) \
	        -ifn_sites $(BASE_SITES) \
	        -ifn_segments $(SEG_BIN_FINAL) \
	        -ofn_counts $(SITE_TRAJECTORY_COUNT) \
	        -ofn_totals $(SITE_TRAJECTORY_TOTAL)

