######################################################################################################
# genome binning
######################################################################################################

# generate coverage matrix 
cov_matrix:
	bin/polypanner cov_matrix \
	        -ifn_cavs $(POP_TABLE) \
	        -ifn_segs $(REFINE_SEGMENTS) \
	        -ifn_fasta $(CONTIG_FASTA) \
	        -actual_nts T \
	        -ofn_fasta $(SEGMENT_FASTA) \
	        -ofn_mat $(SEGMENT_COV_MATRIX)

# call metaBAT 
metaBAT:
	rm -rf /tmp/pp && mkdir -p /tmp/pp
	cp $(SEGMENT_FASTA) $(SEGMENT_COV_MATRIX) /tmp/pp
	docker run -it --rm \
		-v /var/run/docker.sock:/var/run/docker.sock \
		-v /tmp/pp:/work \
	        metabat/metabat:latest \
	        metabat2 -s 1500 -m 1500 --maxP 95 --minS 60 --maxEdges 200 --seed 1 -l --saveCls \
			-i /work/segments.fa \
			-a /work/segments.cov \
			-o /work/out
	cp /tmp/pp/out.MemberMatrix.txt $(SEGMENT_BIN_BASE)

post_metaBAT:
	Rscript utils/bin_summary.r \
		$(SEGMENT_BIN_BASE) \
	        $(REFINE_SEGMENTS) \
	        segment \
		$(BASE_BINS)

refine_bins:
	bin/polypanner refine_bins \
	        -threads 80 \
	        -ifn_libs $(POP_TABLE) \
	        -ifn_segments $(REFINE_SEGMENTS) \
	        -ifn_segments_binned $(SEGMENT_BIN_BASE) \
	        -p_value 0.01 \
	        -margin 10 \
	        -max_side_length 2000 \
	        -pseudo_count 0.1 \
	        -omit_bin 0 \
	        -only_bin NA \
	        -ofn $(SEG_BIN_FINAL)

