######################################################################################################
# get input files
######################################################################################################

SDIR?=/makeshift-mnt/output/val/v1.16/datasets/pp_unit_mini/v6/runs/r1

get_local_assembly:
	mkdir -p $(INPUT_DIR)
	cp $(SDIR)/assembly/v1.02/c1_1/k77_200M/work/contig_table $(CONTIG_TABLE)
	cp $(SDIR)/assembly/v1.02/c1_1/k77_200M/work/contigs $(CONTIG_FASTA)

get_local_sam:
	mkdir -p $(INPUT_LIB_DIR)
	cp $(SDIR)/map/v1.05/c1_1/libs/$(SAMPLE)/chunks/1/R1.sam $(SAM_R1)
	cp $(SDIR)/map/v1.05/c1_1/libs/$(SAMPLE)/chunks/1/R2.sam $(SAM_R2)
	
get_local_sams:
	$(foreach S,$(SAMPLES),$(MAKE) get_local_sam SAMPLE=$S && ) true

get_input: get_local_assembly get_local_sams