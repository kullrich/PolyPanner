all: bin/polypanner

###############################################################################################
# compliation evironment
###############################################################################################

INSTALL_DIR=/makeshift-mnt/bin

HFILES=$(addprefix cpp/,\
Variation.h VariationSet.h cav.h Params.h util.h BinMatrix.h Filter.h Resource.h \
RefineLocal.h)

OBJ=$(addprefix obj/,\
cav_sites.o \
cav_fasta.o \
cav_stats.o \
cav_bin_trajectory.o \
cav_vcluster_trajectory.o \
cav_site_trajectory.o \
Resource.o \
ClusterVariants.o \
cav_variant_cluster.o \
cav_refine_bins.o \
BinMatrix.o \
Filter.o \
Variation.o VariationSet.o Params.o util.o \
Dissolve.o RefineLocal.o \
cav.o cav_construct.o cav_filter.o \
cav_merge.o \
cav_combine.o \
cav_refine_local.o cav_refine_global.o \
cav_dump_local_scores.o \
cav_read_query.o \
cav_dump.o \
cav_restrict.o \
cav_info.o \
cav_cov_mat.o)

CFLAGS=-Wall -Wno-write-strings -std=c++0x -fext-numeric-literals
LDFLAGS=-pthread -lgsl -lgslcblas
CC=g++

obj:
	@mkdir -p $@
bin:
	@mkdir -p $@

obj/%.o: cpp/%.cpp $(HFILES)
	$(CC) $(CFLAGS) -c -o $@ $<

bin/polypanner: obj bin $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDFLAGS)

install: bin/polypanner
	cp bin/polypanner $(INSTALL_DIR)

###############################################################################################
# unit test
###############################################################################################

# input:
# - paired fastq files
# - assembly fasta (optional, requires megaHIT)
# - mapped reads (optional, requires bwa)

# output:
# - sites
# - genomes

# assembly
# mapping
# pp construct
# merge, filter errors
# refine assembly
# metabat2 binning
# trimming bins
# sites

SAMPLES?=t1 t2 t3 t4

construct:
	echo constrcuting polypanner file for sample $(SAMPLE)
	echo polypanner construct

# construct polypanner files
construct_all:
	$(foreach S,$(SAMPLES),$(MAKE) construct SAMPLE=$S)

test:
	$(MAKE) construct_all
