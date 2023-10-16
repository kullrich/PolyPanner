# PolyPanner: Detecting polymorphic sites in complex microbial communities

PolyPanner is a suite of algorithms that fasicilate the detection of polymorphic variants in metagenomic assembled genomes. It leverages dense temporal sampling to improve assembly quality and identify high confidence polymorphic variants, called *dynamic variants*, which are polymorphic variants with frequency variation over time too large to be explained by read sampling stochasticity. The remaining variants, which we call spurious variants, reflect either homology between distinct populations (*ortholog variants*) or homology within a population due to duplicated genes (*paralog variants*), and are considered noise for our purpose of tracking competing alleles within a single bacterial population.

As input, PolyPanner receives a set of shotgun DNA libraries that have been aligned to a metagenomic co-assembly. As a form of data representation, alignments are transformed to single-nucleotide coverage vectors that represent library-specific read counts of perfect and mismatch alignments at each base pair in the co-assembly. 
Tasks performed by PolyPanner are:
- The modelling and removal of sequencing errors.
- The refinement of assembly contigs into segments, through the introduction of breakpoints where there are abrupt transitions in coverage. These transitions likely reflect chimeric assembly breakpoints where the two sides of the breakpoint represent different populations. 
- The refinement of genome bins, by trimming-out segments that differ in their coverage profiles.
- The identification of *dynamic variants*.

PolyPanner was developed by Eitan Yaffe (eitan.yaffe@gmail.com).

## Installation 

You use pre-compiled binaries or compile from code.

### Use pre-compiled binaries

Binaries for macOS (Ventura 13.3.1) or Ubuntu (20.04.1) are part of the v1.0 release. You may need to run ```chmod +x $BINARY``` if they are not executable after downloading them.

MacOS: https://github.com/eitanyaffe/PolyPanner/releases/download/v1.0.0/polypanner-v1.0.0-macOS-i386.macos

Linux: https://github.com/eitanyaffe/PolyPanner/releases/download/v1.0.0/polypanner-v1.0.0-linux-x86-64.linux


### Compile from code

You will need to install gcc, gsl and boost. Compilation should take under an hour.  

On MacOS run:
```
brew install g++ gsl boost
```

On Ubuntu run:
```
apt-get update
apt-get install git build-essential libgsl0-dev libboost-all-dev
```

Clone the repository from github, and compile PolyPanner:
```
git clone https://github.com/eitanyaffe/PolyPanner.git
cd PolyPanner
make
```

### Validate installation

To validate the installation you can run the following unit test which should take under 20 seconds to complete. Docker is required to run metaBAT2 during the unit test. The unit test runs on a small mock dataset composed of 4 DNA libraries, each with ~20k paired reads (2x150nt). 

```
make test
```

A recommended way to start using PolyPanner is by customizing the commands of the unit test. The bash commands can be printed by peforming a dry-run:

```
make test -n
```

## Workflow overview


The standard PolyPanner workflow is detailed below. The syntax of PolyPanner commands and utility scripts, including user-defined parameters, is documented [here](docs/syntax.md).

### 1. Input

The workflow input is composed of:
- A co-assembly with a contig fasta file and a contig table (fields: contig and length).
- A set of mapped paired DNA libraries in the SAM format (mapped separately for R1 and R2).

### 2. POP construction

For each library, SAM files are converted to an internal PolyPanner tabular format using the ```utils/parse_bwa_sam.pl``` script. Read sides are paired using the ```utils/pair_reads.pl``` script. Paired reads are transformed into POP files, which represent mapped reads in an efficient manner that allows quick queries on each library.

### 3. Removal of sequencing errors

Libraries are merged into a single library using the ``polypanner merge`` command. Sequencing errors are identified and removed using the ``polypanner filter`` command. Libraries are restricted to true segregating sites using the ``polypanner restrict`` command. A tab-delimited table with paths to the restricted libraries is generated using the ``utils/make_pop_table.r`` script.

### 4. Assembly refinement

Contigs are refined into segments using the ``polypanner refine`` command.

### 5. Genome binning 

Genome binning is performed by an external binning tool. To use metaBAT2, a segment coverage matrix (input for metaBAT2) is generated using the ``polypanner cov_matrix``. The output of the binning step is required to be a tab-delimited segment-bin table  with 2 columns: segment and bin. The workflow can support other binning methods, just make sure to generate the segment-bin table as specified. 

### 6. Genome refinement 

The segment-bin table is processed using the ``utils/bin_summary.r`` script. Bins are refined using the ``polypanner refine_bins`` command.

### 7. Site annotation

Sites are annotated using the ``polypanner sites`` command.

### 8. Output

Genome bin trajectories are generated using the ``polypanner bin_trajectory`` command. Site trajectories are generated using the ``polypanner site_trajectory``.

Other output files include the refined segment-bin table and the annotated site table. 

## Output formats

### Site table

Fields:
- vid: Variant identifer.
- contig/coord: The contig and coordinate within contig.
- variant: String identifier of the variant.
- segment: The segment.
- bin: The genome bin.
- var_count: Number of reads supporting the variant.
- total_count: Number of reads supporting the position.
- n_samples: Number of samples that have reads supporting the variant.
- is_internal: Is variant close to segment edge.
- variant_p: p-value of the test comparing the variant and the local coverages.
- regional_p: p-value of the test comparing the local and regional coverages.
- comp_p: p-value of the test comparing the complementary and regional coverages.

Dynamic sites are defined as sites for which is_internal equals T, variant_p<X, regional_p>X and comp_p<X, where X=0.05.
