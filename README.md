# PolyPanner

PolyPanner is a suite of algorithms that fascilate the detection of polymorphic variants in metagenomic assembled genomes. It leverages dense temporal sampling to improve assembly quality and identify high confidence polymorphic variants, called *dynamic variants*, which are polymorphic variants with a non-constant frequency over time, after taking into account read sampling stochasticity. The remaining variants, which we call spurious variants, reflect either homology between distinct populations (*ortholog variants*) or homology within a population due to duplicated genes (*paralog variants*), and are considered noise for our purpose of tracking competing alleles within a single bacterial population.

As input, PolyPanner receives a set of shotgun DNA libraries that were aligned to a metagenomic co-assembly. As a form of data representation, alignments are transformed to single-nucleotide coverage vectors that represent library-specific read counts of perfect and mismatch alignments at each base pair in the co-assembly. Tasks performed by PolyPanner are:
- The modelling and removal of sequencing errors in reads.
- The refinement of assembly contigs into segments, through the introduction of breakpoints where there are aprupt transitions in coverage. These transitions likely reflect chimeric assembly breakpoints where the two sides of the breakpoint represent different populations. 
- The refinement of genome bins, by trimming-out segments that differ in their coverage profiles.
- The identification of *dynamic variants*.

PolyPanner was developed by Eitan Yaffe (eitan.yaffe@gmail.com).

## Installation

You will need to install gcc, gsl and boost. 

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

## Unit test

To validate the installation you can run a unit test. Docker is required to run metaBAT2 during the unit test.

```
make test
```

The unit test runs on a small mock dataset generated with these parameters:
- 3 genomes, 10kb each, 10 mutations per genome
- 4 timepoints
- mean coverage of 100x per genome per timepoint

## Workflow overview

The following is a standard workflow. 

Syntanx for all commands is documented [here](docs/syntax.md).

### Input

### PP Construction

### Removal of sequencing errors

### Assembly refinement

### Genome refinement 

### Site annotation

### Output

