# polypanner syntax

All polypanner commands are called through the ``polypanner`` program. For example, to construct a POP file call:

```
bin/polypanner construct \
        -ifn_paired output/test/t1/paired \
        -ifn_R1 output/test/t1/non_paired_R1 \
        -ifn_R2 output/test/t1/non_paired_R2 \
        -contig_table input/test/contig_table \
        -discard_clipped discard_both \
        -min_score 30 \
        -min_length 50 \
        -max_edit 20 \
        -trim 20 \
        -unique_only T \
        -ofn output/test/t1/lib.pop
```

## polypanner construct

Construct POP from mapped reads.

```
usage: polypanner construct [options]
 -contig_table <fn>: contig table (mandatory)
 -ifn_paired <fn>: input file table of mapped paired reads (mandatory)
 -ifn_R1 <fn>: input file table of mapped non-paired reads R1 (mandatory)
 -ifn_R2 <fn>: input file table of mapped non-paired reads R2 (mandatory)
 -discard_clipped <string>: handle clipped reads
 -trim <int>: keep trim from read sides
 -min_length <int>: minimal match length (nts)
 -min_score <int>: minimal score
 -unique_only <T|F>: limit to unique hits
 -max_edit <int>: maximal edit distance
 -ofn <fn>: output variation table (mandatory)
 -query_contig <string>: query contig (debug)
 -query_coord <int>: query coordinate (debug)
 -query_ofn <fn>: output query filename (debug)
discard_clipped options: 
- discard_any: discard reads clipped on either side
- discard_both: discard reads clipped on both side
- keep: discard reads clipped on both side
```

## polypanner merge

Merge multiple CAVs.

```
usage: polypanner merge <ofn> [ifn1 ifn2 ...]
```

## polypanner filter

Infer sequencing error rates and true segregating sites.

```
usage: polypanner filter [options]
 -ifn <fn>: CAV filename (mandatory)
 -threads <int>: Number of threads used
 -alpha <double>: False discovery rate
 -min_variant_count <int>: Minimal count of considered variant
 -min_contig_length <int>: Minimal contig length
 -min_error_pct <double>: Minimal error rate (%)
 -max_error_pct <double>: Maximal error rate (%)
 -seed_error_pct <double>: Seed error rate (%) used at start of optimization
 -ofn_sites <fn>: Output table with segregating sites (mandatory)
 -ofn_error_params <fn>: Output table with segregating sites (mandatory)
```

## polypanner restrict

Restrict CAV to a set of contigs.

```
usage: polypanner restrict [options]
 -ifn_cav <fn>: CAV filename (mandatory)
 -ifn_sites <fn>: Table with variants (mandatory)
 -ofn <fn>: Output CAV file (mandatory)
```

## polypanner refine

Breakdown contigs into segments, on dangle coords with a transition in coverage distribution.

```
usage: polypanner refine [options]
 -ifn <fn>: Table with multiple CAV files (mandatory)
 -threads <int>: Number of threads (mandatory)
 -fdr <double>: False discovery rate
 -ofn_segments <fn>: output segment file (mandatory)
 -ofn_contigs <fn>: output contig file (mandatory)
 -ofn_breakpoints <fn>: output breakpoint file (mandatory)
 -p_value <double>: Min Chi-Square P-value
 -margin <int>: Margin away from breakpoints
 -cand_step <int>: window step size for non-dangle candidates
 -min_contig_length <int>: Min contig length
 -min_length <int>: Min breakpoint side length
 -max_length <int>: Max breakpoint side length
 -pseudo_count <double>: Add pseudo-count
 -max_lib_count <int>: Stop after X libs for debugging (if set to 0 then use all)
 -only_dangles <T|F>: Check only coords which end a trimmed read
 -contig <string>: Only single contig
```

## polypanner cov_matrix

Generate coverage matrix for segments.

```
usage: polypanner cov_matrix [options]
 -ifn_cavs <fn>: Table with multiple CAV files (mandatory)
 -ifn_segs <fn>: Segment table (mandatory)
 -ifn_fasta <fn>: Contig fasta (mandatory)
 -actual_nts <T|F>: Actual nts or random nts used to avoid TNF usage) (mandatory)
 -ofn_fasta <fn>: Output fasta file (mandatory)
 -ofn_mat <fn>: Output matrix with segments, bp-coverage and bp-variance (mandatory)
```

## polypanner refine_bins

Breakdown contigs into segments, on dangle coords with a transition in coverage distribution.

```
usage: polypanner refine_bins [options]
 -ifn_libs <fn>: Table with multiple CAV files (mandatory)
 -ifn_segments <fn>: input segment file (mandatory)
 -ifn_segments_binned <fn>: input table associating segments (1st column) to bins (2nd column) (mandatory)
 -ofn <fn>: output binning table (mandatory)
 -threads <int>: Number of threads used
 -omit_bin <string>: Discard input bin with this value
 -only_bin <string>: Limit to specific bin
 -p_value <double>: Chi-Square P-value clustering threshold
 -margin <int>: Margin away from breakpoints
 -pseudo_count <double>: Add pseudo-count
 -max_side_length <int>: Maximal distance from fragment side used for coverage
```


## polypanner site_trajectory

Extract site trajectories.

```
usage: polypanner site_trajectory [options]
 -ifn_libs <fn>: Table with multiple CAV files (mandatory)
 -ifn_segments <fn>: input segment file, with bin field (mandatory)
 -ifn_sites <fn>: input table with sites (mandatory)
 -ofn_counts <fn>: output matrix with var counts (mandatory)
 -ofn_totals <fn>: output matrix with total counts (mandatory)
```

## polypanner bin_trajectory

Extract average bin trajectories.

```
usage: polypanner bin_trajectory [options]
 -ifn_libs <fn>: Table with multiple CAV files (mandatory)
 -ifn_segments <fn>: input table with binned segments (mandatory)
 -ofn <fn>: output matrix with bin reads counts per sample (mandatory)
```

## polypanner sites

```
usage: polypanner sites [options]
 -ifn_libs <fn>: table with multiple CAV files (mandatory)
 -ifn_sites <fn>: input site table (mandatory)
 -ifn_contigs <fn>: contig table (mandatory)
 -ifn_segments <fn>: segment table (mandatory)
 -ifn_segment_bins <fn>: segment-bin table (mandatory)
 -regional_window <int>: window used for regional coverage (bp)
 -margin <int>: flag sites close to contig end (bp) (mandatory)
 -pseudo_count <double>: pseudo-count for coverage counts
 -ofn <fn>: output site file (mandatory)
```

## polypanner fasta

Create fasta files for bins.

```
usage: polypanner fasta [options]
 -ifn_fasta <fn>: input contig fasta file (mandatory)
 -ifn_bins <fn>: input bin table (mandatory)
 -ifn_segments <fn>: input segment file, with bin field (mandatory)
 -min_length <int>: minimal bin length
 -odir <fn>: output path (mandatory)
```

# Utility scripts

## Rscript utils/bin_summary.r 

Creates summary table of genome bins.

```
usage: utils/bin_summary.r <contig-bin table> <contigs> <contig field> <ofn>
```

## Rscript utils/make_pop_table.r 

Creates table with paths to POP files.

```
usage: utils/make_pop_table.r <base dir> <filename> <ofn> <sample1> <sample2> ...
```

## perl utils/pair_reads.pl 

Pair read sides into single table.

```
usage: utils/pair_reads.pl <read1 ifn> <read2 ifn> <select single mapped read with maximal field value> <ofn> <ofn stats>
```

## perl utils/parse_bwa_sam.pl 

Parse bwa SAM output into the input format supported by PolyPanner.

```
usage: utils/parse_bwa_sam.pl <ifn> <ofn> <ofn stats>
```
