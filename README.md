# nanopore_assembly_and_polishing_assessment
Automation of pipelines that depend on preexisting assembly, polishing, and alignment tools. Performance evaluation and visualization of assembly and polishing

## Requirements

- python3

- pysam (using htslib >=1.7):
  - ubuntu 16.04 - NOT working (htslib version not supported)
  - ubuntu 18.04 - WORKING
  - macOS - UNTESTED
  
- matplotlib

- tkinter


## Examples

### Alignment Visualization

Running contig identity and mapping assessment:
```
align_and_summarize_contigs.py --sequences SEQUENCES --ref REF [--output_dir OUTPUT_DIR]

Required arguments:
  --sequences SEQUENCES   FASTA file path of contig sequences
  --ref REF               FASTA file path of true reference to be compared
                          against

Optional Arguments:
  --output_dir OUTPUT_DIR
                        desired output directory path (will be created during
                        run time if doesn't exist)
```

Using a polished E. Coli assembly aligned to the k12 reference (BAM), this script yields the following plot:

![example output plot](https://github.com/rlorigro/nanopore_assembly_and_polishing_assessment/raw/master/assembled_wtdbg2_r94_ec_rad2_30x-30kb_VS_refEcoli.sorted.png)

It is apparent that the aligner has chosen to break up the contig into multiple supplementary alignments because the chromosome is circular, and there appear to be poorly assembled regions along it. The following (verbose data) is printed to stdout (and exported to a tabular csv format):

```
chromosome_name:         gi
chromosome_length:       4641652

ctg1
reversed:        False
alignment_start:         0
alignment_length:        272186
read_length:             270683
n_initial_clipped_bases: 4345687
n_total_mismatches:      283
n_total_deletes:         1728
n_total_inserts:         225
identity:        0.9953598381930079

... (more output for each alignment omitted for brevity) ...

TOTAL IDENTITY:  0.9954066201845397
```

### NG(A)x Generation and Plotting

![example output plot](https://github.com/rlorigro/nanopore_assembly_and_polishing_assessment/blob/master/shasta_2019_6_28_10_39_18_878109.png)

### Coverage Estimation

![example output plot](https://github.com/rlorigro/nanopore_assembly_and_polishing_assessment/blob/master/squashed_coverage_733_maternal.png)
