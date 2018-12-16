# nanopore_assembly_and_polishing_assessment
Automation of pipelines that depend on preexisting assembly, polishing, and alignment tools. Performance evaluation and visualization of assembly and polishing


## Examples

Running contig identity and mapping assessment:
```
summarize_contig_identities_and_alignments.py --bam BAM --ref REF [--output_dir OUTPUT_DIR]

Required arguments:
  --bam BAM             BAM file path of contigs aligned to true reference
  --ref REF             FASTA file path of true reference to be compared
                        against

Optional Arguments:
  --output_dir OUTPUT_DIR
                        desired output directory path (will be created during
                        run time if doesn't exist)
```

Using a polished E. Coli assembly aligned to the k12 reference (BAM), this script yields the following plot:

![example output plot](https://github.com/rlorigro/nanopore_assembly_and_polishing_assessment/raw/master/assembled_wtdbg2_r94_ec_rad2_30x-30kb_VS_refEcoli.sorted.png)

It is apparent that the aligner has chosen to break up the contig into multiple supplementary alignments because the chromosome is circular, and there appear to be poorly assembled regions along it. The following (verbose data) is printed to stdout:

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
