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

Using a nanopolished E. Coli assembly aligned to the k12 reference (BAM), this script yields the following output:

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

ctg4
reversed:        False
alignment_start:         58747
alignment_length:        8389
read_length:             8136
n_initial_clipped_bases: 35
n_total_mismatches:      524
n_total_deletes:         590
n_total_inserts:         248
identity:        0.8858698940998487

... (more output omitted) ...

TOTAL IDENTITY:  0.9954066201845397
```
