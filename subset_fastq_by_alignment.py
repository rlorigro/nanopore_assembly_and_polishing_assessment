from handlers.FastqHandler import FastqHandler
import argparse
import pysam


def main(fastq_path, bam_path, target_chromosomes):
    """
    Iterate a BAM file and collect read names that belong to target chromosomes, and then subset a fastq based on names
    :param fastq_path:
    :param bam_path:
    :param target_chromosomes:
    :return:
    """
    bam_handler = pysam.AlignmentFile(bam_path, "rb")

    read_names = set()

    for chromosome_name in target_chromosomes:
        reads = bam_handler.fetch(contig=chromosome_name)

        for read in reads:
            if read.mapping_quality >= 5 \
                    and not read.is_secondary \
                    and not read.is_unmapped \
                    and not read.is_qcfail \
                    and not read.is_supplementary:

                read_names.add(read.query_name)

    fastq_handler = FastqHandler(fastq_path)

    fastq_handler.extract_reads_by_id(id_set=read_names, output_path="read_subset.fastq", print_status=True)


def comma_separated_set(s):
    tokens = s.strip().split(",")
    s = set(tokens)
    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.register("type", "comma_separated_set", comma_separated_set)

    parser.add_argument("--bam",
                        required=True,
                        type=str,
                        help="Path to BAM file containing reads aligned to a reference")

    parser.add_argument("--fastq",
                        required=True,
                        type=str,
                        help="Path to FASTQ file containing sequences")

    parser.add_argument("--contigs",
                        required=True,
                        type=comma_separated_set,
                        help="A comma-separted list of contig names. Only reads aligned to these contigs will be retained.")

    args = parser.parse_args()

    main(fastq_path=args.fastq, bam_path=args.bam, target_chromosomes=args.contigs)
