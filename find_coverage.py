from matplotlib import pyplot
import argparse
import pysam
import math
import sys

pyplot.switch_backend('agg')


def get_coverage(bam_file_path, chromosome_name, start, stop):
    sam_file = pysam.AlignmentFile(bam_file_path, "rb")

    pileup_columns = sam_file.pileup(chromosome_name, start, stop)
    coverages = list()

    for pileup_column in pileup_columns:
        position = pileup_column.pos
        coverage = pileup_column.nsegments

        if start < position < stop:
            coverages.append(coverage)

    sam_file.close()

    return coverages


def test():
    # BAM of reads aligned to some sequence
    bam_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_sequence_2019_3_27_13_13_54_735626/sequence_subset_train_60x_10kb_rle_VS_refEcoli_rle.sorted.bam"

    # Fasta containing the sequences that are the reference in the BAM
    fasta_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_sequence_2019_3_27_11_15_11_736904/refEcoli_rle.fasta"

    fasta_handler = pysam.FastaFile(fasta_path)
    chromosome_names = fasta_handler.references
    print(chromosome_names)

    n_samples = 100

    coverages = list()
    for name in chromosome_names:
        length = fasta_handler.get_reference_length(name)

        start = 0
        stop = length
        step_size = int(round(length/n_samples))

        steps = range(start, stop, step_size)

        for c,coord in enumerate(steps):
            coord = int(math.floor(coord))
            coverage = get_coverage(bam_file_path=bam_path, chromosome_name=name, start=coord, stop=coord+2)[0]
            coverages.append(coverage)
            sys.stderr.write("\r %.2f%%" % (c+1/n_samples*100))

        sys.stderr.write("\n")

    pyplot.plot(coverages)
    pyplot.show()
    pyplot.close()


def main(bam_path, chromosome_name, start, stop, n_samples):
    """
    :param bam_path: BAM of reads aligned to some sequence
    :param fasta_path: Fasta containing the sequences that are the reference in the BAM
    :param chromosome_name:
    :param start:
    :param stop:
    :param n_samples: How many times to check coverage within region (regularly spaced)
    :return:
    """
    print("Testing chromosome '%s' from %d to %d" % (chromosome_name, start, stop))

    coverages = list()

    step_size = int(round((stop-start)/n_samples))
    steps = range(start, stop, step_size)

    print(steps)

    max_coverage = 0
    for c,coord in enumerate(steps):
        coord = int(math.floor(coord))
        coverage = get_coverage(bam_file_path=bam_path, chromosome_name=chromosome_name, start=coord, stop=coord+2)[0]
        coverages.append(coverage)
        sys.stderr.write("\r %.2f%%" % (c+1/n_samples*100))

        if coverage > max_coverage:
            max_coverage = coverage

    sys.stderr.write("\n")

    axes = pyplot.axes()
    axes.plot(coverages)
    axes.set_ylim([0, round(max_coverage*1.1)])
    axes.set_title("Coverage on chromosome '%s' from %d to %d" % (chromosome_name, start, stop))
    axes.set_ylabel("Coverage (# reads)")
    axes.set_xlabel("Interval")

    pyplot.savefig("coverage.png")

    # pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file path of contigs aligned to true reference"
    )
    parser.add_argument(
        "--contig", "-c",
        type=str,
        required=True,
        help="Name of reference chromosome/contig to query"
    )
    parser.add_argument(
        "--start", "-a",
        type=int,
        required=False,
        help="Region start coordinate"
    )
    parser.add_argument(
        "--stop", "-b",
        type=int,
        required=False,
        help="Region stop coordinate"
    )
    parser.add_argument(
        "--samples", "-n",
        type=int,
        required=False,
        default=100,
        help="How many times to sample the region"
    )

    args = parser.parse_args()

    main(bam_path=args.bam, chromosome_name=args.contig, start=args.start, stop=args.stop, n_samples=args.samples)
