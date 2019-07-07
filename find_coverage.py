import matplotlib
from matplotlib import pyplot
import argparse
import platform
import pysam
import math
import sys
import os


if os.environ.get("DISPLAY", "") == "":
   print("no display found. Using non-interactive Agg backend")
   matplotlib.use("Agg")

if platform.system() == "Darwin":
   matplotlib.use("macosx")


def get_coverage_range(bam_file_path, chromosome_name, start, stop):
    sam_file = pysam.AlignmentFile(bam_file_path, "rb")

    pileup_columns = sam_file.pileup(chromosome_name, start, stop)
    coverages = list()
    positions = list()

    max_coverage = 0
    for pileup_column in pileup_columns:
        position = pileup_column.pos
        position_coverage = pileup_column.nsegments

        if start < position < stop:
            coverages.append(position_coverage)
            positions.append(position)

            if position_coverage > max_coverage:
                max_coverage = position_coverage

    sam_file.close()

    return positions, coverages, max_coverage


def get_coverage(bam_file_path, chromosome_name, coordinate):
    sam_file = pysam.AlignmentFile(bam_file_path, "rb")

    pileup_columns = sam_file.pileup(chromosome_name, coordinate, coordinate+2)
    coverage = 0

    for pileup_column in pileup_columns:
        position = pileup_column.pos
        position_coverage = pileup_column.nsegments

        if position == coordinate:
            coverage = position_coverage

    sam_file.close()

    return coverage


def get_coverage_subsample(bam_path, chromosome_name, start, stop, n_samples):
    coverages = list()
    positions = list()

    step_size = int(round((stop-start)/n_samples))
    steps = range(start, stop, step_size)

    print(steps)

    max_coverage = 0
    for c,coord in enumerate(steps):
        coord = int(math.floor(coord))
        coverage = get_coverage(bam_file_path=bam_path, chromosome_name=chromosome_name, coordinate=coord)
        coverages.append(coverage)
        positions.append(coord)
        sys.stderr.write("\r %.2f%%" % ((c+1)/n_samples*100))

        if coverage > max_coverage:
            max_coverage = coverage

    return positions, coverages, max_coverage


def test():
    # BAM of reads aligned to some sequence
    bam_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_guppy_fasta_wg_60x_10kb/sequence_subset_train_60x_10kb_rle_VS_refEcoli_rle.sorted.bam"

    # Fasta containing the sequences that are the reference in the BAM
    fasta_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_guppy_fasta_wg_60x_10kb/refEcoli_rle.fasta"

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
            coverage = get_coverage(bam_file_path=bam_path, chromosome_name=name, coordinate=coord)
            coverages.append(coverage)
            sys.stderr.write("\r %.2f%%" % ((c+1)/n_samples*100))

        sys.stderr.write("\n")

    pyplot.plot(coverages)
    pyplot.show()
    pyplot.close()


def main(bam_path, chromosome_name, start, stop, n_samples, y_max=None):
    """
    :param bam_path: BAM of reads aligned to some sequence
    :param chromosome_name:
    :param start:
    :param stop:
    :param n_samples: How many times to check coverage within region (regularly spaced)
    :return:
    """
    print("Testing chromosome '%s' from %d to %d" % (chromosome_name, start, stop))

    if n_samples is None:
        # Sample every position
        positions, coverages, max_coverage = get_coverage_range(bam_file_path=bam_path,
                                                                chromosome_name=chromosome_name,
                                                                start=start,
                                                                stop=stop)

    else:
        # Sample only the specified number of times
        positions, coverages, max_coverage = get_coverage_subsample(bam_path=bam_path,
                                                                    chromosome_name=chromosome_name,
                                                                    start=start,
                                                                    stop=stop,
                                                                    n_samples=n_samples)

    print("\n", len(coverages))

    print(positions[:10])
    print(coverages[:10])

    sys.stderr.write("\n")

    axes = pyplot.axes()
    # axes.scatter(positions, coverages, s=0.5)
    axes.plot(positions, coverages)
    axes.set_title("Coverage on %s" % chromosome_name)
    axes.set_ylabel("Coverage (# reads)")
    axes.set_xlabel("Coordinate")

    if y_max is None:
        axes.set_ylim([0, round(max_coverage*1.1)])
    else:
        axes.set_ylim([0, y_max])

    pyplot.xticks(rotation=45)
    axes.ticklabel_format(useOffset=False, style='plain')
    pyplot.tight_layout()

    pyplot.savefig("coverage.png", dpi=300)
    pyplot.savefig("coverage.pdf", dpi=300)

    pyplot.show()
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
        help="How many times to sample the region"
    )
    parser.add_argument(
        "--y_max", "-y",
        type=int,
        required=False,
        help="Optionally fix the y limit for plotting (coverage axis)"
    )

    args = parser.parse_args()

    main(bam_path=args.bam, chromosome_name=args.contig, start=args.start, stop=args.stop, n_samples=args.samples, y_max=args.y_max)
