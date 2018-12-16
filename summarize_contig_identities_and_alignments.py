from handlers.BamHandler import BamHandler
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from matplotlib import pyplot, patches
import argparse
import os

'''
Generate stats/plots on contig identity and alignment given a BAM of contigs VS true reference
'''


READ_ID = 0
REVERSAL_STATUS = 1
REF_ALIGNMENT_START = 2
ALIGNMENT_LENGTH = 3
READ_LENGTH = 4
CONTIG_LENGTH = 5
N_INITIAL_CLIPPED_BASES = 6
N_TOTAL_MISMATCHES = 7
N_TOTAL_DELETES = 8
N_TOTAL_INSERTS = 9
IDENTITY = 10


def plot_contig_blocks(axes, y, ref_alignment_start, n_initial_clipped_bases, alignment_length, scale, color, contig_length):
    """
    Given relevant data about a contig, generate a plot showing how it aligns to the chromosome
    :param axes:
    :param ref_alignment_start:
    :param n_initial_clipped_bases:
    :param alignment_length:
    :param scale:
    :param color:
    :param contig_length:
    :return:
    """
    # ---- plot alignment block ----
    x = ref_alignment_start / scale
    width = alignment_length / scale
    height = 0.2

    rect = patches.Rectangle((x, y), width, height, color=color)
    axes.add_patch(rect)

    # ---- plot left side clipped block ----

    x = ref_alignment_start / scale - n_initial_clipped_bases / scale
    width = n_initial_clipped_bases / scale
    height = 0.2

    rect = patches.Rectangle((x, y), width, height, color=color, alpha=0.2)
    axes.add_patch(rect)

    # ---- plot right side clipped block ----

    x = ref_alignment_start / scale + alignment_length / scale
    width = contig_length / scale - alignment_length / scale - n_initial_clipped_bases / scale
    height = 0.2

    rect = patches.Rectangle((x, y), width, height, color=color, alpha=0.2)
    axes.add_patch(rect)


def plot_contigs(output_dir, read_data, chromosome_name, chromosome_length, total_identity, bam_path, y_min=None, y_max=None, save=True, show=False):
    """
    Create a figure showing how contigs align to the reference chromosome and what their identities are
    :param output_dir: where to save figures
    :param read_data: summary object constructed by "parse_reads" function containing all the alignment data
    :param chromosome_name:
    :param chromosome_length:
    :param total_identity: weighted average of all contig's identities
    :param bam_path: path to BAM containing alignment of contigs to true reference
    :param y_min: force a cutoff for y (useful for comparing multiple BAMs with varying contigs)
    :param y_max: force a cutoff for y (useful for comparing multiple BAMs with varying contigs)
    :param show: open a window to display the figure
    :param save: save the figure using the same name prefix as the BAM
    :return:
    """
    figure = pyplot.figure()
    axes = pyplot.axes()

    read_data = sorted(read_data, key=lambda x: x[READ_ID])

    left_column_offset = -6
    column_width = 3

    scale = 1000000

    if chromosome_length < scale:
        scale = chromosome_length

    # ---- plot reference chromosome and REF text legend ----

    width = chromosome_length/scale
    height = 0.2

    rect = patches.Rectangle((0, 0), width, height, color="gray")
    axes.add_patch(rect)

    x = -(chromosome_length / scale) + left_column_offset
    pyplot.text(x, 0, "REF", weight="bold")

    n_forward = 0
    n_reverse = 0
    for r,read in enumerate(read_data):
        read_id = read[READ_ID]
        reversal_status = read[REVERSAL_STATUS]
        ref_alignment_start = read[REF_ALIGNMENT_START]
        alignment_length = read[ALIGNMENT_LENGTH]
        contig_length = read[CONTIG_LENGTH]
        n_initial_clipped_bases = read[N_INITIAL_CLIPPED_BASES]
        identity = read[IDENTITY]

        if reversal_status:
            n_reverse += 1
            color = [77.6, 39.6, 7.1]
            color = [c/100 for c in color]
            y = - n_reverse/2
        else:
            n_forward += 1
            color = [4.3, 45.1, 47.1]
            color = [c/100 for c in color]
            y = n_forward/2

        plot_contig_blocks(axes=axes,
                           y=y,
                           ref_alignment_start=ref_alignment_start,
                           n_initial_clipped_bases=n_initial_clipped_bases,
                           alignment_length=alignment_length,
                           scale=scale,
                           color=color,
                           contig_length=contig_length)

        # ---- plot divider lines and text data ----

        axes.axhline(y + 0.33, linestyle="--", linewidth=0.6)

        x = -(chromosome_length/scale) + left_column_offset
        pyplot.text(x, y, read_id)

        x = -(chromosome_length/scale) + left_column_offset + column_width
        pyplot.text(x, y, str(round(identity*100, 3)))

    # ---- plot leftover divider lines ----

    axes.axhline(-(n_reverse+1)/2 + 0.33, linestyle="--", linewidth=0.6)
    axes.axhline(0 + 0.33, linestyle="--", linewidth=0.6)

    # ---- find limits for axes ----

    x_min = -(chromosome_length/scale) + left_column_offset
    x_max = 2*(chromosome_length/scale)

    axes.set_xlim(x_min, x_max)

    if y_min is None:
        y_min = -(n_reverse / 2 + 1)

    if y_max is None:
        y_max = n_forward / 2 + 0.5

    axes.set_ylim(y_min, y_max)

    # ---- add TOTAL identity text and ----

    x = -(chromosome_length / scale) + left_column_offset
    y = y_min
    pyplot.text(x, y, "TOTAL:", weight="bold")

    x = -(chromosome_length / scale) + left_column_offset + column_width
    pyplot.text(x, y, str(round(total_identity*100, 3)))

    # ---- chromosome title ----

    pyplot.title("Chromosome: "+chromosome_name)

    if save:
        filename = os.path.basename(bam_path)
        filename_prefix = ".".join(filename.split(".")[:-1])
        output_filename = filename_prefix + ".png"
        output_file_path = os.path.join(output_dir, output_filename)
        pyplot.axis("off")

        figure.savefig(output_file_path)

    if show:
        pyplot.show()
        pyplot.close()


def parse_match(alignment_position, length, read_sequence, ref_sequence):
    """
    Process a cigar operation that is a match
    :param alignment_position: Position where this match happened
    :param read_sequence: Read sequence
    :param ref_sequence: Reference sequence
    :param length: Length of the operation
    :return:

    This method updates the candidates dictionary.
    """
    start = alignment_position
    stop = start + length

    n_mismatches = 0

    for i in range(start, stop):
        allele = read_sequence[i-alignment_position]
        ref = ref_sequence[i-alignment_position]

        if allele != ref:
            # mismatch
            n_mismatches += 1

    return n_mismatches


def parse_cigar_tuple(cigar_code, length, alignment_position, read_sequence, ref_sequence):
    """
    Parse through a cigar operation to find possible candidate variant positions in the read
    :param cigar_code: Cigar operation code
    :param length: Length of the operation
    :param alignment_position: Alignment position corresponding to the reference
    :param read_sequence: Read sequence
    :return:
    cigar key map based on operation.
    details: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    0: "MATCH",
    1: "INSERT",
    2: "DELETE",
    3: "REFSKIP",
    4: "SOFTCLIP",
    5: "HARDCLIP",
    6: "PAD"
    """
    ref_index_increment = length
    read_index_increment = length

    n_mismatches = 0
    n_deletes = 0
    n_inserts = 0

    # deal different kinds of operations
    if cigar_code == 0:
        # match
        n_mismatches = parse_match(alignment_position=alignment_position,
                                   length=length,
                                   read_sequence=read_sequence,
                                   ref_sequence=ref_sequence)

    elif cigar_code == 1:
        # insert
        ref_index_increment = 0
        n_inserts = length

    elif cigar_code == 2 or cigar_code == 3:
        # delete or ref_skip
        read_index_increment = 0
        n_deletes = length

    elif cigar_code == 4:
        # soft clip
        ref_index_increment = 0

    elif cigar_code == 5:
        # hard clip
        ref_index_increment = 0
        read_index_increment = 0

    elif cigar_code == 6:
        # pad
        ref_index_increment = 0
        read_index_increment = 0

    else:
        raise ("INVALID CIGAR CODE: %s" % cigar_code)

    return ref_index_increment, read_index_increment, n_mismatches, n_deletes, n_inserts


def get_read_stop_position(read):
    """
    Returns the stop position of the reference to where the read stops aligning
    :param read: The read
    :return: stop position of the reference where the read last aligned
    """
    ref_alignment_stop = read.reference_end

    # only find the position if the reference end is fetched as none from pysam API
    if ref_alignment_stop is None:
        positions = read.get_reference_positions()

        # find last entry that isn't None
        i = len(positions) - 1
        ref_alignment_stop = positions[-1]
        while i > 0 and ref_alignment_stop is None:
            i -= 1
            ref_alignment_stop = positions[i]

    return ref_alignment_stop


def parse_reads(reads, chromosome_name, fasta_handler):
    """
    Given a set of pysam read objects, generate data for matches/mismatches/inserts/deletes and contig size/position for
    each read
    :param reads: pysam aligned segment objects
    :param chromosome_name:
    :param fasta_handler: fasta_handler object that can retrieve substrings from the reference sequence
    :return:
    """
    read_data = list()

    for read in reads:
        ref_alignment_start = read.reference_start
        ref_alignment_stop = get_read_stop_position(read)
        ref_length = ref_alignment_stop - ref_alignment_start

        reversal_status = read.is_reverse

        ref_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                  start=ref_alignment_start,
                                                  stop=ref_alignment_stop + 10)

        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence
        read_length = len(read_sequence)
        contig_length = read.infer_read_length()

        read_id = read.query_name
        # read_quality = read.query_qualities

        # read_index: index of read sequence
        # ref_index: index of reference sequence
        read_index = 0
        ref_index = 0
        found_valid_cigar = False

        n_total_mismatches = 0
        n_total_deletes = 0
        n_total_inserts = 0

        n_initial_clipped_bases = 0

        for c, cigar in enumerate(cigar_tuples):
            cigar_code = cigar[0]
            length = cigar[1]

            # get the sequence segments that are effected by this operation
            read_sequence_segment = read_sequence[read_index:read_index + length]
            ref_sequence_segment = ref_sequence[ref_index:ref_index+length]

            # skip parsing the first segment if it is not a match
            if cigar_code != 0 and found_valid_cigar is False:
                # only increment the read index if the non-match cigar code is INS or SOFTCLIP
                if cigar_code == 1 or cigar_code == 4:
                    read_index += length
                if cigar_code == 5 or cigar_code == 4:
                    n_initial_clipped_bases = length
                continue

            found_valid_cigar = True

            # send the cigar tuple to get attributes we got by this operation
            ref_index_increment, read_index_increment, n_mismatches, n_deletes, n_inserts = \
                parse_cigar_tuple(cigar_code=cigar_code,
                                  length=length,
                                  alignment_position=ref_alignment_start + ref_index,
                                  read_sequence=read_sequence_segment,
                                  ref_sequence=ref_sequence_segment)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment
            n_total_mismatches += n_mismatches
            n_total_deletes += n_deletes
            n_total_inserts += n_inserts

        # total_non_matches = n_total_mismatches + n_total_deletes + n_total_inserts
        total_non_matches = 2*n_total_mismatches + n_total_deletes + n_total_inserts

        # identity = (ref_length - total_non_matches) / ref_length
        identity = (ref_length + read_length - total_non_matches) / (ref_length + read_length)

        data = [read_id,
                reversal_status,
                ref_alignment_start,
                ref_length,
                read_length,
                contig_length,
                n_initial_clipped_bases,
                n_total_mismatches,
                n_total_deletes,
                n_total_inserts,
                identity]

        read_data.append(data)

    return read_data


def process_bam(bam_path, reference_path, output_dir=None):
    """
    Find useful summary data from a bam that can be represented as a table of identities, and a plot of alignments
    :param bam_path: path to a bam containing contigs aligned to a true reference
    :param reference_path: the true reference that contigs were aligned to
    :param output_dir: where to save plots
    :return:
    """
    print("\n" + bam_path + "\n")

    if output_dir is None:
        output_dir = "plots/"

    FileManager.ensure_directory_exists(output_dir)

    bam_handler = BamHandler(bam_file_path=bam_path)
    fasta_handler = FastaHandler(reference_path)

    chromosome_names = fasta_handler.get_contig_names()

    for chromosome_name in chromosome_names:
        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        start = 0
        stop = chromosome_length

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=start, stop=stop)

        read_data = parse_reads(reads=reads, fasta_handler=fasta_handler, chromosome_name=chromosome_name)

        print("chromosome_name:\t", chromosome_name)
        print("chromosome_length:\t", chromosome_length)

        for data in read_data:
            read_id, reversal_status, ref_alignment_start, alignment_length, read_length, contig_length, n_initial_clipped_bases, n_total_mismatches, n_total_deletes, n_total_inserts, identity = data
            print()
            print(read_id)
            print("reversed:\t", reversal_status)
            print("alignment_start:\t", ref_alignment_start)
            print("alignment_length:\t", alignment_length)
            print("read_length:\t\t", read_length)
            print("n_initial_clipped_bases:", n_initial_clipped_bases)
            print("n_total_mismatches:\t", n_total_mismatches)
            print("n_total_deletes:\t", n_total_deletes)
            print("n_total_inserts:\t", n_total_inserts)
            print("identity:\t", identity)

        total_weighted_identity = sum([x[ALIGNMENT_LENGTH] * x[IDENTITY] for x in read_data])
        total_alignment_bases = sum([x[ALIGNMENT_LENGTH] for x in read_data])
        total_identity = total_weighted_identity/total_alignment_bases

        print("\nTOTAL IDENTITY:\t", total_identity)

        plot_contigs(output_dir=output_dir,
                     read_data=read_data,
                     chromosome_name=chromosome_name,
                     chromosome_length=chromosome_length,
                     total_identity=total_identity,
                     bam_path=bam_path,
                     y_min=-1,
                     y_max=4,
                     show=False)


def main(bam_path, reference_path, output_dir):

    process_bam(bam_path=bam_path, reference_path=reference_path, output_dir=output_dir)


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file path of contigs aligned to true reference"
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="FASTA file path of true reference to be compared against"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )

    args = parser.parse_args()

    main(bam_path=args.bam, reference_path=args.ref, output_dir=args.output_dir)