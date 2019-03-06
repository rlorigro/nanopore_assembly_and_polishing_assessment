from handlers.BamHandler import BamHandler
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from multiprocessing import cpu_count, Pool
from matplotlib import pyplot, patches
from collections import defaultdict
import argparse
import csv
import os
pyplot.switch_backend('agg')

'''
Generate stats/plots on contig identity and alignment given a BAM of contigs VS true reference
'''

# read data indexes
READ_ID = 0
REVERSAL_STATUS = 1
REF_ALIGNMENT_START = 2
REF_ALIGNMENT_STOP = 3
ALIGNMENT_LENGTH = 4
READ_LENGTH = 5
CONTIG_LENGTH = 6
N_INITIAL_CLIPPED_BASES = 7
N_MATCHES = 8
N_TOTAL_MISMATCHES = 9
N_TOTAL_DELETES = 10
N_TOTAL_INSERTS = 11
IDENTITY = 12


def read_centromere_table(centromere_table_path, target_chromosome_name):
    """
    Read tsv file describing centromere locations with format convention:
        #bin	chrom	chromStart	chromEnd	name
        23	    chr1	122503247	124785432	GJ212202.1
        ...

    :param centromere_table_path:
    :param chromosome_name:
    :return:
    """
    centromere_coordinates = list()

    with open(centromere_table_path, "r") as file:
        for l,line in enumerate(file):
            if l == 0:
                continue
            if line.isspace():
                continue

            line = line.strip().split("\t")
            _, chromosome_name, start, end, event = line

            start = int(start)
            end = int(end)

            if chromosome_name == target_chromosome_name:
                centromere_coordinates.append([start,end])

    return centromere_coordinates


def read_gap_table(table_path, target_chromosome_name, size_cutoff=0):
    """
    Read tsv file describing gap locations with format convention:
        #bin	chrom	chromStart	chromEnd	...
        23	    chr1	122503247	124785432	...
        ...

    :param table_path:
    :param target_chromosome_name:
    :return:
    """
    coordinates = list()

    with open(table_path, "r") as file:
        for l,line in enumerate(file):
            if l == 0:
                continue
            if line.isspace():
                continue

            line = line.strip().split("\t")
            chromosome_name = line[1]
            start = line[2]
            end = line[3]
            type = line[-2]

            start = int(start)
            end = int(end)
            size = end - start

            if chromosome_name == target_chromosome_name:
                if size > size_cutoff:
                    coordinates.append([start,end])

    return coordinates


def plot_contig_blocks(axes, y, ref_alignment_start, n_initial_clipped_bases, alignment_length, scale, color, contig_length, draw_divider):
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

    last_match_position = x + width

    rect = patches.Rectangle((x, y), width, height, color=color)
    axes.add_patch(rect)

    # ---- If the block is a supplementary extension to another block ----
    # ---- draw a line to show this ----

    if draw_divider:
        print("supplementary!", x, y)
        axes.plot([x,x], [y-0.1,y+height+0.1], linestyle=(0, (1, 1)), color=[0.5,0.5,0.5])

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

    return last_match_position


def plot_windows(figure, axes, coordinates, scale, y_min, y_max, color):
    """
    add a vertical patch to the plot showing where the centromere falls (supposedly, according to UCSC table browser)
    :param figure:
    :param axes:
    :param centromere_coordinates:
    :param scale:
    :param y_min:
    :param y_max:
    :return:
    """

    for coordinate in coordinates:
        x_min, x_max = coordinate

        x_min = x_min/scale
        x_max = x_max/scale

        width = x_max - x_min
        height = y_max - y_min

        rect = patches.Rectangle((x_min, y_min), width, height, color=color, zorder=0)
        axes.add_patch(rect)


def plot_contigs(output_dir, read_data, chromosome_name, chromosome_length, total_identity, bam_path,
                 centromere_coordinates, gap_coordinates, segdup_coordinates, y_min=None, y_max=None,
                 save=True, show=False, group_supplementaries=False):
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

    read_data = sorted(read_data, key=lambda x: x[REF_ALIGNMENT_START])

    scale = 1000000

    if chromosome_length < scale:
        scale = chromosome_length

    left_column_offset = -0.5*chromosome_length/scale
    column_width = 0.5*chromosome_length/scale

    # ---- plot reference chromosome and REF text legend ----

    width = chromosome_length/scale
    height = 0.2

    rect = patches.Rectangle((0, 0), width, height, color="gray")
    axes.add_patch(rect)

    x = -(chromosome_length / scale) + left_column_offset
    pyplot.text(x, 0, "REF", weight="bold")

    forward_contig_positions = defaultdict(list)
    reverse_contig_positions = defaultdict(list)

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
            contig_positions = reverse_contig_positions
        else:
            contig_positions = forward_contig_positions

        found_suitable_lane = False
        l = 0

        if group_supplementaries:
            if read_id in contig_positions:
                for l,lane in enumerate(contig_positions[read_id]):
                    x,y = lane

                    x_current = ref_alignment_start / scale

                    if x_current > x:
                        position = [x_current, y]
                        contig_positions[read_id][l] = position
                        found_suitable_lane = True
                        break

        if reversal_status:
            color = [77.6, 39.6, 7.1]
            color = [c/100 for c in color]

            if not found_suitable_lane:
                n_reverse += 1
                y = - n_reverse/2

        else:
            color = [4.3, 45.1, 47.1]
            color = [c/100 for c in color]

            if not found_suitable_lane:
                n_forward += 1
                y = n_forward/2

        last_match_position = plot_contig_blocks(axes=axes,
                                                 y=y,
                                                 ref_alignment_start=ref_alignment_start,
                                                 n_initial_clipped_bases=n_initial_clipped_bases,
                                                 alignment_length=alignment_length,
                                                 scale=scale,
                                                 color=color,
                                                 contig_length=contig_length,
                                                 draw_divider=False)

        # print(l, contig_positions)
        if not found_suitable_lane:
            contig_positions[read_id].append([last_match_position,y])

            # ---- plot divider lines and text data ----

        # axes.axhline(y + 0.33, linestyle="--", linewidth=0.6)

        x = -(chromosome_length/scale) + left_column_offset
        pyplot.text(x, y, read_id)

        x = -(chromosome_length/scale) + left_column_offset + column_width
        pyplot.text(x, y, str(round(identity*100, 3)))

    # ---- plot leftover divider lines ----

    # axes.axhline(-(n_reverse+1)/2 + 0.33, linestyle="--", linewidth=0.6)
    # axes.axhline(0 + 0.33, linestyle="--", linewidth=0.6)

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

    # ---- fig size ----

    figure.set_size_inches(max(6, (x_max-x_min)/2 + 1), (y_max - y_min)/2 + 1)

    # ---- chromosome title ----

    pyplot.title("Chromosome: "+chromosome_name)

    if centromere_coordinates is not None:
        plot_windows(figure=figure,
                     axes=axes,
                     coordinates=centromere_coordinates,
                     scale=scale,
                     y_min=y_min,
                     y_max=y_max,
                     color=[0.9, 0.9, 0.9])

    if gap_coordinates is not None:
        plot_windows(figure=figure,
                     axes=axes,
                     coordinates=gap_coordinates,
                     scale=scale,
                     y_min=y_min,
                     y_max=y_max,
                     color=[0.7, 0.9, 0.7])

    if segdup_coordinates is not None:
        plot_windows(figure=figure,
                     axes=axes,
                     coordinates=segdup_coordinates,
                     scale=scale,
                     y_min=y_min,
                     y_max=y_max,
                     color=[0.9, 0.7, 0.7])

    if save:
        filename = os.path.basename(bam_path+"_"+chromosome_name)
        filename_prefix = ".".join(filename.split(".")[:-1])
        output_filename = filename_prefix + ".png"
        output_file_path = os.path.join(output_dir, output_filename)
        pyplot.axis("off")

        figure.savefig(output_file_path)

    if show:
        pyplot.show()
        pyplot.close()

    return figure, axes


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
    Parse through a cigar operation to find mismatches, inserts, deletes
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

    n_matches = 0
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

        n_matches = length - n_mismatches

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

    return ref_index_increment, read_index_increment, n_matches, n_mismatches, n_deletes, n_inserts


def get_read_stop_position(read):
    """
    Returns the stop position of the reference to where the read stops aligning
    :param read: The read
    :return: stop position of the reference where the read last aligned
    """
    ref_alignment_stop = read.reference_end

    # only find the position if the reference end is fetched as none from pysam API
    # From pysam docs:
    #   "None values will be included for any soft-clipped or unaligned positions within the read."
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

    n_secondary = 0

    for read in reads:
        if read.is_secondary:
            n_secondary += 1

        if read.mapping_quality > 0 and not read.is_secondary:
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

            n_total_matches = 0
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
                ref_index_increment, read_index_increment, n_matches, n_mismatches, n_deletes, n_inserts = \
                    parse_cigar_tuple(cigar_code=cigar_code,
                                      length=length,
                                      alignment_position=ref_alignment_start + ref_index,
                                      read_sequence=read_sequence_segment,
                                      ref_sequence=ref_sequence_segment)

                # increase the read index iterator
                read_index += read_index_increment
                ref_index += ref_index_increment
                n_total_matches += n_matches
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
                    ref_alignment_stop,
                    ref_length,
                    read_length,
                    contig_length,
                    n_initial_clipped_bases,
                    n_total_matches,
                    n_total_mismatches,
                    n_total_deletes,
                    n_total_inserts,
                    identity]

            read_data.append(data)

    return read_data


def print_read_summaries(read_data, total_identity, chromosome_name, chromosome_length):
    print("chromosome_name:\t", chromosome_name)
    print("chromosome_length:\t", chromosome_length)

    total_forward_alignment_length = 0
    total_reverse_alignment_length = 0

    for data in read_data:
        read_id, reversal_status, ref_alignment_start, ref_alignment_stop, alignment_length, read_length, contig_length, n_initial_clipped_bases, n_total_matches, n_total_mismatches, n_total_deletes, n_total_inserts, identity = data
        print()
        print(read_id)
        print("reversed:\t", reversal_status)
        print("alignment_start:\t", ref_alignment_start)
        print("alignment_stop:\t\t", ref_alignment_stop)
        print("alignment_length:\t", alignment_length)
        print("read_length:\t\t", read_length)
        print("n_initial_clipped_bases:", n_initial_clipped_bases)
        print("n_total_matches:\t", n_total_matches)
        print("n_total_mismatches:\t", n_total_mismatches)
        print("n_total_deletes:\t", n_total_deletes)
        print("n_total_inserts:\t", n_total_inserts)
        print("identity:\t", identity)

        if not reversal_status:
            total_forward_alignment_length += alignment_length
        else:
            total_reverse_alignment_length += alignment_length

    print("\nTOTAL IDENTITY:\t", total_identity)

    print("\nTOTAL FORWARD ALIGNMENT LENGTH:\t", total_forward_alignment_length)
    print("ROUGH FORWARD COVERAGE ESTIMATE:\t", total_forward_alignment_length / chromosome_length)

    print("\nTOTAL REVERSE ALIGNMENT LENGTH:\t", total_reverse_alignment_length)
    print("ROUGH REVERSE COVERAGE ESTIMATE:\t", total_reverse_alignment_length / chromosome_length)


def export_summaries_to_csv(read_data, total_identity, chromosome_length, output_dir, bam_path, chromosome_name):
    total_forward_alignment_length = 0
    total_reverse_alignment_length = 0
    total_forward_matches = 0
    total_reverse_matches = 0
    total_forward_mismatches = 0
    total_reverse_mismatches = 0
    total_forward_inserts = 0
    total_reverse_inserts = 0
    total_forward_deletes = 0
    total_reverse_deletes = 0

    csv_rows = list()
    csv_rows.append(["contig_name"])
    csv_rows.append(["reversal_status"])
    csv_rows.append(["ref_alignment_start"])
    csv_rows.append(["ref_alignment_stop"])
    csv_rows.append(["alignment_length"])
    csv_rows.append(["read_length"])
    csv_rows.append(["contig_length"])
    csv_rows.append(["n_initial_clipped_bases"])
    csv_rows.append(["n_total_matches"])
    csv_rows.append(["n_total_mismatches"])
    csv_rows.append(["n_total_deletes"])
    csv_rows.append(["n_total_inserts"])

    for data in sorted(read_data, key=lambda x: x[REF_ALIGNMENT_START]):
        for i in range(IDENTITY):
            csv_rows[i].append(data[i])

        if not data[REVERSAL_STATUS]:
            total_forward_alignment_length += data[ALIGNMENT_LENGTH]
            total_forward_matches += data[N_MATCHES]
            total_forward_mismatches += data[N_TOTAL_MISMATCHES]
            total_forward_inserts += data[N_TOTAL_INSERTS]
            total_forward_deletes += data[N_TOTAL_DELETES]

        else:
            total_reverse_alignment_length += data[ALIGNMENT_LENGTH]
            total_reverse_matches += data[N_MATCHES]
            total_reverse_mismatches += data[N_TOTAL_MISMATCHES]
            total_reverse_inserts += data[N_TOTAL_INSERTS]
            total_reverse_deletes += data[N_TOTAL_DELETES]

    # Transpose
    csv_rows = list(map(list, zip(*csv_rows)))

    csv_rows.append(["total_reverse_matches",total_reverse_matches])
    csv_rows.append(["total_reverse_mismatches",total_reverse_mismatches])
    csv_rows.append(["total_reverse_inserts",total_reverse_inserts])
    csv_rows.append(["total_reverse_deletes",total_reverse_deletes])
    csv_rows.append(["total_forward_matches",total_forward_matches])
    csv_rows.append(["total_forward_mismatches",total_forward_mismatches])
    csv_rows.append(["total_forward_inserts",total_forward_inserts])
    csv_rows.append(["total_forward_deletes",total_forward_deletes])
    csv_rows.append(["total_identity",total_identity])
    csv_rows.append(["total_forward_alignment_length",total_forward_alignment_length])
    csv_rows.append(["forward_coverage_estimate",total_forward_alignment_length / chromosome_length])
    csv_rows.append(["total_reverse_alignment_length",total_reverse_alignment_length])
    csv_rows.append(["reverse_coverage_estimate",total_reverse_alignment_length / chromosome_length])

    total_forward_conventional_identity = total_forward_matches / (total_forward_matches + total_forward_mismatches + total_forward_inserts + total_forward_deletes)
    total_reverse_conventional_identity = total_reverse_matches / (total_reverse_matches + total_reverse_mismatches + total_reverse_inserts + total_reverse_deletes)

    csv_rows.append(["total_forward_conventional_identity", total_forward_conventional_identity])
    csv_rows.append(["total_reverse_conventional_identity", total_reverse_conventional_identity])

    filename = os.path.basename(bam_path+"_"+chromosome_name)
    filename_prefix = ".".join(filename.split(".")[:-1])
    output_filename = "summary_" + filename_prefix + "_" + chromosome_name + ".csv"
    output_file_path = os.path.join(output_dir, output_filename)

    print("\nSAVING CSV: %s" % output_file_path)

    with open(output_file_path, "w") as file:
        writer = csv.writer(file)
        for row in csv_rows:
            writer.writerow(row)


def get_chromosome_data(bam_path, reference_path, chromosome_name, output_dir, centromere_table_path, gap_table_path, segdup_table_path):
    fasta_handler = FastaHandler(reference_path)
    bam_handler = BamHandler(bam_file_path=bam_path)

    chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

    start = 0
    stop = chromosome_length

    reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=start, stop=stop)

    read_data = parse_reads(reads=reads, fasta_handler=fasta_handler, chromosome_name=chromosome_name)

    total_weighted_identity = sum([x[ALIGNMENT_LENGTH] * x[IDENTITY] for x in read_data])
    total_alignment_bases = sum([x[ALIGNMENT_LENGTH] for x in read_data])

    # Calculate total identity, and approximate 0 if denominator is zero
    total_identity = total_weighted_identity / max(1e-9, total_alignment_bases)
    total_identity = round(total_identity, 6)

    # print_read_summaries(read_data=read_data,
    #                      total_identity=total_identity,
    #                      chromosome_name=chromosome_name,
    #                      chromosome_length=chromosome_length)

    export_summaries_to_csv(read_data=read_data,
                            total_identity=total_identity,
                            chromosome_length=chromosome_length,
                            output_dir=output_dir,
                            bam_path=bam_path,
                            chromosome_name=chromosome_name)

    if centromere_table_path is not None:
        centromere_coordinates = read_centromere_table(centromere_table_path=centromere_table_path,
                                                       target_chromosome_name=chromosome_name)
    else:
        centromere_coordinates = None

    if gap_table_path is not None:
        gap_coordinates = read_gap_table(table_path=gap_table_path,
                                         target_chromosome_name=chromosome_name)
    else:
        gap_coordinates = None

    if segdup_table_path is not None:
        segdup_coordinates = read_gap_table(table_path=segdup_table_path,
                                            target_chromosome_name=chromosome_name,
                                            size_cutoff=10000)
    else:
        segdup_coordinates = None

    figure, axes = plot_contigs(output_dir=output_dir,
                                read_data=read_data,
                                chromosome_name=chromosome_name,
                                chromosome_length=chromosome_length,
                                total_identity=total_identity,
                                bam_path=bam_path,
                                centromere_coordinates=centromere_coordinates,
                                gap_coordinates=gap_coordinates,
                                segdup_coordinates=segdup_coordinates,
                                show=False)

    pyplot.close(figure)


def process_bam(bam_path, reference_path, output_dir=None, centromere_table_path=None, gap_table_path=None, segdup_table_path=None, max_threads=None):
    """
    Find useful summary data from a bam that can be represented as a table of identities, and a plot of alignments
    :param bam_path: path to a bam containing contigs aligned to a true reference
    :param reference_path: the true reference that contigs were aligned to
    :param output_dir: where to save plots
    :return:
    """
    print("\n" + bam_path)

    if max_threads is None:
        max_threads = cpu_count() - 2

    if output_dir is None:
        output_dir = "plots/"

    FileManager.ensure_directory_exists(output_dir)

    fasta_handler = FastaHandler(reference_path)

    chromosome_names = fasta_handler.get_contig_names()

    arguments = list()

    for chromosome_name in chromosome_names:
        # chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        arguments.append([bam_path, reference_path, chromosome_name, output_dir, centromere_table_path, gap_table_path, segdup_table_path])

    if len(arguments) < max_threads:
        max_threads = len(arguments)

    with Pool(processes=max_threads) as pool:
        pool.starmap(get_chromosome_data, arguments)


def main(bam_path, reference_path, output_dir, centromere_table_path, gap_table_path, segdup_table_path):

    process_bam(bam_path=bam_path,
                reference_path=reference_path,
                output_dir=output_dir,
                centromere_table_path=centromere_table_path,
                gap_table_path=gap_table_path,
                segdup_table_path=segdup_table_path)


if __name__ == "__main__":
    '''
    Processes arguments
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
        "--centromeres",
        type=str,
        required=False,
        help="UCSC table browser output for centromere locations"
    )
    parser.add_argument(
        "--gap",
        type=str,
        required=False,
        help="UCSC table browser output for gap locations"
    )
    parser.add_argument(
        "--seg_dups",
        type=str,
        required=False,
        help="UCSC table browser output for segmental duplication locations"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )

    args = parser.parse_args()

    main(bam_path=args.bam, reference_path=args.ref, output_dir=args.output_dir, centromere_table_path=args.centromeres,
         gap_table_path=args.gap, segdup_table_path=args.seg_dups)
