from modules.entropy import calculate_shannon_entropy, find_longest_repeat
from modules.sort import sort_chromosome_names, sort_summary_data
from handlers.BamHandler import BamHandler
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from collections import defaultdict
import argparse
import math
import csv
import os

'''
Generate stats/plots on contig identity and alignment given a BAM of contigs VS true reference
'''


# headers and indexes for indel data vectors
DATA_INDEXES = {"sequence_name": 0,
                "chromosome_name": 1,
                "cigar_type": 2,
                "ref_start": 3,
                "ref_stop": 4,
                "ref_allele": 5,
                "ref_allele_context": 6,
                "read_start": 7,
                "read_stop": 8,
                "read_allele": 9,
                "read_allele_context": 10,
                "reversal_status": 11,
                "ref_window": 12,
                "entropy": 13,
                "max_repeat": 14,
                "is_runlength_error": 15}


MISMATCH_INDEXES = {"ref_allele": 0,
                    "read_allele": 1,
                    "ref_start": 2,
                    "ref_stop": 3,
                    "read_start": 4,
                    "read_stop": 5}


def parse_match(ref_index, read_index, length, read_sequence, ref_sequence):
    """
    Process a cigar operation that is a match
    :param alignment_position: Position where this match happened
    :param read_sequence: Read sequence
    :param ref_sequence: Reference sequence
    :param length: Length of the operation
    :return:
    This method updates the candidates dictionary.
    """
    n_mismatches = 0
    mismatches = list()
    for i in range(0, length):
        allele = read_sequence[i]
        ref = ref_sequence[i]

        if allele != ref:
            # mismatch
            n_mismatches += 1
            mismatches.append([ref, allele, ref_index+i, ref_index+i+1, read_index+i, read_index+i+1])

    return n_mismatches, mismatches


def parse_cigar_tuple(cigar_code, length, ref_index, read_index, read_sequence, ref_sequence):
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

    mismatches = list()

    n_mismatches = 0
    n_deletes = 0
    n_inserts = 0

    # deal different kinds of operations
    if cigar_code == 0:
        # match
        n_mismatches, mismatches = parse_match(ref_index=ref_index,
                                               read_index=read_index,
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

    return ref_index_increment, read_index_increment, n_mismatches, n_deletes, n_inserts, mismatches


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


def parse_reads(reads, chromosome_name, fasta_handler, homopolymer_window_size=11):
    """
    Given a set of pysam read objects, generate data for matches/mismatches/inserts/deletes and contig size/position for
    each read
    :param reads: pysam aligned segment objects
    :param chromosome_name:
    :param fasta_handler: fasta_handler object that can retrieve substrings from the reference sequence
    :return:
    """
    left_pad = math.floor((homopolymer_window_size - 1)/2)
    right_pad = math.ceil((homopolymer_window_size - 1)/2) + 1

    inserts = defaultdict(list)
    deletes = defaultdict(list)
    mismatches = defaultdict(list)

    n_secondary = 0

    for read in reads:
        if read.is_secondary:
            n_secondary += 1
            # print(read.query_name, n_secondary)

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
                ref_index_increment, read_index_increment, n_mismatches, n_deletes, n_inserts, segment_mismatches = \
                    parse_cigar_tuple(cigar_code=cigar_code,
                                      length=length,
                                      ref_index=ref_index,
                                      read_index=read_index,
                                      read_sequence=read_sequence_segment,
                                      ref_sequence=ref_sequence_segment)

                if cigar_code == 0:
                    for mismatch in segment_mismatches:
                        # mismatch
                        cigar_type = "SNP"

                        ref_start = ref_alignment_start + mismatch[MISMATCH_INDEXES["ref_start"]]
                        ref_stop = ref_alignment_start + mismatch[MISMATCH_INDEXES["ref_stop"]]
                        read_start = mismatch[MISMATCH_INDEXES["read_start"]]
                        read_stop = mismatch[MISMATCH_INDEXES["read_stop"]]

                        ref_allele = mismatch[MISMATCH_INDEXES["ref_allele"]]
                        read_allele = mismatch[MISMATCH_INDEXES["read_allele"]]

                        left_index = mismatch[MISMATCH_INDEXES["ref_start"]] - left_pad
                        right_index = mismatch[MISMATCH_INDEXES["ref_start"]] + right_pad

                        left_index = max(0, left_index)
                        right_index = min(len(ref_sequence), right_index)

                        ref_window = ref_sequence[left_index:right_index]

                        entropy = round(calculate_shannon_entropy(ref_window),3)
                        max_repeat = find_longest_repeat(ref_window)

                        is_runlength_error = False

                        ref_allele_context = ref_sequence[mismatch[MISMATCH_INDEXES["ref_start"]] - 1:mismatch[MISMATCH_INDEXES["ref_start"]] + 2]
                        read_allele_context = read_sequence[mismatch[MISMATCH_INDEXES["read_start"]] - 1:mismatch[MISMATCH_INDEXES["read_start"]] + 2]

                        data = [chromosome_name, cigar_type, ref_start, ref_stop, ref_allele, ref_allele_context, read_start, read_stop,
                                read_allele, read_allele_context, reversal_status, ref_window, entropy, max_repeat, is_runlength_error]

                        mismatches[read_id].append(data)

                elif cigar_code == 1:
                    # insert
                    cigar_type = "INS"

                    ref_start = ref_alignment_start + ref_index
                    ref_stop = ref_alignment_start + ref_index + ref_index_increment
                    read_start = read_index
                    read_stop = read_index + read_index_increment

                    read_allele = read_sequence[read_start:read_stop]
                    ref_allele = ref_sequence[ref_index:ref_index + ref_index_increment]

                    left_index = max(0, ref_index - left_pad)
                    right_index = min(len(ref_sequence), ref_index + right_pad)

                    ref_window = ref_sequence[left_index:right_index]

                    entropy = round(calculate_shannon_entropy(ref_window), 3)
                    max_repeat = find_longest_repeat(ref_window)

                    is_runlength_error = False

                    characters = set(read_allele)
                    if len(characters) == 1:
                        if read_allele[0] == ref_sequence[ref_index-1] or read_allele[-1] == ref_sequence[ref_index]:
                            is_runlength_error = True

                    # print("INSERT")
                    # print("REF\t",ref_sequence[ref_index-1:ref_index + 1])
                    # print("READ\t", read_sequence[read_index-1:read_index+read_index_increment+1])
                    # print(is_runlength_error)
                    # print()

                    ref_allele_context = ref_sequence[ref_index-1:ref_index + 1]
                    read_allele_context = read_sequence[read_index-1:read_index+read_index_increment+1]

                    data = [chromosome_name, cigar_type, ref_start, ref_stop, ref_allele, ref_allele_context, read_start, read_stop,
                            read_allele, read_allele_context, reversal_status, ref_window, entropy, max_repeat, is_runlength_error]

                    inserts[read_id].append(data)

                elif cigar_code == 2 or cigar_code == 3:
                    # delete or refskip
                    cigar_type = "DEL"

                    ref_start = ref_alignment_start + ref_index
                    ref_stop = ref_alignment_start + ref_index + ref_index_increment
                    read_start = read_index
                    read_stop = read_index + read_index_increment

                    read_allele = read_sequence[read_start:read_stop]
                    ref_allele = ref_sequence[ref_index:ref_index + ref_index_increment]

                    left_index = max(0, ref_index - left_pad)
                    right_index = min(len(ref_sequence), ref_index + right_pad)

                    ref_window = ref_sequence[left_index:right_index]

                    entropy = round(calculate_shannon_entropy(ref_window), 3)
                    max_repeat = find_longest_repeat(ref_window)

                    is_runlength_error = False

                    characters = set(ref_allele)
                    if len(characters) == 1:
                        if ref_allele[0] == read_sequence[read_index-1] or ref_allele[-1] == read_sequence[read_stop]:
                            is_runlength_error = True

                    # print("DELETE")
                    # print("REF\t",ref_sequence[ref_index-1:ref_index+ref_index_increment+1])
                    # print("READ\t",read_sequence[read_start-1:read_stop+1])
                    # print(is_runlength_error)
                    # print()

                    ref_allele_context = ref_sequence[ref_index-1:ref_index+ref_index_increment+1]
                    read_allele_context = read_sequence[read_start-1:read_stop+1]

                    data = [chromosome_name, cigar_type, ref_start, ref_stop, ref_allele, ref_allele_context, read_start, read_stop,
                            read_allele, read_allele_context, reversal_status, ref_window, entropy, max_repeat, is_runlength_error]

                    deletes[read_id].append(data)

                # increase the read/ref index iterator
                read_index += read_index_increment
                ref_index += ref_index_increment
                n_total_mismatches += n_mismatches
                n_total_deletes += n_deletes
                n_total_inserts += n_inserts

    return inserts, deletes, mismatches


def merge_variants(inserts, deletes, mismatches):
    """
    :param inserts: dictionary of read_id:[coord1, coord2, ...]
    :param deletes: dictionary of read_id:[coord1, coord2, ...]
    :return:
    """
    keys = set(inserts.keys()) | set(deletes.keys() | set(mismatches.keys()))    # union of insert and delete read IDs

    variants = dict()
    for key in keys:
        if key in inserts:
            read_inserts = inserts[key]
        else:
            read_inserts = list()

        if key in deletes:
            read_deletes = deletes[key]
        else:
            read_deletes = list()

        if key in mismatches:
            read_mismatches = mismatches[key]
        else:
            read_mismatches = list()

        read_variants = read_deletes + read_inserts + read_mismatches
        read_variants = sorted(read_variants, key=lambda x: x[DATA_INDEXES["ref_start"]])   # Sort by start pos

        variants[key] = read_variants

    return variants


def write_hashed_lists_to_csv(data, output_path):
    if not output_path.endswith(".csv"):
        output_path += ".csv"

    with open(output_path, "w") as file:
        headers = [item[0] for item in sorted(DATA_INDEXES.items(), key=lambda x: x[1])]

        header_line = ",".join(headers) + "\n"
        file.write(header_line)

        for key in sort_chromosome_names(data.keys()):
            for element in data[key]:
                line = [str(key)] + list(map(str,element))  # prepend the key name
                line = ",".join(line) + "\n"

                file.write(line)

    return


def export_variants_to_csv(output_dir, chromosome_name, mismatches, inserts, deletes, merge=True):
    if merge:
        # Make one output only (combine inserts and deletes)
        indels = merge_variants(mismatches=mismatches, inserts=inserts, deletes=deletes)

        # Prepare file paths
        filename = "variants_" + chromosome_name + ".csv"
        output_path = os.path.join(output_dir, filename)
        print("Saving variants: %s" % output_path)

        write_hashed_lists_to_csv(data=indels, output_path=output_path)

    else:
        # Prepare file paths
        filename = "mismatches_" + chromosome_name + ".csv"
        output_path = os.path.join(output_dir, filename)
        print("Saving mismatches: %s" % output_path)

        write_hashed_lists_to_csv(data=mismatches, output_path=output_path)

        # Prepare file paths
        filename = "deletes_" + chromosome_name + ".csv"
        output_path = os.path.join(output_dir, filename)
        print("Saving deletes: %s" % output_path)

        write_hashed_lists_to_csv(data=deletes, output_path=output_path)

        # Prepare file paths
        filename = "inserts_" + chromosome_name + ".csv"
        output_path = os.path.join(output_dir, filename)
        print("Saving deletes: %s" % output_path)

        write_hashed_lists_to_csv(data=inserts, output_path=output_path)


def process_bam(bam_path, reference_path, output_dir=None):
    """
    Find useful summary data from a bam that can be represented as a table of identities, and a plot of alignments
    :param bam_path: path to a bam containing contigs aligned to a true reference
    :param reference_path: the true reference that contigs were aligned to
    :param output_dir: where to save plots
    :return:
    """
    print("\n" + bam_path)

    if output_dir is None:
        output_dir = "variants/"

    # Make a subdirectory to contain everything
    datetime_string = FileManager.get_datetime_string()
    output_subdirectory = "variants_" + datetime_string
    output_dir = os.path.join(output_dir, output_subdirectory)
    FileManager.ensure_directory_exists(output_dir)

    bam_handler = BamHandler(bam_file_path=bam_path)
    fasta_handler = FastaHandler(reference_path)

    chromosome_names = fasta_handler.get_contig_names()
    chromosome_names = sort_chromosome_names(names=chromosome_names, prefix="chr")

    print("ref contig names:", chromosome_names)

    for chromosome_name in chromosome_names:
        print("Parsing alignments for ref contig:", chromosome_name)

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        start = 0
        stop = chromosome_length

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=start, stop=stop)

        inserts, deletes, mismatches = parse_reads(reads=reads, fasta_handler=fasta_handler, chromosome_name=chromosome_name)

        export_variants_to_csv(output_dir=output_dir,
                               chromosome_name=chromosome_name,
                               mismatches=mismatches,
                               inserts=inserts,
                               deletes=deletes,
                               merge=True)


def main(bam_path, reference_path, output_dir):

    process_bam(bam_path=bam_path,
                reference_path=reference_path,
                output_dir=output_dir)


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
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )

    args = parser.parse_args()

    main(bam_path=args.bam, reference_path=args.ref, output_dir=args.output_dir)
