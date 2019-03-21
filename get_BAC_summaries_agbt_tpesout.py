from summarize_contig_identities_and_alignments import parse_reads
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from handlers.BamHandler import BamHandler
from collections import defaultdict
from multiprocessing import Pool
import argparse
import csv
import os


'''
Generate stats on contig identity and alignment given a BAM of contigs VS true reference
'''

# read data indexes
CHROMOSOME_NAME = 0
READ_ID = 1
REVERSAL_STATUS = 2
REF_ALIGNMENT_START = 3
REF_ALIGNMENT_STOP = 4
ALIGNMENT_LENGTH = 5
READ_LENGTH = 6
CONTIG_LENGTH = 7
N_INITIAL_CLIPPED_BASES = 8
N_MATCHES = 9
N_TOTAL_MISMATCHES = 10
N_TOTAL_DELETES = 11
N_TOTAL_INSERTS = 12
IDENTITY = 13


def export_bac_data_to_csv(read_data, output_dir, bam_path):
    total_alignment_length = 0

    csv_rows = list()
    csv_rows.append(["chromosome_name"])
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
    csv_rows.append(["identity"])

    for data in read_data.values():
        print(data)
        for i in range(IDENTITY+1):
            csv_rows[i].append(data[i])

        total_alignment_length += data[ALIGNMENT_LENGTH]

    # total_matches = sum([x[N_MATCHES] for x in read_data])
    # total_mismatches = sum([x[N_TOTAL_MISMATCHES] for x in read_data])
    # total_deletes = sum([x[N_TOTAL_INSERTS] for x in read_data])
    # total_inserts = sum([x[N_TOTAL_DELETES] for x in read_data])
    #
    # weighted_identity = total_matches/(total_matches + total_mismatches + total_deletes + total_inserts)
    #
    # print("___________")
    # print(total_alignment_length)
    # print(weighted_identity)
    # print("___________")

    # Transpose
    csv_rows = list(map(list, zip(*csv_rows)))

    filename = os.path.basename(bam_path)
    filename_prefix = ".".join(filename.split(".")[:-1])
    output_filename = "summary_" + filename_prefix + ".csv"
    output_file_path = os.path.join(output_dir, output_filename)

    print("\nSAVING CSV: %s" % output_file_path)

    with open(output_file_path, "w") as file:
        writer = csv.writer(file)
        for row in csv_rows:
            writer.writerow(row)


def filter_supplementaries_by_largest(data_per_bac):
    filtered_data = list()

    for bac_name in data_per_bac:
        bac_data = data_per_bac[bac_name]

        max_score = 0
        max_index = None
        max_identity = None
        for d,data in enumerate(sorted(bac_data, key=lambda x: x[REF_ALIGNMENT_START])):

            identity = data[N_MATCHES]/(data[N_MATCHES]+data[N_TOTAL_MISMATCHES]+data[N_TOTAL_INSERTS]+data[N_TOTAL_DELETES])
            score = data[N_MATCHES]

            if score > max_score:
                max_index = d
                max_score = score
                max_identity = round(identity,3)

        filtered_data.append(bac_data[max_index] + [max_identity])

    print(filtered_data)

    return filtered_data

def aggregate_bac_data(data_per_bac):
    filtered_data = list()

    for bac_name in data_per_bac:
        bac_data = data_per_bac[bac_name]

        for d,data in enumerate(sorted(bac_data, key=lambda x: x[REF_ALIGNMENT_START])):

            identity = data[N_MATCHES]/(data[N_MATCHES]+data[N_TOTAL_MISMATCHES]+data[N_TOTAL_INSERTS]+data[N_TOTAL_DELETES])
            score = data[N_MATCHES]

            filtered_data.append(data + [identity])

    print(filtered_data)

    return filtered_data

def process_bam(bam_path, reference_path, output_dir=None):
    """
    Find useful summary data from a bam that can be represented as a table of identities/matches/mismatches/indels
    :param bam_path: path to a bam containing contigs aligned to a true reference
    :param reference_path: the true reference that contigs were aligned to
    :param output_dir: where to save stats
    :return:
    """
    if output_dir is None:
        output_dir = "stats/"

    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_handler = FastaHandler(reference_path)
    # bac_fasta_handler = FastaHandler(bac_path)

    chromosome_names = ref_fasta_handler.get_contig_names()
    # bac_names = bac_fasta_handler.get_contig_names()
    #
    # print(chromosome_names)
    # print(bac_names)

    bam_reads = dict()
    # data_per_bac = defaultdict(list)

    for chromosome_name in chromosome_names:
        chromosome_length = ref_fasta_handler.get_chr_sequence_length(chromosome_name)

        start = 0
        stop = chromosome_length

        ref_fasta_handler = FastaHandler(reference_file_path=reference_path)
        bam_handler = BamHandler(bam_file_path=bam_path)

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=start, stop=stop)

        read_data = parse_reads(reads=reads, fasta_handler=ref_fasta_handler, chromosome_name=chromosome_name)

        for read in read_data:
            bac_name = read[1]
            if bac_name in bam_reads:
                print("{} appeared in {} ({}) multiple times!".format(bac_name, bam_path, chromosome_name))
                bam_reads[bac_name] = None
            else:
                bam_reads[bac_name] = read

    unique_well_aligned_bac = list(filter(lambda x: bam_reads[x] is not None, bam_reads.keys()))
    filtered_bam_reads = {bac:bam_reads[bac] for bac in unique_well_aligned_bac}

    export_bac_data_to_csv(read_data=filtered_bam_reads,
                           output_dir=output_dir,
                           bam_path=bam_path)

    return filtered_bam_reads


def main(bam_glob, ref_glob, output_dir):
    import glob
    bams = glob.glob(bam_glob)
    if len(bams) == 0: raise Exception("No files matching {}".format(bam_glob))

    refs = glob.glob(ref_glob)
    if len(refs) == 0: raise Exception("No files matching {}".format(ref_glob))

    all_bam_bacs = dict()
    for bam_path in bams:
        print("Handling BAM {}".format(bam_path))

        # guess at ref glob
        our_ref = None
        for ref in refs:
            if os.path.basename(ref).replace(".fa", "") in bam_path:
                if our_ref is not None:
                    raise Exception("Multiple refs ({}, {}) match for {}".format(our_ref, ref, bam_path))
                our_ref = ref
        if our_ref is None: raise Exception("Found no ref of {} for {}".format(refs, bam_path))

        bam_output = process_bam(bam_path=bam_path,
                                 reference_path=our_ref,
                                 output_dir=output_dir)
        all_bam_bacs[bam_path] = bam_output

    all_present_bacs = None
    for bam_data in all_bam_bacs.values():
        if all_present_bacs is None:
            all_present_bacs = set(bam_data.keys())
        else:
            all_present_bacs = all_present_bacs.intersection(bam_data.keys())

    bac_to_skip = ['AC275306.1', 'AC277731.1', 'AC275449.1', 'AC275257.1',]
    for back in bac_to_skip:
        if back in all_present_bacs: all_present_bacs.remove(back)
    print("All present bacs: {}\n\n".format(all_present_bacs))

    idents = dict()

    output_for_all_present_bacs = dict()
    output_dir_for_apb = os.path.join(output_dir, "apb")
    if not os.path.isdir(output_dir_for_apb):
        os.mkdir(output_dir_for_apb)

    # bam bacs
    for bam_path in all_bam_bacs.keys():
        subset_bam_data = dict()
        for bac in all_present_bacs:
            subset_bam_data[bac] = all_bam_bacs[bam_path][bac]

        output_for_all_present_bacs[bam_path] = subset_bam_data

        export_bac_data_to_csv(read_data=subset_bam_data,
                               output_dir=output_dir_for_apb,
                               bam_path=bam_path)

        total_match = 0
        total_mismatch = 0
        total_insert = 0
        total_delete = 0
        for bam_data in subset_bam_data.values():
            total_match += bam_data[N_MATCHES]
            total_mismatch += bam_data[N_TOTAL_MISMATCHES]
            total_insert += bam_data[N_TOTAL_INSERTS]
            total_delete += bam_data[N_TOTAL_DELETES]
        identity = total_match / (total_match + total_mismatch + total_insert + total_delete)

        idents[bam_path] = identity

    for bam_path in sorted(list(idents.keys())):
        print("{} total iden: {}".format(bam_path, idents[bam_path]))

    pass





if __name__ == "__main__":
    '''
    Processes arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam_glob",
        type=str,
        required=True,
        help="bam glob files BAM file path of contigs aligned to true reference"
    )
    parser.add_argument(
        "--ref_glob",
        type=str,
        required=True,
        help="FASTA file glob that should (less the .fa) be present in the path of the bam because fuck you"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )
    args = parser.parse_args()

    main(bam_glob=args.bam_glob, ref_glob=args.ref_glob, output_dir=args.output_dir)
