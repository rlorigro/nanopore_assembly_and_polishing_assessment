from summarize_contig_identities_and_alignments import parse_reads, export_chromosome_summary_to_csv, export_genome_summary_to_csv
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from handlers.BamHandler import BamHandler
from multiprocessing import Manager, Pool, cpu_count
import argparse

'''
Generate stats on contig identity and alignment given a BAM of contigs VS true reference
'''


def get_chromosome_stats(genome_data, reference_path, chromosome_name, start, stop, output_dir, bam_path):
    fasta_handler = FastaHandler(reference_file_path=reference_path)
    bam_handler = BamHandler(bam_file_path=bam_path)

    chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

    reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=start, stop=stop)

    read_data, chromosome_data = parse_reads(reads=reads,
                                             chromosome_name=chromosome_name,
                                             chromosome_length=chromosome_length,
                                             fasta_handler=fasta_handler)

    genome_data.append(chromosome_data)

    export_chromosome_summary_to_csv(read_data=read_data,
                                     chromosome_data=chromosome_data,
                                     output_dir=output_dir,
                                     bam_path=bam_path,
                                     chromosome_name=chromosome_name)


def process_bam(bam_path, reference_path, max_threads, output_dir=None):
    """
    Find useful summary data from a bam that can be represented as a table of identities/matches/mismatches/indels
    :param bam_path: path to a bam containing contigs aligned to a true reference
    :param reference_path: the true reference that contigs were aligned to
    :param output_dir: where to save stats
    :return:
    """
    if output_dir is None:
        output_dir = "stats/"

    if max_threads is None:
        max_threads = max(1, cpu_count() - 2)

    process_manager = Manager()
    genome_data = process_manager.list()

    FileManager.ensure_directory_exists(output_dir)

    fasta_handler = FastaHandler(reference_path)

    chromosome_names = fasta_handler.get_contig_names()

    arguments = list()

    for chromosome_name in chromosome_names:
        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        start = 0
        stop = chromosome_length

        arguments.append([genome_data, reference_path, chromosome_name, start, stop, output_dir, bam_path])

    if len(arguments) < max_threads:
        print("Fewer jobs than threads")
        max_threads = len(arguments)

    print("Using %d threads..." % max_threads)

    with Pool(processes=max_threads) as pool:
        pool.starmap(get_chromosome_stats, arguments)

    print("genome_data", genome_data)

    export_genome_summary_to_csv(bam_path=bam_path, output_dir=output_dir, genome_data=genome_data)


def main(bam_path, reference_path, output_dir, max_threads):
    process_bam(bam_path=bam_path,
                reference_path=reference_path,
                output_dir=output_dir,
                max_threads=max_threads)


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
        "--max_threads", "-t",
        type=int,
        required=False,
        help="FASTA file path of true reference to be compared against"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )
    args = parser.parse_args()

    main(bam_path=args.bam, reference_path=args.ref, output_dir=args.output_dir, max_threads=args.max_threads)
