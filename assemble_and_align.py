from modules.align import *
from modules.assemble import *
from handlers.FileManager import FileManager
import argparse


def main(reads_file_path, genome_size=None, output_dir=None):
    if output_dir is None:
        output_dir = "./output/"
    else:
        FileManager.ensure_directory_exists(output_dir)

    if genome_size is None:
        genome_size = "3g"
        print("WARNING: genome size flag not specified, defaulting to human size (3g)")

    assembly_sequence_path = assemble_wtdbg2(output_dir=output_dir,
                                             input_file_path=reads_file_path,
                                             genome_size=genome_size)

    reads_vs_ref_sam_path, reads_vs_ref_bam_path = align_minimap(output_dir=output_dir,
                                                                 ref_sequence_path=assembly_sequence_path,
                                                                 reads_sequence_path=reads_file_path)


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences",
        type=str,
        required=True,
        help="file path of FASTQ or FASTA sequence file"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )
    parser.add_argument(
        "--genome_size", "-g",
        type=str,
        required=False,
        help="override value for wtdbg2's -g flag for expected genome size. Default is human size (3g)"
    )

    args = parser.parse_args()

    main(reads_file_path=args.sequences, output_dir=args.output_dir)
