from modules.align import *
from modules.assemble import *
from modules.polish import *
from handlers.FileManager import FileManager
import argparse


def main(input_file_path, true_ref_sequence_path=None, output_dir=None):
    if output_dir is None:
        output_dir = "./"
    else:
        FileManager.ensure_directory_exists(output_dir)


    polished_ref_sequence_filename = polish_racon(output_dir=output_dir,
                                                  reads_file_path=input_file_path,
                                                  reads_vs_ref_sam_path=reads_vs_ref_sam_filename,
                                                  ref_sequence_path=assembly_sequence_filename)

    polished_vs_true_ref_sam_filename = align_minimap(output_dir=output_dir,
                                                      ref_sequence_path=true_ref_sequence_path,
                                                      reads_sequence_path=polished_ref_sequence_filename)


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
        "--true_ref",
        type=str,
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

    main(args.sequences, true_ref_sequence_path=args.true_ref, output_dir=args.output_dir)
