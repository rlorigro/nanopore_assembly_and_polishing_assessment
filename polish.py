from modules.align import *
from modules.polish import *
from handlers.FileManager import FileManager
import argparse


def polish(reads_file_path, assembly_sequence_path, true_ref_sequence_path=None, output_dir=None, n_passes=False):
    if output_dir is None:
        output_dir = "./"
    else:
        FileManager.ensure_directory_exists(output_dir)

    if true_ref_sequence_path is not None:
        assembled_vs_true_ref_sam_path, assembled_vs_true_ref_bam_path = align_minimap(output_dir=output_dir,
                                                                                       ref_sequence_path=true_ref_sequence_path,
                                                                                       reads_sequence_path=assembly_sequence_path)

    polished_ref_paths = list()

    for i in range(n_passes):
        suffix = str(i+1) + "x"
        polish_output_dir = join(output_dir,suffix)
        FileManager.ensure_directory_exists(polish_output_dir)

        if i == 0:
            ref_sequence_path = assembly_sequence_path
        else:
            ref_sequence_path = polished_ref_paths[i-1]

        reads_vs_polished_ref_sam_path = align_minimap(output_dir=polish_output_dir,
                                                       ref_sequence_path=ref_sequence_path,
                                                       reads_sequence_path=reads_file_path,
                                                       sam_only=True)

        repolished_ref_sequence_path = polish_racon(output_dir=polish_output_dir,
                                                    reads_file_path=reads_file_path,
                                                    reads_vs_ref_sam_path=reads_vs_polished_ref_sam_path,
                                                    ref_sequence_path=ref_sequence_path,
                                                    suffix=suffix)

        remove(reads_vs_polished_ref_sam_path)  # delete last sam

        polished_ref_paths.append(repolished_ref_sequence_path)

        if true_ref_sequence_path is not None:
            repolished_vs_true_ref_sam_path, repolished_vs_true_ref_bam_path = \
                align_minimap(output_dir=polish_output_dir,
                              ref_sequence_path=true_ref_sequence_path,
                              reads_sequence_path=repolished_ref_sequence_path)


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contigs",
        type=str,
        required=True,
        help="FASTA file path of contigs to be polished"
    )
    parser.add_argument(
        "--sequences",
        type=str,
        required=True,
        help="file path of FASTQ or FASTA sequence file containing reads"
    )
    parser.add_argument(
        "--true_ref",
        type=str,
        required=False,
        help="FASTA file path of true reference to be compared against for QC"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )
    parser.add_argument(
        "--n_passes",
        type=int,
        default=1,
        required=False,
        help="Rerun polisher with previously polished reference as input for n times. Default=1"
    )

    args = parser.parse_args()

    polish(args.sequences, args.contigs, true_ref_sequence_path=args.true_ref, output_dir=args.output_dir, n_passes=args.n_passes)
