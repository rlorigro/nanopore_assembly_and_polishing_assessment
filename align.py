from modules.align import *
from handlers.FileManager import FileManager
import argparse


def main(ref_sequence_path, reads_sequence_path, output_dir=None, minimap_preset="map-ont"):
    if output_dir is None:
        output_dir = "./"
    else:
        FileManager.ensure_directory_exists(output_dir)

    reads_vs_ref_bam_path = align_minimap(output_dir=output_dir,
                                          ref_sequence_path=ref_sequence_path,
                                          reads_sequence_path=reads_sequence_path,
                                          preset=minimap_preset)


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
    parser.add_argument(
        "--minimap_preset",
        type=str,
        default="map-ont",
        choices=["map-ont", "asm5", "asm10", "asm20"],
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )

    args = parser.parse_args()

    main(reads_sequence_path=args.sequences, ref_sequence_path=args.ref, output_dir=args.output_dir, minimap_preset=args.minimap_preset)
