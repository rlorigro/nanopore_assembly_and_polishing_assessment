from modules.align import *
from handlers.FileManager import FileManager
import argparse


def main(ref_sequence_path, reads_sequence_path, max_threads=None, output_dir=None, minimap_preset="map-ont", k=15):
    if output_dir is None:
        output_dir = "./"
    else:
        FileManager.ensure_directory_exists(output_dir)

    reads_vs_ref_bam_path = align_minimap(output_dir=output_dir,
                                          ref_sequence_path=ref_sequence_path,
                                          reads_sequence_path=reads_sequence_path,
                                          preset=minimap_preset,
                                          max_threads=max_threads,
                                          k=k)


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
        help="which of the minimap alignment presets to use: 'map-ont', 'asm5', 'asm10', 'asm20'"
    )
    parser.add_argument(
        "--k",
        type=int,
        default=15,
        required=False,
        help="what size k-mer to use for minimizers"
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=None,
        required=False,
        help="how many vCPU to allocate for minimap"
    )

    args = parser.parse_args()

    main(reads_sequence_path=args.sequences,
         ref_sequence_path=args.ref,
         output_dir=args.output_dir,
         minimap_preset=args.minimap_preset,
         max_threads=args.max_threads,
         k=args.k)
