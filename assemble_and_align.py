from os.path import exists, dirname, basename, join
from handlers.FileManager import FileManager
from os import mkdir, walk
from subprocess import run, PIPE
import argparse


def assemble_wtdbg2(input_file_path, output_dir=None):
    print("\n-------- ASSEMBLING --------\n")

    input_filename_prefix = basename(input_file_path)
    input_filename_prefix = "_".join(input_filename_prefix.split(".")[:-1])

    output_filename_prefix = "_".join(["assembled", "wtdbg2", input_filename_prefix])

    print("output dir:\t%s" % output_dir)
    print("input:\t%s" % input_file_path)
    print("output:\t%s" % output_filename_prefix)

    arguments = ["wtdbg2", "-t", "30", "-x", "ont", "-i", input_file_path, "-o", output_filename_prefix]

    print("\nRUNNING: ", " ".join(arguments))
    run(arguments, cwd=output_dir)

    input_filename = output_filename_prefix + ".ctg.lay.gz"
    output_filename = output_filename_prefix + ".fa"
    arguments = ["wtpoa-cns", "-t", "30", "-i", input_filename, "-f", "-o", output_filename]

    print("\nRUNNING: ", " ".join(arguments))
    run(arguments, cwd=output_dir)

    return output_filename


def align_minimap(ref_sequence_filename, reads_sequence_path, output_dir=None):
    """
    Given a reference file and reads file align using minimap, generating a
    :param ref_sequence_filename:
    :param reads_sequence_path:
    :param output_dir:
    :return:
    """

    print("\n-------- ALIGNING --------\n")

    ref_sequence_filename_prefix = basename(ref_sequence_filename)
    ref_sequence_filename_prefix = "_".join(ref_sequence_filename_prefix.split(".")[:-1])

    input_filename_prefix = basename(reads_sequence_path)
    input_filename_prefix = "_".join(input_filename_prefix.split(".")[:-1])
    output_filename_prefix = input_filename_prefix + "_VS_" + ref_sequence_filename_prefix

    # ---- Minimap -----------

    output_filename = output_filename_prefix + ".sam"
    arguments = ["minimap2", "-a", "-t", "30", "-x", "map-ont", ref_sequence_filename, reads_sequence_path]

    print("\nRUNNING: ", " ".join(arguments))
    output_file_path = join(output_dir, output_filename)
    with open(output_file_path, "w") as output_file:
        print("REDIRECTING TO: ", output_file_path, "\n")
        run(arguments, cwd=output_dir, stdout=output_file, check=True)

    # ---- Sort SAM ----------

    input_filename = output_filename
    output_filename = output_filename_prefix + ".sorted.sam"
    arguments = ["samtools", "sort", input_filename, "-@", "30", "-O", "SAM", "-o", output_filename]

    print("\nRUNNING: ", " ".join(arguments))
    run(arguments, cwd=output_dir, check=True)

    # ---- Convert to BAM ----

    input_filename = output_filename
    output_filename = output_filename_prefix + ".sorted.bam"
    output_file_path = join(output_dir, output_filename)
    arguments = ["samtools", "view", input_filename, "-O", "BAM", "-@", "30"]

    print("\nRUNNING: ", " ".join(arguments))
    with open(output_file_path, "w") as output_file:
        print("REDIRECTING TO: ", output_file_path, "\n")
        run(arguments, cwd=output_dir, stdout=output_file)

    # ---- Index --------------

    input_filename = output_filename
    arguments = ["samtools", "index", input_filename]

    print("\nRUNNING: ", " ".join(arguments))
    run(arguments, cwd=output_dir, check=True)

    return output_filename_prefix + ".sorted.sam"


def polish_racon(output_dir, reads_file_path, reads_vs_ref_sam_path, ref_sequence_path):
    print("\n-------- POLISHING --------\n")

    output_filename_prefix = basename(reads_file_path)
    output_filename_prefix = ".".join(output_filename_prefix.split(".")[:-1])

    output_filename = "polished_racon_" + output_filename_prefix + ".fasta"
    output_file_path = join(output_dir, output_filename)
    arguments = ["racon", "-t", "30", reads_file_path, reads_vs_ref_sam_path, ref_sequence_path]

    print("\nRUNNING: ", " ".join(arguments))
    print("REDIRECTING TO: ", output_file_path, "\n")

    with open(output_file_path, "w") as output_file:
        run(arguments, cwd=output_dir, stdout=output_file, check=True)

    return output_filename


def main(input_file_path, true_ref_sequence_path=None, output_dir=None):
    if output_dir is None:
        output_dir = "./"
    else:
        FileManager.ensure_directory_exists(output_dir)

    assembly_sequence_filename = assemble_wtdbg2(output_dir=output_dir,
                                                 input_file_path=input_file_path)

    reads_vs_ref_sam_filename = align_minimap(output_dir=output_dir,
                                              ref_sequence_filename=assembly_sequence_filename,
                                              reads_sequence_path=input_file_path)

    polished_ref_sequence_filename = polish_racon(output_dir=output_dir,
                                                  reads_file_path=input_file_path,
                                                  reads_vs_ref_sam_path=reads_vs_ref_sam_filename,
                                                  ref_sequence_path=assembly_sequence_filename)

    if true_ref_sequence_path is not None:
        assembled_vs_true_ref_sam_filename = align_minimap(output_dir=output_dir,
                                                           ref_sequence_filename=true_ref_sequence_path,
                                                           reads_sequence_path=assembly_sequence_filename)

        polished_vs_true_ref_sam_filename = align_minimap(output_dir=output_dir,
                                                          ref_sequence_filename=true_ref_sequence_path,
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
