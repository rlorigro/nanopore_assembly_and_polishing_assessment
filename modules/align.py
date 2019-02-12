from os.path import exists, dirname, basename, join, abspath
from multiprocessing import cpu_count
from subprocess import run, PIPE


def align_minimap(ref_sequence_path, reads_sequence_path, max_threads=None, output_dir=None, preset="map-ont"):
    """
    Given a reference file and reads file align using minimap, generating a
    :param ref_sequence_path:
    :param reads_sequence_path:
    :param output_dir:
    :return:
    """

    if max_threads is None:
        max_threads = cpu_count() - 2
        max_threads = str(max_threads)

    ref_sequence_path = abspath(ref_sequence_path)
    reads_sequence_path = abspath(reads_sequence_path)

    print("\n-------- ALIGNING --------\n")

    ref_sequence_filename_prefix = basename(ref_sequence_path)
    ref_sequence_filename_prefix = "_".join(ref_sequence_filename_prefix.split(".")[:-1])

    input_filename_prefix = basename(reads_sequence_path)
    input_filename_prefix = "_".join(input_filename_prefix.split(".")[:-1])
    output_filename_prefix = input_filename_prefix + "_VS_" + ref_sequence_filename_prefix

    # ---- Minimap -----------

    output_filename = output_filename_prefix + ".sam"
    output_file_path = join(output_dir, output_filename)

    arguments = ["minimap2", "-a", "-t", max_threads, "-x", preset, ref_sequence_path, reads_sequence_path]

    print("\nRUNNING: ", " ".join(arguments))

    with open(output_file_path, "w") as output_file:
        print("REDIRECTING TO: ", output_file_path, "\n")
        run(arguments, cwd=output_dir, stdout=output_file, check=True)

    # ---- Sort SAM ----------

    input_filename = output_filename
    output_filename = output_filename_prefix + ".sorted.sam"

    arguments = ["samtools", "sort", input_filename, "-@", max_threads, "-O", "SAM", "-o", output_filename]

    print("\nRUNNING: ", " ".join(arguments))

    run(arguments, cwd=output_dir, check=True)

    # ---- Convert to BAM ----

    input_filename = output_filename
    output_filename = output_filename_prefix + ".sorted.bam"
    output_file_path = join(output_dir, output_filename)

    arguments = ["samtools", "view", input_filename, "-O", "BAM", "-@", max_threads]

    print("\nRUNNING: ", " ".join(arguments))

    with open(output_file_path, "w") as output_file:
        print("REDIRECTING TO: ", output_file_path, "\n")
        run(arguments, cwd=output_dir, stdout=output_file, check=True)

    # ---- Index --------------

    input_filename = output_filename
    arguments = ["samtools", "index", input_filename]

    print("\nRUNNING: ", " ".join(arguments))

    run(arguments, cwd=output_dir, check=True)

    output_sam_file_path = abspath(join(output_dir, output_filename_prefix + ".sorted.sam"))
    output_bam_file_path = abspath(join(output_dir, output_filename_prefix + ".sorted.bam"))

    return output_sam_file_path, output_bam_file_path
