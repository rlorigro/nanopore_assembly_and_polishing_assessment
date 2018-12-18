from os.path import exists, dirname, basename, join, abspath
from multiprocessing import cpu_count
from subprocess import run, PIPE


def polish_racon(output_dir, reads_file_path, reads_vs_ref_sam_path, ref_sequence_path, max_threads=None, suffix=""):
    if max_threads is None:
        max_threads = cpu_count() - 2
        max_threads = str(max_threads)

    # Get absolute path if path not already absolute
    reads_file_path = abspath(reads_file_path)
    reads_vs_ref_sam_path = abspath(reads_vs_ref_sam_path)
    ref_sequence_path = abspath(ref_sequence_path)

    print("\n-------- POLISHING --------\n")

    if len(suffix) > 0:
        suffix = "_" + suffix

    output_filename_prefix = basename(reads_file_path)
    output_filename_prefix = ".".join(output_filename_prefix.split(".")[:-1])
    output_filename = "polished_racon_" + output_filename_prefix + suffix + ".fasta"
    output_file_path = join(output_dir, output_filename)

    arguments = ["racon", "-t", max_threads, reads_file_path, reads_vs_ref_sam_path, ref_sequence_path]

    print("\nRUNNING: ", " ".join(arguments))
    print("REDIRECTING TO: ", output_file_path, "\n")

    with open(output_file_path, "w") as output_file:
        run(arguments, cwd=output_dir, stdout=output_file, check=True)

    return output_file_path
