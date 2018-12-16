from os.path import exists, dirname, basename, join, abspath
from multiprocessing import cpu_count
from subprocess import run, PIPE


def assemble_wtdbg2(input_file_path, output_dir=None, max_threads=None):
    if max_threads is None:
        max_threads = cpu_count() - 2
        max_threads = str(max_threads)

    input_file_path = abspath(input_file_path)

    print("\n-------- ASSEMBLING --------\n")

    # ---- Step 1 ----

    input_filename_prefix = basename(input_file_path)
    input_filename_prefix = "_".join(input_filename_prefix.split(".")[:-1])
    output_filename_prefix = "_".join(["assembled", "wtdbg2", input_filename_prefix])

    arguments = ["wtdbg2", "-t", max_threads, "-x", "ont", "-i", input_file_path, "-o", output_filename_prefix]

    print("output dir:\t%s" % output_dir)
    print("input:\t%s" % input_file_path)
    print("output:\t%s" % output_filename_prefix)
    print("\nRUNNING: ", " ".join(arguments))

    run(arguments, cwd=output_dir)

    # ---- Step 2 ----

    input_filename = output_filename_prefix + ".ctg.lay.gz"
    output_filename = output_filename_prefix + ".fa"

    arguments = ["wtpoa-cns", "-t", max_threads, "-i", input_filename, "-f", "-o", output_filename]

    print("\nRUNNING: ", " ".join(arguments))

    run(arguments, cwd=output_dir)

    output_file_path = abspath(join(output_dir, output_filename))

    return output_file_path
