from summarize_contig_identities_and_alignments import get_read_stop_position
from multiprocessing import Manager, Pool, cpu_count
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from handlers.BamHandler import BamHandler
from matplotlib import pyplot
import argparse
import numpy
import os

NAME = 0
LENGTH = 1


def get_aligned_contig_lengths(bam_path, aligned_assembly_contigs):
    bam_handler = BamHandler(bam_file_path=bam_path)

    reads = bam_handler.get_reads(chromosome_name=None, start=None, stop=None)

    aligned_lengths = list()

    n_secondary = 0
    for read in reads:
        if read.is_secondary:
            n_secondary += 1

        if read.mapping_quality > 5 and not read.is_secondary:
            read_id = read.query_name
            ref_alignment_start = read.reference_start
            ref_alignment_stop = get_read_stop_position(read)
            ref_length = ref_alignment_stop - ref_alignment_start

            aligned_lengths.append([read_id, ref_length])

            print(read_id, ref_length)

    aligned_lengths = sorted(aligned_lengths, key=lambda x: x[LENGTH], reverse=True)

    aligned_assembly_contigs[bam_path] = aligned_lengths


def get_contig_lengths(assembly_path, assembly_contigs):
    handler = FastaHandler(assembly_path)

    contig_names = handler.get_contig_names()
    contigs = list()

    for name in sorted(contig_names):
        length = handler.get_chr_sequence_length(name)

        contigs.append([name, length])

    contigs = sorted(contigs, key=lambda x: x[LENGTH], reverse=True)

    print("Assembly parsed: %s" % assembly_path)

    assembly_contigs[assembly_path] = contigs


def get_contig_lengths_from_quast_log(quast_log_path, assembly_contigs):
    assembly_path = None
    contigs = None

    with open(quast_log_path, "r") as file:
        for l,line in enumerate(file):
            if l == 0 and line[0] == ">":
                # Header 1
                pass

            if l == 1 and line[0] != ">":
                # Field 1
                assembly_path = line.strip("'[]\n")

            if l == 2 and line[0] == ">":
                # Header 2
                pass

            if l == 3 and line[0] != ">":
                # Field 2
                lengths = line.strip("[]\n").split(", ")
                lengths = list(map(int, lengths))
                lengths = sorted(lengths, reverse=True)
                lengths = [l for l in lengths if l > 0]

                names = ["ctg%d" % i for i in range(len(lengths))]

                contigs = zip(names,lengths)

                assembly_contigs[assembly_path] = contigs

    return assembly_contigs


def get_all_contig_lengths_from_quast_logs_dir(parent_directory, recursive=True):
    quast_log_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=parent_directory,
                                                             file_extension="NAx.txt",
                                                             recursive=recursive)

    assembly_contigs = dict()
    for path in sorted(quast_log_paths):
        assembly_contigs = get_contig_lengths_from_quast_log(quast_log_path=path, assembly_contigs=assembly_contigs)

    return assembly_contigs


def get_all_lengths(assembly_path, recursive=False):
    if os.path.isdir(assembly_path):
        assembly_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=assembly_path,
                                                                file_extension=".fasta",
                                                                recursive=recursive)

    else:
        assembly_paths = [assembly_path]

    manager = Manager()
    assembly_contigs = manager.dict()

    max_threads = max(1, cpu_count()-2)

    arguments = list()

    for path in assembly_paths:

        arguments.append([path, assembly_contigs])

    if len(arguments) < max_threads:
        print("Fewer jobs than threads")
        max_threads = len(arguments)

    print("Using %d threads..." % max_threads)

    with Pool(processes=max_threads) as pool:
        pool.starmap(get_contig_lengths, arguments)

    return assembly_contigs


def get_all_aligned_lengths(bam_path, recursive=False):
    if os.path.isdir(bam_path):
        bam_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=bam_path,
                                                           file_extension=".bam",
                                                           recursive=recursive)

        print(bam_paths)
    else:
        bam_paths = [bam_path]

    manager = Manager()
    assembly_contigs = manager.dict()

    max_threads = max(1, cpu_count()-2)

    arguments = list()

    for path in bam_paths:
        arguments.append([path, assembly_contigs])

    if len(arguments) < max_threads:
        print("Fewer jobs than threads")
        max_threads = len(arguments)

    print("Using %d threads..." % max_threads)

    with Pool(processes=max_threads) as pool:
        pool.starmap(get_aligned_contig_lengths, arguments)

    return assembly_contigs


def generate_ngx_plot(assembly_contigs, input_dir, genome_size=None, y_max=180, title="NGx", figure=None, axes=None):
    samples = ["03492", "03098", "02723", "02080", "02055", "01243", "01109", "00733", "24385", "24149", "24143", "CHM13", "hg38_no_alts"]

    colors = [(175/256.0,   48/256.0,   51/256.0),  # red
              (224/256.0,   99/256.0,   58/256.0),  # orange
              (215/256.0,   219/256.0,  84/256.0),  # yellow
              (110/256.0,   170/256.0,  100/256.0),  # light green
              (80/256.0,    180/256.0,  150/256.0),  # green
              (100/256.0,   189/256.0,  197/256.0),  # green-blue
              (0/256.0,     170/256.0,  231/256.0),  # turquoise
              (51/256.0,    87/256.0,   182/256.0),  # blue
              (37/256.0,    36/256.0,   93/256.0),  # indigo
              (95/256.0,    51/256.0,   139/256.0),  # purple
              (200/256.0,   53/256.0,   93/256.0),  # pink
              (224/256.0,   99/256.0,   58/256.0),
              (110/256.0,   170/256.0,  100/256.0)]

    alphas = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1.0, 0.3, 0.3, 0.3, 1.0, 1.0]
    zorders = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1]

    labels = {}

    # ---------------------------------------------------------------------------

    # samples = ["shasta", "wtdbg2", "canu", "flye"]
    #
    # colors = [(0.890,0.120,0.031),
    #           (0.999,0.696,0.031),  # (112/256, 37/256, 163/256)
    #           (0.039,0.463,0.58),
    #           (0.024,0.69,0.224)]
    #
    # zorders = [1,0,0,0]
    # alphas = [1,0.9,1,1]
    #
    # labels = {}

    # ---------------------------------------------------------------------------
    #
    # samples = ["shasta", "hifi"]
    #
    # colors = [(0.933,0.153,0.031),
    #           (112/256, 37/256, 163/256),
    #           (0.039,0.463,0.58),
    #           (0.024,0.69,0.224)]
    #
    # zorders = [1,1]
    # alphas = [1,1]
    #
    # labels = {}

    # ---------------------------------------------------------------------------

    # samples = ["assembly_GM24385",
    #            "assembly_HG00733",
    #            "scaffold_GM24385",
    #            "scaffold_HG00733"]
    #
    # labels = {}
    #
    # colors = [(51/256.0,    87/256.0,   182/256.0),     # blue
    #           (51/256.0,    87/256.0,   182/256.0),     # green-blue
    #           # (200/256.0,   200/256.0,  200/256.0),     # grey
    #           (100/256.0,   189/256.0,  197/256.0),      # orange
    #           (100/256.0,   189/256.0,  197/256.0)]  # light green
    #
    # zorders = [1,1,1,1]
    # alphas = [0.5,1,0.5,1]

    # ---------------------------------------------------------------------------

    if genome_size is None:
        print("WARNING: genome_size unspecified, using human as default")
        genome_size = 3.23 * 1000**3

    if y_max is None:
        print("WARNING: y_max unspecified, using 180Mbp as default")
        y_max = 180

    if figure is None and axes is None:
        figure = pyplot.figure()
        axes = pyplot.axes()

    legend_names = list()
    for path,contigs in sorted(assembly_contigs.items(), key=lambda x: x[0]):
        print("Plotting assembly: %s" % path)

        sample_matched = False
        for name in samples:
            if name.lower() in path.lower():
                sample_index = samples.index(name)
                color = colors[sample_index]
                alpha = alphas[sample_index]
                zorder = zorders[sample_index]
                sample_name = name
                sample_matched = True

        if not sample_matched:
            print("ERROR: color not found for %s" % path)
            sample_index = 0
            color = colors[sample_index]
            alpha = alphas[sample_index]
            zorder = zorders[sample_index]
            sample_name = os.path.basename(path).split(".")[0]

        if sample_name in labels:
            label = labels[sample_name]
        else:
            label = sample_name

        x1 = 0
        y_prev = None

        x_coords = list()
        y_coords = list()

        for contig in contigs:
            y = contig[LENGTH]
            width = contig[LENGTH]/genome_size
            x2 = x1 + width

            if y_prev is not None:
                x_coords.extend([x1, x1])
                y_coords.extend([y_prev, y])

            x_coords.extend([x1,x2])
            y_coords.extend([y,y])

            x1 = x2
            y_prev = y

        if y_coords[-1] != 0:
            y_coords.append(0)
            x_coords.append(x_coords[-1])

        dashes = [1,0,1,0]

        if "hifi" in path.lower():
            label = "Canu CCS"

        if "shasta" in path:
            label = "Shasta Nanopore"

        if label not in legend_names:
            legend_names.append(label)

        axes.plot(x_coords, y_coords, color=color, alpha=alpha, zorder=zorder, dashes=dashes, linewidth=0.6)

    axes.legend(legend_names)

    axes.axvline(0.5, linestyle="--", alpha=0.3, linewidth=0.7, zorder=-1)

    max_size = y_max

    step_size = 20
    if step_size >= y_max:
        step_size = 1

    scale = 1_000_000

    axes.set_xlim([0,1])
    axes.set_ylim([0,max_size*scale])
    axes.set_yticks(numpy.arange(0,max_size+step_size,step_size)*scale)
    axes.set_yticklabels(numpy.arange(0,max_size+step_size,step_size))

    axes.set_title(title)
    axes.set_ylabel("Contig/scaffold size (Mbp)")
    axes.set_xlabel("Cumulative coverage")

    FileManager.ensure_directory_exists("output")

    output_dir = "output/"
    filename = input_dir.rstrip("/").split("/")[-1] + "_" + FileManager.get_datetime_string()
    file_path = os.path.abspath(os.path.join(output_dir, filename))

    print("SAVING FIGURE: %s" % file_path)
    figure.savefig(file_path + ".png", dpi=300)
    figure.savefig(file_path + ".pdf", dpi=300)

    # pyplot.show()
    pyplot.close()


def main(assembly_path, alignment_path, quast_path, genome_size, share_axes, y_max, recursive):
    if share_axes:
        figure = pyplot.figure()
        axes = pyplot.axes()
    else:
        figure = None
        axes = None

    if assembly_path is not None:
        assembly_contigs = get_all_lengths(assembly_path=assembly_path, recursive=recursive)

        generate_ngx_plot(assembly_contigs=assembly_contigs,
                          input_dir=assembly_path,
                          genome_size=genome_size,
                          figure=figure,
                          axes=axes,
                          y_max=y_max)

    if alignment_path is not None:
        assembly_contigs = get_all_aligned_lengths(bam_path=alignment_path, recursive=recursive)

        generate_ngx_plot(assembly_contigs=assembly_contigs,
                          input_dir=alignment_path,
                          genome_size=genome_size,
                          title="NGAx",
                          figure=figure,
                          axes=axes,
                          y_max=y_max)

    if quast_path is not None:
        assembly_contigs = get_all_contig_lengths_from_quast_logs_dir(parent_directory=quast_path, recursive=recursive)

        generate_ngx_plot(assembly_contigs=assembly_contigs,
                          input_dir=quast_path,
                          genome_size=genome_size,
                          title="NGAx",
                          figure=figure,
                          axes=axes,
                          y_max=y_max)

    if assembly_path is None and alignment_path is None and quast_path is None:
        exit("No input directory provided, please provide one of the following: \n"
             "\t- An assembly FASTA file path"
             "\t- A BAM alignment of assembly to ref ")


def string_as_bool(s):
    s = s.lower()
    boolean = None

    if s in {"t", "true", "1", "y", "yes"}:
        boolean = True
    elif s in {"f", "false", "0", "n", "no"}:
        boolean = False
    else:
        exit("Error: invalid argument specified for boolean flag: %s" % s)

    return boolean


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", string_as_bool)  # add type keyword to registries

    parser.add_argument(
        "--assembly",
        type=str,
        required=False,
        help="path of folder or file containing FASTA sequence file"
    )
    parser.add_argument(
        "--alignment",
        type=str,
        required=False,
        help="path of folder or file containing BAM alignment of assembly vs ref"
    )
    parser.add_argument(
        "--quast",
        type=str,
        required=False,
        help="path of parent directory containing quast logs of NAx/NGAx anywhere within"
    )
    parser.add_argument(
        "--genome_size",
        type=int,
        required=False,
        help="length of genome size for NGX"
    )
    parser.add_argument(
        "--y_max",
        type=int,
        required=False,
        help="Y limit of plot (ideally corresponding to Mbp of largest contig)"
    )
    parser.add_argument(
        "--share",
        type="bool",
        default="False",
        required=False,
        help="Whether to share the same axes for all plots"
    )
    parser.add_argument(
        "--recursive",
        type="bool",
        default="False",
        required=False,
        help="Whether to search all sub directories for relevant files"
    )
    args = parser.parse_args()

    main(assembly_path=args.assembly,
         alignment_path=args.alignment,
         genome_size=args.genome_size,
         quast_path=args.quast,
         share_axes=args.share,
         y_max=args.y_max,
         recursive=args.recursive)
