from multiprocessing import Manager, Pool, cpu_count
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from matplotlib import pyplot
import argparse
import numpy

NAME = 0
LENGTH = 1


def get_contig_lengths(assembly_path, assembly_contigs):
    handler = FastaHandler(assembly_path)

    contig_names = handler.get_contig_names()
    contigs = list()

    for name in sorted(contig_names):
        length = handler.get_chr_sequence_length(name)
        # print(name, length)

        contigs.append([name, length])

    contigs = sorted(contigs, key=lambda x: x[LENGTH], reverse=True)

    print("Assembly parsed: %s" % assembly_path)

    assembly_contigs[assembly_path] = contigs


def main(assembly_dir):
    samples = ["03492", "03098", "02723", "02080", "02055", "01243", "01109", "00733", "24385", "24149", "24143"]

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
              (200/256.0,   53/256.0,   93/256.0)]  # pink

    alphas = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 1.0, 0.3, 0.3, 0.3]
    zorders = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

    genome_size = 3.3 * 1000**3

    # assembly_dir = "/Users/saureous/data/nanopore/agbt/assemblies"
    paths = FileManager.get_all_file_paths_by_type(parent_directory_path=assembly_dir, file_extension=".fasta")

    # print(paths)
    figure = pyplot.figure()
    axes = pyplot.axes()

    manager = Manager()
    assembly_contigs = manager.dict()

    # max_threads = max(1, cpu_count()-2)
    max_threads = 1

    arguments = list()

    for path in paths:
        arguments.append([path, assembly_contigs])

    if len(arguments) < max_threads:
        print("Fewer jobs than threads")
        max_threads = len(arguments)

    print("Using %d threads..." % max_threads)

    with Pool(processes=max_threads) as pool:
        pool.starmap(get_contig_lengths, arguments)

    legend_names = list()
    for path,contigs in assembly_contigs.items():
        print("Plotting assembly: %s" % path)

        sample_matched = False
        for name in samples:
            if name in path:
                sample_index = samples.index(name)
                color = colors[sample_index]
                alpha = alphas[sample_index]
                zorder = zorders[sample_index]
                sample_name = name
                sample_matched = True

        if not sample_matched:
            print("ERROR: color not found for %s" % path)
            sample_index = 4
            color = colors[sample_index]
            alpha = alphas[sample_index]
            zorder = zorders[sample_index]
            sample_name = "Unknown"

        legend_names.append(sample_name)

        # # color = "blue"

        x1 = 0
        y_prev = None

        x_coords = list()
        y_coords = list()

        for contig in contigs:
            y = contig[LENGTH]
            width = contig[LENGTH]/genome_size
            x2 = x1 + width

            if y_prev is not None:
                # axes.plot([x1, x1], [y_prev, y], color=color)
                x_coords.extend([x1, x1])
                y_coords.extend([y_prev, y])

            # axes.plot([x1,x2],[y,y], color=color)
            x_coords.extend([x1,x2])
            y_coords.extend([y,y])

            x1 = x2
            y_prev = y

        print(x_coords[:10])
        print(y_coords[:10])
        axes.plot(x_coords, y_coords, color=color, alpha=alpha, zorder=zorder)

    axes.legend(legend_names)

    max_size = 150
    axes.set_xlim([0,1])
    axes.set_ylim([0,100000000])
    axes.set_yticks(numpy.arange(0,max_size,20)*1000000)
    axes.set_yticklabels(numpy.arange(0,max_size,20))

    axes.set_title("NGx")
    axes.set_ylabel("Contig/scaffold size (Mbp)")
    axes.set_xlabel("Cumulative genome coverage")

    FileManager.ensure_directory_exists("output")
    figure.savefig("output/" + assembly_dir.split("/")[-1] + ".png", dpi=300)

    # pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir",
        type=str,
        required=True,
        help="file path of FASTQ or FASTA sequence file"
    )
    args = parser.parse_args()

    main(assembly_dir=args.input_dir)
