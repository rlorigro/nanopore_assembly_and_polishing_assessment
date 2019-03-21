from handlers.FileManager import FileManager
from matplotlib import pyplot
from collections import defaultdict
import argparse
import numpy
import csv
import os


# read data indexes
READ_ID = 0
REVERSAL_STATUS = 1
REF_ALIGNMENT_START = 2
REF_ALIGNMENT_STOP = 3
ALIGNMENT_LENGTH = 4
READ_LENGTH = 5
CONTIG_LENGTH = 6
N_INITIAL_CLIPPED_BASES = 7
N_MATCHES = 8
N_TOTAL_MISMATCHES = 9
N_TOTAL_DELETES = 10
N_TOTAL_INSERTS = 11
IDENTITY = 12


def plot_distribution(step, bins, frequencies):
    axes = pyplot.axes()

    center = (bins[:-1] + bins[1:]) / 2

    axes.bar(center, frequencies, width=step, align="center")

    axes.set_xlabel("Identity (matches/length)")
    axes.set_ylabel("Frequency")

    axes.set_xticks(numpy.arange(0,1,0.05))

    pyplot.show()
    pyplot.close()


def write_chromosomal_summary_data_to_csv(summary_headers, summary_data, output_dir, sample_name=None):
    if sample_name is None:
        sample_name = FileManager.get_datetime_string()

    filename = "aggregate_summary_" + sample_name + ".csv"
    file_path = os.path.join(output_dir, filename)

    print("Saving aggregate data to: %s" % os.path.abspath(file_path))

    with open(file_path, "w") as file:
        writer = csv.writer(file)

        writer.writerow(summary_headers)

        for data in summary_data:
            writer.writerow(data)


def get_chromosome_name_from_path(path, chromosome_name_prefix="chr"):
    """
    Split a filename based on assumptions regarding the chromosome name:
        - All names have a common prefix
        - Chromosome name has been appended to the filename as the final suffix before the file type extension ".csv"
    :param path:
    :return:
    """
    filename = os.path.basename(path)
    filename = filename.replace(".csv", "")

    chromosome_name_index = filename.find(chromosome_name_prefix)
    chromosome_name = filename[chromosome_name_index:]

    return chromosome_name


def aggregate_summary_data(summary_file_paths):
    """
    Read each file and
    :param summary_file_paths:
    :return:
    """
    header_length = None
    summary_headers = None
    summary_data = list()

    identities = list()

    global_summary_dict = defaultdict(int)

    for p,path in enumerate(summary_file_paths):
        filename = os.path.basename(path)
        chromosome_name = get_chromosome_name_from_path(path)

        file_summary_headers = ["chromosome_name", "summary_filename"]  # These are not part of the file contents
        file_summary_data = [chromosome_name, filename]

        print("Parsing file: %s" % filename)

        with open(path, "r") as file:
            reader = csv.reader(file)

            for l,line in enumerate(reader):
                if l == 0:
                    # Header line
                    header_length = len(line)

                elif len(line) == header_length:
                    # Raw data
                    n_matches = int(line[N_MATCHES])
                    n_mismatches = int(line[N_TOTAL_MISMATCHES])
                    n_deletes = int(line[N_TOTAL_DELETES])
                    n_inserts = int(line[N_TOTAL_INSERTS])
                    identity = n_matches / (n_matches + n_mismatches + n_deletes + n_inserts)

                    identities.append(identity)

                if len(line) != header_length:
                    # Get the names of the summary data from the first file
                    if p == 0:
                        file_summary_headers.append(line[0])

                    # Summary data (the good part)
                    file_summary_data.append(line[-1])
                    global_summary_dict[line[0]] += int(float((line[-1].strip())))

        if summary_headers is None:
            summary_headers = file_summary_headers

        summary_data.append(file_summary_data)

    # calculate global stats
    g_n_matches = int(global_summary_dict['total_reverse_matches'] + global_summary_dict['total_forward_matches'])
    g_n_mismatches = int(global_summary_dict['total_reverse_mismatches'] + global_summary_dict['total_forward_mismatches'])
    g_n_deletes = int(global_summary_dict['total_reverse_deletes'] + global_summary_dict['total_forward_deletes'])
    g_n_inserts = int(global_summary_dict['total_reverse_inserts'] + global_summary_dict['total_forward_inserts'])
    g_identity = 1.0 * g_n_matches / (g_n_matches + g_n_mismatches + g_n_deletes + g_n_inserts)
    global_summary_dict['total_identity'] = g_identity
    global_summary_dict['forward_coverage_estimate'] = -1
    global_summary_dict['reverse_coverage_estimate'] = -1

    global_summary_data = ['genome', 'genome']
    for k in summary_headers:
        global_summary_data.append(global_summary_dict[k])
    summary_data.append(global_summary_data)

    return summary_headers, summary_data, identities


def filter_decoys_from_paths(summary_file_paths):
    filtered_paths = list()

    for path in summary_file_paths:
        chromosome_name = get_chromosome_name_from_path(path)

        # if chromosome name contains no nonsense (is not a decoy), then let it pass
        if len(chromosome_name.split("_")) == 1:
            filtered_paths.append(path)

    return filtered_paths


def get_ordering(data, prefix):
    ordering = list()

    suffix = get_chromosome_suffix(data, prefix)

    if suffix.isdigit():
        index = int(suffix)
    else:
        index = -1

    ordering.append(index)
    ordering.append(suffix)
    ordering.append(data[1])

    print(suffix, ordering)
    return ordering


def get_chromosome_suffix(data, prefix):
    chromosome_name = data[0]
    suffix = chromosome_name.split(prefix)[-1]

    return suffix


def sort_summary_data(summary_data, prefix="chr"):
    summary_data = sorted(summary_data, key=lambda x: get_ordering(x, prefix))

    return summary_data


def main(summary_dir, output_dir, filter_decoys):
    FileManager.ensure_directory_exists(output_dir)

    summary_file_paths = FileManager.get_file_paths_from_directory(summary_dir)

    if filter_decoys:
        print("Filtering decoy chromosomes")
        summary_file_paths = filter_decoys_from_paths(summary_file_paths)

    summary_headers, summary_data, identities = aggregate_summary_data(summary_file_paths)

    # ---- Plotting ----

    # step = 0.005        # bin size
    # max_length = 1      # end of histogram
    #
    # bins = numpy.arange(0, max_length + step, step=step)
    # frequencies, _ = numpy.histogram(identities, bins=bins)
    #
    # plot_distribution(step=step, bins=bins, frequencies=frequencies)
    #
    # summary_data = sort_summary_data(summary_data)

    # ------------------

    sample_name = os.path.basename(summary_dir.rstrip('/'))  # replace this with sample name extractor function?
    if sample_name == '': sample_name=None

    write_chromosomal_summary_data_to_csv(summary_headers=summary_headers,
                                          summary_data=summary_data,
                                          sample_name=sample_name,
                                          output_dir=output_dir)


def stringAsBool(s):
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
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", stringAsBool)  # add type keyword to registries

    parser.add_argument(
        "--summary_dir", "-i",
        type=str,
        required=True,
        help="directory where the output of get_summary_stats.py is located"
    )
    parser.add_argument(
        "--output_dir", "-o",
        type=str,
        required=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )
    parser.add_argument(
        "--filter_decoys", "-f",
        type=stringAsBool,
        required=False,
        default=False,
        help="desired output directory path (will be created during run time if doesn't exist)"
    )

    args = parser.parse_args()

    main(summary_dir=args.summary_dir, output_dir=args.output_dir, filter_decoys=args.filter_decoys)

