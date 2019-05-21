from handlers.FileManager import FileManager
from matplotlib import pyplot
from collections import defaultdict
import argparse
import numpy
import csv
import os
import seaborn as sns
import pandas as pd
import operator as op
import random
import glob
import sys
from scipy import stats


# read data indexes
READ_ID = 0
REVERSAL_STATUS = 1
REF_ALIGNMENT_START = 2
REF_ALIGNMENT_STOP = 3
ALIGNMENT_LENGTH = 4
READ_LENGTH = 5
CONTIG_LENGTH = 6
N_INITIAL_CLIPPED_BASES = 7
N_MATCHES = 9
N_TOTAL_MISMATCHES = 10
N_TOTAL_DELETES = 11
N_TOTAL_INSERTS = 12
IDENTITY = 14

# chromosome_name,
# read_id,
# reversal_status,
# ref_alignment_start,
# ref_alignment_stop,
# ref_length,
# read_length,
# contig_length,
# n_initial_clipped_bases,
# n_total_matches,
# n_total_mismatches,
# n_total_deletes,
# n_total_inserts,
# sequence_identity,
# alignment_identity


def plot_distribution(step, bins, frequencies, save_location=None, title=None):
    axes = pyplot.axes()

    center = (bins[:-1] + bins[1:]) / 2

    if type(frequencies) == dict:
        keys = list(frequencies.keys())
        keys.sort()
        for k in keys:
            axes.bar(center, frequencies[k], label=k, width=step, align="center", alpha=0.2)
        axes.legend()
    else:
        axes.bar(center, frequencies, width=step, align="center")


    axes.set_xlabel("Identity (matches/length)")
    axes.set_ylabel("Frequency")

    axes.set_xticks(numpy.arange(0,1,0.1))
    axes.set_xlim(.6, 1)

    if title is not None:
        pyplot.title(title)

    if save_location is not None:
        pyplot.savefig(save_location)

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


def aggregate_summary_data(summary_file_paths, args):
    """
    Read each file and
    :param summary_file_paths:
    :return:
    """
    header_length = None
    summary_headers = None
    summary_data = list()

    read_len_to_identity = list()
    identities = list()
    identities_per_file = defaultdict(list)
    read_lengths_per_file = defaultdict(list)

    global_summary_dict = defaultdict(int)

    onemb_reads = list()

    for p,path in enumerate(summary_file_paths):
        filename = os.path.basename(path)

        if "whole_genome" in filename: continue
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

                    if int(line[READ_LENGTH]) < args.min_read_length: continue

                    # Raw data
                    n_matches = int(line[N_MATCHES])
                    n_mismatches = int(line[N_TOTAL_MISMATCHES])
                    n_deletes = int(line[N_TOTAL_DELETES])
                    n_inserts = int(line[N_TOTAL_INSERTS])
                    identity = n_matches / (n_matches + n_mismatches + n_deletes + n_inserts)

                    # if int(line[READ_LENGTH]) >= 1000000:
                    #     line.append(identity)
                    #     onemb_reads.append(line)

                    identities.append(identity)
                    read_len_to_identity.append((int(line[READ_LENGTH]), identity))
                    identities_per_file[path].append(identity)
                    read_lengths_per_file[path].append(line[READ_LENGTH])

                if len(line) != header_length:
                    # Get the names of the summary data from the first file
                    if p == 0:
                        file_summary_headers.append(line[0])

                    # Summary data (the good part)
                    file_summary_data.append(line[-1])

                    #TODO
                    # global_summary_dict[line[0]] += int(float((line[-1].strip())))

        if summary_headers is None:
            summary_headers = file_summary_headers

        summary_data.append(file_summary_data)

    # if len(onemb_reads) != 0:
    #     print("\n1MB+ reads:")
    # for mb in onemb_reads:
    #     print("{}".format(mb))

    # calculate global stats
    g_n_matches = int(global_summary_dict['total_reverse_matches'] + global_summary_dict['total_forward_matches'])
    g_n_mismatches = int(global_summary_dict['total_reverse_mismatches'] + global_summary_dict['total_forward_mismatches'])
    g_n_deletes = int(global_summary_dict['total_reverse_deletes'] + global_summary_dict['total_forward_deletes'])
    g_n_inserts = int(global_summary_dict['total_reverse_inserts'] + global_summary_dict['total_forward_inserts'])
    g_identity = 0 #TODO 1.0 * g_n_matches / (g_n_matches + g_n_mismatches + g_n_deletes + g_n_inserts)
    global_summary_dict['total_identity'] = g_identity
    global_summary_dict['forward_coverage_estimate'] = -1
    global_summary_dict['reverse_coverage_estimate'] = -1

    global_summary_data = ['genome', 'genome']
    for k in summary_headers:
        global_summary_data.append(global_summary_dict[k])
    summary_data.append(global_summary_data)

    return summary_headers, summary_data, identities, identities_per_file, read_lengths_per_file, read_len_to_identity


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


def plot_identity_histogram(identities, output_location=None, title=None):
    step = 0.0025        # bin size
    max_length = 1      # end of histogram

    bins = numpy.arange(0, max_length + step, step=step)
    frequencies, _ = numpy.histogram(identities, bins=bins)

    plot_distribution(step=step, bins=bins, frequencies=frequencies, save_location=output_location, title=title)


def plot_read_len_to_identity(read_len_to_identity, output_base=None, title=None):
    lengths = list(map(lambda x: x[0], read_len_to_identity))
    log_lengths = list(map(numpy.log2, lengths))
    identities = list(map(lambda x: x[1], read_len_to_identity))

    df = pd.DataFrame({"Length":lengths, "Identity": identities, "Log2 Length": log_lengths})

    # standard len
    sns.jointplot(x="Length", y="Identity", data=df, alpha=.005, s=.5, xlim=[-10000,155000], ylim=[.6,1],
                  marginal_kws=dict(bins=400))
    if output_base is not None:
        pyplot.savefig("{}.dot.png".format(output_base))
    pyplot.show()
    pyplot.close()

    # log len
    sns.jointplot(x="Log2 Length", y="Identity", data=df, alpha=.01, s=.5)
    if output_base is not None:
        pyplot.savefig("{}.loglen-dot.png".format(output_base))
    pyplot.show()
    pyplot.close()

    # downsampling (kde is expeeeeensive)
    lengths = lengths[:100000]
    identities = identities[:100000]
    ds_len = list()
    ds_loglen = list()
    ds_iden = list()
    ds_likelihood = 50000.0 / len(identities)
    for l, ll, i in zip(lengths, log_lengths, identities):
        if random.random() < ds_likelihood:
            ds_len.append(l)
            ds_loglen.append(ll)
            ds_iden.append(i)
    print("Downsampled for kde: {}, {}, {}".format(len(ds_len), len(ds_loglen), len(ds_iden)))

    # kde
    df = pd.DataFrame({"Length":ds_len, "Log2 Length": ds_loglen, "Identity": ds_iden})
    sns.jointplot(x="Length", y="Identity", data=df, kind="kde", xlim=[-10000,155000], ylim=[.6,1])
    if output_base is not None:
        pyplot.savefig("{}.kde.png".format(output_base))
    pyplot.show()
    pyplot.close()

    # kde loglen
    sns.jointplot(x="Log2 Length", y="Identity", data=df, kind="kde")
    if output_base is not None:
        pyplot.savefig("{}.loglen-kde.png".format(output_base))
    pyplot.show()
    pyplot.close()


def plot_per_file_identity_violin(raw_identities_per_file, output_base=None, title=None):
    # step = 0.0025        # bin size
    # max_length = 1       # end of histogram
    #
    # frequencies = dict()
    # bins = numpy.arange(0, max_length + step, step=step)
    # for key in identities_per_file.keys():
    #     f, _ = numpy.histogram(identities_per_file[key], bins=bins)
    #     frequencies[key] = f
    #
    # plot_distribution(step=step, bins=bins, frequencies=frequencies, save_location=output_location, title=title)

    identities_per_file = merge_dicts_by_key_idx(raw_identities_per_file, key_start_pos=8, key_end_pos=15, basename=True)

    # sort keys and values together
    sorted_keys, sorted_vals = zip(*sorted(identities_per_file.items(), key=op.itemgetter(0)))

    samples = ["03492", "03098", "02723", "02080", "02055", "01243", "01109", "00733", "24385", "24149", "24143", 'CHM13']

    colors = [(175 / 256.0, 48 / 256.0, 51 / 256.0),  # red
              (224 / 256.0, 99 / 256.0, 58 / 256.0),  # orange
              (215 / 256.0, 219 / 256.0, 84 / 256.0),  # yellow
              (110 / 256.0, 170 / 256.0, 100 / 256.0),  # light green
              (80 / 256.0, 180 / 256.0, 150 / 256.0),  # green
              (100 / 256.0, 189 / 256.0, 197 / 256.0),  # green-blue
              (0 / 256.0, 170 / 256.0, 231 / 256.0),  # turquoise
              (51 / 256.0, 87 / 256.0, 182 / 256.0),  # blue
              (37 / 256.0, 36 / 256.0, 93 / 256.0),  # indigo
              (95 / 256.0, 51 / 256.0, 139 / 256.0),  # purple
              (200 / 256.0, 53 / 256.0, 93 / 256.0),  # pink
    ]

    # print mmm
    for key in sorted_keys:
        mmm(identities_per_file[key], identifier=key)

    # sns.set_palette(sns.husl_palette(len(sorted_keys)))
    sns.set(style='white')
    ax = sns.violinplot(data=sorted_vals, inner=None, linewidth=0, cut=0, palette=sns.color_palette("husl", 8))
    pyplot.setp(ax.collections, alpha=.8)

    # category labels
    pyplot.xticks(pyplot.xticks()[0], sorted_keys)
    pyplot.ylabel("Identities")
    pyplot.xlabel("Samples")
    pyplot.ylim(0.4, 1.05)

    if output_base is not None:
        pyplot.savefig("{}.identity_violin.png".format(output_base))
    pyplot.show()
    pyplot.close()


def plot_identity_comparison_violin(left_identities_per_file, right_identities_per_file, left_lengths_per_file,
                                    right_lengths_per_file, title=None, output_base=None):
    left_key = "Flipflop"
    right_key = "Non-Flipflop"
    # left_key = "No Filter"
    # right_key = "GTE Q7"
    name_key_start=8
    name_key_end=19

    left_identities_per_file = merge_dicts_by_key_idx(left_identities_per_file, key_start_pos=name_key_start,
                                                      key_end_pos=name_key_end, basename=True)
    right_identities_per_file = merge_dicts_by_key_idx(right_identities_per_file, key_start_pos=name_key_start,
                                                       key_end_pos=name_key_end, basename=True)

    l_keys = set(left_identities_per_file.keys())
    r_keys = set(right_identities_per_file.keys())
    keys = l_keys.intersection(r_keys)
    if len(l_keys) != len(r_keys) or len(l_keys) != len(keys) or not any([k in keys for k in l_keys]):
        print("Divergent keys for identity comparison violin")

    # how to order these
    ordered_keys = list(keys)
    ordered_keys.sort()
    ordered_hues = [left_key, right_key]


    # print mmm
    for key in ordered_keys:
        mmm(left_identities_per_file[key], identifier="{} {}".format(left_key, key))
        mmm(right_identities_per_file[key], identifier="{} {}".format(right_key, key))

    # make datafarame
    hands = list()
    idents = list()
    samples = list()

    # for scaling the plot by total length
    left_lengths_per_file = merge_dicts_by_key_idx(left_lengths_per_file, key_start_pos=name_key_start,
                                                   key_end_pos=name_key_end, basename=True)
    right_lengths_per_file = merge_dicts_by_key_idx(right_lengths_per_file, key_start_pos=name_key_start,
                                                    key_end_pos=name_key_end, basename=True)
    lenght_scale_per_file = dict()

    for i, k in enumerate(ordered_keys):
        llpf = sum(map(int, left_lengths_per_file[k]))
        rlpf = sum(map(int, right_lengths_per_file[k]))
        print("{}:\t {} {}, {} {}".format(k, ordered_hues[0], llpf, ordered_hues[1], rlpf))
        lenght_scale_per_file[i] = { #this order must match ordered_hues
            0: llpf / max(llpf, rlpf),
            1: rlpf / max(llpf, rlpf)
        }

    for k in keys:
        for v in left_identities_per_file[k]:
            hands.append(left_key)
            idents.append(v)
            samples.append(k)
        for v in right_identities_per_file[k]:
            hands.append(right_key)
            idents.append(v)
            samples.append(k)
    df = pd.DataFrame({"Sample": samples, "Identity": idents, "Filter Type": hands})

    sns.violinplot(x="Sample", y="Identity", hue="Filter Type", order=ordered_keys, hue_order=ordered_hues,
                   data=df, split=True, inner='quartile',
                   scale='width')
                   # scale='custom', custom_hue_scale=lenght_scale_per_file)
    pyplot.legend(loc=4)

    if output_base is not None:
        pyplot.savefig("{}.identity_comparison_violin.png".format(output_base))
    pyplot.show()
    pyplot.close()


def plot_per_file_identity_curve(identities_per_file, title=None, output_base=None):
    name_key_start=8
    name_key_end=19
    basename = True

    identities_per_file = merge_dicts_by_key_idx(identities_per_file, key_start_pos=name_key_start,
                                                 key_end_pos=name_key_end, basename=basename)

    # how to order these
    ordered_keys = list(identities_per_file.keys())
    ordered_keys.sort()

    # print mmm
    for key in ordered_keys:
        mmm(identities_per_file[key], identifier=key)

    sns.set(style="white", palette="bright", color_codes=True)
    for key in ordered_keys:
        sns.distplot(identities_per_file[key], hist=False, rug=False, kde_kws={"shade": True, "linewidth": 1},
                     color='grey' if 'n' in key else 'blue')
    sns.despine(left=True, top=True, right=True)
    pyplot.yticks([])
    pyplot.xlim([.5, 1])


    if output_base is not None:
        pyplot.savefig("{}.identity_comparison_curve.png".format(output_base))
    pyplot.show()
    pyplot.close()


def mmm(value_list, identifier, key_fcn=None, print_it=True):

    if key_fcn is not None:
        value_list = list(map(key_fcn, list))

    mean = numpy.mean(value_list)
    median = numpy.median(value_list)
    bucketed_value_list = list(map(lambda x: int(x * 200) / 200.0, value_list))
    mode = stats.mode(bucketed_value_list)

    if print_it:
        print("{}: \tmean: {:.5f}\tmedian: {:.5f}\tmode: {:.5f}".format(identifier, mean, median, mode[0][0]))

    return mean, median, mode


def merge_dicts_by_key_idx(identities_per_file, key_start_pos=8, key_end_pos=15, basename=True):
    merged_identities = defaultdict(list)
    for key in identities_per_file.keys():
        new_key = (lambda x : os.path.basename(x) if basename else x)(key)[key_start_pos:key_end_pos]
        # new_key = (lambda x : os.path.basename(x) if basename else x)(key).replace("_non_flipflop","")[key_start_pos:key_end_pos]
        merged_identities[new_key].extend(identities_per_file[key])
    return merged_identities


def main(summary_glob, output_dir, filter_decoys, args):
    FileManager.ensure_directory_exists(output_dir)

    summary_file_paths = glob.glob(summary_glob)
    if len(summary_file_paths) == 0:
        print("No files matched '{}'".format(summary_glob))
        sys.exit(1)

    if filter_decoys:
        print("Filtering decoy chromosomes")
        summary_file_paths = filter_decoys_from_paths(summary_file_paths)

    summary_headers, summary_data, identities, identities_per_file, read_lengths_per_file, read_len_to_identity = \
        aggregate_summary_data(summary_file_paths, args)

    # all_read_lengths = list()
    # for rli in read_len_to_identity:
    #     all_read_lengths.append(rli[0])
    # all_read_lengths.sort()
    # print("top 15 read lengths: {}".format(all_read_lengths[:-15]))


    for file in identities_per_file.keys():
        mmm(identities_per_file[file], file)
    mmm(identities, "All Data")

    sample_name = args.sample
    if sample_name is None:
        sample_name = summary_glob.rstrip('/').replace('/', "_").replace('*', "_")  # replace this with sample name extractor function?

    # save to file
    # write_chromosomal_summary_data_to_csv(summary_headers=summary_headers,
    #                                       summary_data=summary_data,
    #                                       sample_name=sample_name,
    #                                       output_dir=output_dir)

    # plots
    if args.plot:
        pass

        # plot_identity_histogram(identities, title=sample_name, output_location=os.path.join(output_dir, "{}.all_identities.png".format(sample_name)))
        # plot_read_len_to_identity(read_len_to_identity, title=sample_name, output_base=os.path.join(output_dir, "{}.read_len_to_identity".format(sample_name)))
        plot_per_file_identity_violin(identities_per_file, title=sample_name,
                                      output_base=os.path.join(output_dir, sample_name))
        # plot_per_file_identity_curve(identities_per_file, output_base=os.path.join(output_dir, sample_name))
        # if args.comparison_glob is not None:
        #     comparison_paths = glob.glob(args.comparison_glob)
        #     if len(comparison_paths) == 0:
        #         raise Exception("No comparison files found for '{}'".format(args.comparison_glob))
        #     _, _, _, comparison_identities_per_file, comparison_lengths_per_file, _ = aggregate_summary_data(comparison_paths, args)
        #     plot_identity_comparison_violin(identities_per_file, comparison_identities_per_file,
        #                                     read_lengths_per_file, comparison_lengths_per_file,
        #                                     title=sample_name, output_base=os.path.join(output_dir, sample_name))



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
        "--summary_glob", "-i",
        type=str,
        required=True,
        help="GLOB matching input files"
    )
    parser.add_argument(
        "--output_dir", "-o",
        type=str,
        required=False,
        default='.',
        help="desired output directory path (will be created during run time if doesn't exist)"
    )
    parser.add_argument(
        "--filter_decoys", "-f",
        type=stringAsBool,
        required=False,
        default=False,
        help="filter decoy contigs"
    )
    parser.add_argument(
        "--sample", "-s",
        type=str,
        required=False,
        default=None,
        help="sample name (for use in output naming)"
    )
    parser.add_argument(
        "--min_read_length", "-l",
        type=int,
        required=False,
        default=0,
        help="discard reads with length below this threshold"
    )
    parser.add_argument(
        "--plot", "-p",
        type=stringAsBool,
        required=False,
        default=True,
        help="produce plots (default TRUE)"
    )

    parser.add_argument(
        "--comparison_glob", "-c",
        type=str,
        required=False,
        default=None,
        help="GLOB matching all files to compare to"
    )

    args = parser.parse_args()

    main(summary_glob=args.summary_glob, output_dir=args.output_dir, filter_decoys=args.filter_decoys, args=args)

