import sys
from collections import defaultdict
from handlers.FastqReader import FastqReader
from matplotlib import pyplot
import numpy

'''
Iterate a fastq file and find the read length distribution, as well as cumulative number of reads above thresholds
'''


# READS_PATH = "/home/ryan/data/Nanopore/ecoli/miten/guppy/subsampled/30x/r94_ec_rad2.30x.fastq"
# READS_PATH = "/home/ryan/data/Nanopore/ecoli/miten/guppy/r94_ec_guppy.first50k.fastq"
READS_PATH = "/home/ryan/data/Nanopore/ecoli/miten/guppy/r94_ec_guppy.fastq"
# READS_PATH = "/home/ryan/data/Nanopore/ecoli/flapppie/03_22_19_R941_gEcoli_first_410k.fastq"
# READS_PATH = "/home/ryan/data/Nanopore/ecoli/flapppie/03_22_19_R941_gEcoli_last_410k.fastq"
# READS_PATH = "/home/ryan/Downloads/r94_ec_rad2.30x.fastq"


def plot_length_distribution(step, bins, frequencies):
    axes = pyplot.axes()

    center = (bins[:-1] + bins[1:]) / 2

    axes.bar(center, frequencies, width=step, align="center")

    axes.set_xlabel("Read length (bp)")
    axes.set_ylabel("Frequency")

    pyplot.show()
    pyplot.close()


def print_stats(step, frequencies, n_reads):
    print(READS_PATH.strip().split("/")[-1])
    print("\t\t\t\tn\tproportion")

    for threshold in [10000,20000,30000]:
        index = threshold/step - 1
        index = int(index)

        right_side_sum = numpy.sum(frequencies[index:])

        proportion = right_side_sum/n_reads

        # number of reads greater than "threshold"
        print("reads greater than %d:\t%d\t%.3f" % (threshold, right_side_sum, proportion))

    print("n reads total (all lengths):\t%d"%n_reads)


def main():
    parser = FastqReader()

    n_reads = 0
    lengths = list()
    length_sum = 0

    for i, item in enumerate(parser.iterate_file(path=READS_PATH)):
        n_reads += 1

        header, sequence, quality = item

        # print()
        # print(header)
        # print(sequence[:30])
        # print(quality[:30])

        lengths.append(len(sequence))
        length_sum += len(sequence)

        sys.stdout.write("\r%d"%i)

    print()
    
    # ---- Plotting ----

    step = 500          # bin size
    max_length = 50000  # end of histogram

    bins = numpy.arange(0, max_length + step, step=step)
    frequencies, _ = numpy.histogram(lengths, bins=bins)

    print(bins, frequencies)

    plot_length_distribution(step=step, bins=bins, frequencies=frequencies)

    # ---- Printing ----

    print_stats(step=step, frequencies=frequencies, n_reads=n_reads)
    print("total bp:\t%d" % length_sum)
    print("coverage (E. Coli):\t%f" % (length_sum/(5.4*1000000)))


if __name__ == "__main__":
    main()
