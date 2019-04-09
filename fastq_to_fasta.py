from handlers.FastqReader import FastqReader
from handlers.FastaWriter import FastaWriter
import os

HEADER = 0
SEQUENCE = 1
QUALITY = 2


def main():
    fastq_path = "/home/ryan/data/Nanopore/ecoli/miten/guppy/subsampled/11-29/r94_ec_rad2.30x-30kb.fastq"
    output_directory = "output/"
    output_filename_prefix = ".".join(os.path.basename(fastq_path).split(".")[:-1])
    output_filename = output_filename_prefix + ".fasta"

    output_path = os.path.join(output_directory, output_filename)

    print("WRITING:", output_path, "...")

    fastq_reader = FastqReader()
    fasta_writer = FastaWriter(output_file_path=output_path)

    for read in fastq_reader.iterate_file(fastq_path):
        header = read[HEADER].split(" ")[0]
        fasta_writer.write_entry(sequence=read[SEQUENCE], label=header)


if __name__ == "__main__":
    main()