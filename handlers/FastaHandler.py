from pysam import FastaFile
"""
This class handles fasta reference files, ensuring that the sequence is not a terminal 'N' and that the end of the
sequence has not been reached
"""


class FastaHandler:
    """
    Handles fasta files using pyfaidx API
    """
    def __init__(self, reference_file_path):
        """
        create fasta file object given file path to a fasta reference file
        :param fasta_file_path: full path to a fasta reference file
        """

        self.fasta_file_path = reference_file_path

        try:
            self.fasta = FastaFile(self.fasta_file_path)
        except:
            raise IOError("FASTA FILE READ ERROR")

    def get_sequence(self, chromosome_name, start, stop):
        """
        Return the sequence of a query region
        :param chromosome_name: Chromosome name
        :param start: Region start
        :param stop: Region end
        :return: Sequence of the region
        """
        return self.fasta.fetch(region=chromosome_name, start=start, end=stop).upper()

    def get_chr_sequence_length(self, chromosome_name):
        """
        Get sequence length of a chromosome. This is used for selecting windows of parallel processing.
        :param chromosome_name: Chromosome name
        :return: Length of the chromosome reference sequence
        """
        return self.fasta.get_reference_length(chromosome_name)

    def get_contig_names(self):
        return self.fasta.references

    def get_ref_of_region(self, contig, site):
        """
        Return a string containing reference of a site
        :param contig: Contig [ex chr3]
        :param site: Site [ex 100000-200000]
        :return:
        """
        ret_val = ""
        error_val = 0
        try:
            ret_val = self.fasta.fetch(region=contig+site).upper()
        except:
            print("ERROR IN REF FETCH: ", contig, site)
            error_val = 1
        return ret_val, error_val

    def close(self):
        self.fasta.close()
