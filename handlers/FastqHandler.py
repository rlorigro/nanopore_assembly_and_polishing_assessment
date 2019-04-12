import sys


class Read:
    def __init__(self, read_id, sequence, quality):
        self.id = read_id
        self.sequence = ''.join(sequence)
        self.quality = quality


class FastqHandler:
    def __init__(self, path, strip_header=False):
        """
        Handles Fastq sequence files.
        :param path:
        :param strip_header: If True, reads yielded by this handler will have stripped headers, in which
        only the first element is returned
        """
        self.path = path
        self.start = False
        self.index = 0

        self.strip_header = strip_header

        self.header = None
        self.sequence = None
        self.quality = None

    def extract_reads_by_id(self, id_set, output_path, print_status=False, sequence_cutoff=sys.maxsize):
        file = open(self.path, "r")
        out_file = open(output_path, "w")

        valid_read = False
        read_lines = list()
        n_sequences = 0
        n_found = 0
        for l, line in enumerate(file):
            # Header
            if line[0] == "@":
                n_sequences += 1

                if valid_read and len(read_lines) > 0:
                    n_found += 1

                    # Reset flag
                    valid_read = False

                    # Write the buffer
                    for read_line in read_lines:
                        out_file.write(read_line)

                    read_lines = list()

                if print_status and l%1000:
                    sys.stderr.write("\r %d %d" % (n_sequences,n_found))

                read_id = line.strip()[1:].split(" ")[0]

                if read_id in id_set:
                    valid_read = True

                if n_found == sequence_cutoff:
                    break

            if valid_read:
                read_lines.append(line)

        if valid_read and len(read_lines) > 0:
            # Write the buffer
            for read_line in read_lines:
                out_file.write(read_line)

        sys.stderr.write("\n")
        file.close()
        out_file.close()

    def iterate_read_names(self, print_status=False, line_number=False):
        reads_file = open(self.path, "r")
        n_reads = 0

        for l,line in enumerate(reads_file):
            # skip lines until the first header is found
            if not self.start and line[0] == "@":
                self.start = True

            if self.start:
                # Header
                if self.index % 4 == 0:
                    header = self.parse_header(line)
                    n_reads += 1

                    if print_status:
                        sys.stderr.write("\r %d" % n_reads)

                    if line_number:
                        yield l, header

                    else:
                        yield header

                self.index += 1

        sys.stderr.write("\n")
        reads_file.close()
        self.__init__(self.path)

    def iterate_file(self):
        reads_file = open(self.path, "r")

        for l,line in enumerate(reads_file):

            # skip lines until the first header is found
            if not self.start and line[0] == "@":
                self.start = True

            if self.start:
                # Header
                if self.index % 4 == 0:
                    self.header = self.parse_header(line)

                # Sequence
                elif self.index % 4 == 1:
                    self.sequence = self.parse_sequence(line)

                # Header B
                elif self.index % 4 == 2:
                    # should be empty or redundant
                    pass

                # Quality
                elif self.index % 4 == 3:
                    self.quality = self.parse_quality(line)

                    # Item complete, return read data
                    yield Read(read_id=self.header, sequence=self.sequence, quality=self.quality)

                self.index += 1

        reads_file.close()
        self.__init__(self.path)

    def parse_header(self, line):
        line = line.strip()[1:]

        if self.strip_header:
            line = line.split(" ")[0]

        return line

    def parse_sequence(self, line):
        return line.strip()

    def parse_quality(self, line):
        return line.strip()
