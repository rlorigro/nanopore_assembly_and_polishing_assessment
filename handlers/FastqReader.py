

class FastqReader:
    def __init__(self):
        self.start = False
        self.index = 0

        self.header = None
        self.sequence = None
        self.quality = None

    def iterate_file(self, path):
        reads_file = open(path, "r")

        for l,line in enumerate(reads_file):

            # skip lines until the first header is found
            if not self.start and line[0] == "@":
                self.start = True

            if self.start:
                # Header
                if self.index % 4 == 0:
                    self.header = self.parse_header(line, l)

                # Sequence
                elif self.index % 4 == 1:
                    self.sequence = self.parse_sequence(line, l)

                # Header B
                elif self.index % 4 == 2:
                    # should be empty or redundant
                    pass

                # Quality
                elif self.index % 4 == 3:
                    self.quality = self.parse_quality(line, l)

                    # Item complete, return read data
                    yield [self.header, self.sequence, self.quality]

            self.index += 1

        reads_file.close()

    def parse_header(self, line, line_index):
        if not line[0] == "@":
            exit("ERROR: non fastq symbol found at line: %d"%line_index)
        return line.strip()[1:]

    def parse_sequence(self, line, line_index):
        return line.strip()

    def parse_quality(self, line, line_index):
        return line.strip()
