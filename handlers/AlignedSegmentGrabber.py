from handlers.BamHandler import BamHandler
from handlers.FastaHandler import FastaHandler
from collections import defaultdict


MAX_COVERAGE = 50
SEQUENCE_LENGTH_CUTOFF_FACTOR = 3   # exclude aligned segments that exceed window_size by this multiple
DEFAULT_MIN_MAP_QUALITY = 5


class SegmentGrabber:
    def __init__(self, chromosome_name, start_position, end_position, ref_sequence, reads, padding=None, padding_end_offset=1, exclude_loose_ends=True):
        # candidate finder includes end position, so should the reference sequence
        self.chromosome_name = chromosome_name
        self.start_position = start_position
        self.end_position = end_position
        self.window_size = end_position - start_position
        self.ref_sequence = ref_sequence
        self.max_coverage = MAX_COVERAGE
        self.padding = padding
        self.padding_end_offset = padding_end_offset
        self.exclude_loose_ends = exclude_loose_ends

        self.reads = reads

        self.read_start_indices = dict()
        self.read_end_indices = dict()
        self.read_alignment_starts = dict()
        self.read_alignment_ends = dict()

        self.sequences = defaultdict(list)
        self.reversal_status = defaultdict()
        self.qualities = defaultdict(list)
        self.cigars = defaultdict(list)

    def get_read_segments(self):
        for r,read in enumerate(self.reads):
            if read.mapping_quality >= DEFAULT_MIN_MAP_QUALITY and read.is_secondary is False \
                    and read.is_supplementary is False and read.is_unmapped is False and read.is_qcfail is False:

                read.query_name = read.query_name + '_1' if read.is_read1 else read.query_name + '_2'
                self.get_aligned_segment_from_read(read)

            if r == self.max_coverage:
                break

        return self.sequences, self.reversal_status

    def get_aligned_segment_from_read(self, read):
        """
        This method finds candidates given a read. We walk through the cigar string to find these candidates.
        :param read: Read from which we need to find the variant candidate positions.
        :return:
        Read candidates use a set data structure to find all positions in the read that has a possible variant.
        """

        read_alignment_start = read.reference_start
        # read_alignment_stop = self.get_read_stop_position(read)

        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence
        read_id = read.query_name
        # read_quality = read.query_qualities

        # read_index: index of read sequence
        # ref_index: index of reference sequence
        read_index = 0
        ref_index = 0
        found_valid_cigar = False
        completion_status = False

        if read_id in self.read_start_indices:
            print("WARNING: read_id hash conflict", read_id)

        for c,cigar in enumerate(cigar_tuples):
            cigar_code = cigar[0]
            length = cigar[1]

            # get the sequence segments that are effected by this operation
            # read_quality_segment = read_quality[read_index:read_index + length]
            read_sequence_segment = read_sequence[read_index:read_index + length]

            # skip parsing the first segment if it is not a match
            if cigar_code != 0 and found_valid_cigar is False:
                # only increment the read index if the non-match cigar code is INS or SOFTCLIP
                if cigar_code == 1 or cigar_code == 4:
                    read_index += length
                continue
            found_valid_cigar = True

            # send the cigar tuple to get attributes we got by this operation
            ref_index_increment, read_index_increment, completion_status = \
                self.parse_cigar_tuple(read_index=read_index,
                                       cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=read_alignment_start + ref_index,
                                       read_sequence=read_sequence_segment,
                                       read_id=read_id,
                                       completion_status=completion_status)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment

            if completion_status or c == len(cigar_tuples) - 1:
                start_index = self.read_start_indices[read_id]
                end_index = self.read_end_indices[read_id]

                segment_alignment_start = self.read_alignment_starts[read_id]
                segment_alignment_end = self.read_alignment_ends[read_id]

                # to simulate Paolo Carnevali's data, all reads should span the full region, match on start and end pos.
                if self.exclude_loose_ends:
                    if segment_alignment_start == self.start_position and segment_alignment_end == self.end_position:
                        if self.padding is not None and self.padding_end_offset is not None:
                            # print(start_index - self.padding, end_index + self.padding + self.padding_end_offset)
                            # print(start_index, end_index)
                            # print(read_sequence[start_index - self.padding:end_index + self.padding + self.padding_end_offset + 1])
                            # print(self.padding*" "+read_sequence[start_index:end_index + 1])

                            start_index = start_index - self.padding
                            end_index = end_index + self.padding + self.padding_end_offset

                        sequence = read_sequence[start_index:end_index + 1]

                        if len(sequence) < SEQUENCE_LENGTH_CUTOFF_FACTOR*self.window_size:
                            self.sequences[read_id] = sequence

                else:
                    if self.padding is not None and self.padding_end_offset is not None:

                        start_index = start_index - self.padding
                        end_index = end_index + self.padding + self.padding_end_offset

                    sequence = read_sequence[start_index:end_index + 1]

                    if len(sequence) < SEQUENCE_LENGTH_CUTOFF_FACTOR * self.window_size:
                        self.sequences[read_id] = sequence
                    else:
                        print("excessive read length found for region", len(sequence), self.window_size)

                # if the read segment has been obtained then fetch its directionality (Forward/Reverse), True if Reverse
                self.reversal_status[read_id] = read.is_reverse

                # else:
                #     print("incomplete read segment")
                #     print("expected interval:", self.start_position, self.end_position)
                #     print("segment interval:", segment_alignment_start, segment_alignment_end)
                # if len(sequence) == 0:
                #     print()
                #     print("***WARNING***: EMPTY SEQUENCE!")
                #     print(read_id)
                #     # print(cigar_tuples)
                #     print("start i\t", start_index)
                #     print("end i\t", end_index)
                #     print("start pos\t\t", read.reference_start)
                #     print("len\t\t\t\t", len(read_sequence))
                #     print("start + length\t", read.reference_start + len(read_sequence))
                #     print(sequence)
                #     # print(read_sequence)
                #     # print(''.join([str(c[0])*c[1] for c in cigar_tuples]))
                #     print()
                # else:
                #     print()
                #     print("GOOD SEQUENCE!")
                #     print(read_id)
                #     # print(cigar_tuples)
                #     print("start i\t",start_index)
                #     print("end i\t", end_index)
                #     print("start pos\t\t", read.reference_start)
                #     print("len\t\t\t\t", len(read_sequence))
                #     print("start + length\t", read.reference_start + len(read_sequence))
                #     print(sequence)
                #     # print(read_sequence)
                #     # print(''.join([str(c[0])*c[1] for c in cigar_tuples]))
                #     print()

                break

        return True

    def parse_match(self, read_index, read_id, alignment_position, length, read_sequence):
        index = read_index
        read_complete = False
        start = alignment_position
        stop = start + length

        for i in range(start, stop):
            # update start position
            if i >= self.start_position and read_id not in self.read_start_indices:
                self.read_start_indices[read_id] = index
                self.read_alignment_starts[read_id] = i

            # update end position
            if i <= self.end_position:
                self.read_end_indices[read_id] = index
                self.read_alignment_ends[read_id] = i
            else:
                read_complete = True
                break

            index += 1

        return read_complete

    def parse_cigar_tuple(self, read_index, cigar_code, length, alignment_position, read_sequence, read_id, completion_status):
        """
        Parse through a cigar operation to find possible candidate variant positions in the read
        :param cigar_code: Cigar operation code
        :param length: Length of the operation
        :param alignment_position: Alignment position corresponding to the reference
        :param read_sequence: Read sequence
        :return:
        cigar key map based on operation.
        details: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        0: "MATCH",
        1: "INSERT",
        2: "DELETE",
        3: "REFSKIP",
        4: "SOFTCLIP",
        5: "HARDCLIP",
        6: "PAD"
        """
        ref_index_increment = length
        read_index_increment = length

        # deal different kinds of operations
        if cigar_code == 0:
            # match
            completion_status = self.parse_match(read_index=read_index,
                                                 read_id=read_id,
                                                 alignment_position=alignment_position,
                                                 length=length,
                                                 read_sequence=read_sequence)

        elif cigar_code == 1:
            # insert
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            ref_index_increment = 0

        elif cigar_code == 2 or cigar_code == 3:
            # delete or ref_skip
            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            read_index_increment = 0

        elif cigar_code == 4:
            # soft clip
            ref_index_increment = 0
            # print("CIGAR CODE ERROR SC")

        elif cigar_code == 5:
            # hard clip
            ref_index_increment = 0
            read_index_increment = 0
            # print("CIGAR CODE ERROR HC")

        elif cigar_code == 6:
            # pad
            ref_index_increment = 0
            read_index_increment = 0
            # print("CIGAR CODE ERROR PAD")

        else:
            raise ("INVALID CIGAR CODE: %s" % cigar_code)

        return ref_index_increment, read_index_increment, completion_status