

class FastaWriter:
    def __init__(self, output_file_path, label_prefix="read"):
        self.file_path = output_file_path
        self.label_prefix = label_prefix + "_"

        open(output_file_path, 'w').close()

    def write_file(self, sequences, labels=None, reversal_statuses=None, append=False):
        file_operation = "w"

        if append:
            file_operation = "a"

        with open(self.file_path, file_operation) as file:
            for i in range(len(sequences)):
                sequence = sequences[i]

                if labels is None:
                    label = None
                else:
                    label = labels[i]

                if reversal_statuses is None:
                    reversal = None
                else:
                    reversal = reversal_statuses[i]

                string = self.generate_sequence_string(sequence=sequence,
                                                       label=label,
                                                       reversal=reversal,
                                                       i=i)

                file.write(string)

            if append:
                file.write("\n")

    def write_entry(self, sequence, label=None, reversal=None, i=None):
        if label is None and i is None:
            exit("If no label is given, index must be provided")

        file_operation = "a"

        with open(self.file_path, file_operation) as file:
            string = self.generate_sequence_string(sequence=sequence,
                                                   label=label,
                                                   reversal=reversal,
                                                   i=i)

            file.write(string)

    def generate_label_string(self, label, reversal, i):
        if label is None:
            label = str(i)

        if reversal is None:
            reversal_flag = ""
        else:
            if reversal:
                reversal_flag = "_R"
            else:
                reversal_flag = "_F"

        label_string = ">" + self.label_prefix + label + reversal_flag

        return label_string

    def generate_sequence_string(self, sequence, label, reversal, i):
        label = self.generate_label_string(label, reversal, i)

        if i == 0:
            entry = [label,sequence]
        else:
            entry = ['',label,sequence]

        string = '\n'.join(entry)

        return string


def test_fasta_writer():
    path = "output/fasta_test.txt"
    sequences = ["ACTG","CTGA","TGAC","GACT"]

    writer = FastaWriter(output_file_path=path)

    writer.write_sequences(sequences)


if __name__ == "__main__":
    test_fasta_writer()
