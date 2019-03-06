

def get_ordering(chromsome_name, prefix):
    ordering = list()

    suffix = get_chromosome_suffix(chromsome_name, prefix)

    if suffix.isdigit():
        index = int(suffix)
    else:
        index = -1

    ordering.append(index)
    ordering.append(suffix)

    # print(suffix, ordering)
    return ordering


def get_chromosome_suffix(chromosome_name, prefix):
    suffix = chromosome_name.split(prefix)[-1]

    return suffix


def sort_chromosome_names(names, prefix="chr"):
    summary_data = sorted(names, key=lambda x: get_ordering(x, prefix))

    return summary_data


def sort_summary_data(data, chromosome_name_index, prefix="chr"):
    summary_data = sorted(data, key=lambda x: get_ordering(x[chromosome_name_index], prefix))

    return summary_data


def test_sort():
    names = ["chr2", "chr11", "chrY", "chr1",  "chr20", "chrX", "chr10"]

    sorted_names = sort_chromosome_names(names)

    print("unsorted\t", names)
    print("sorted\t\t", sorted_names)

    data = [["ctg75","chr2","SNP","27790350","27790351","A","15","16","G","False","CTTTCAGGCGT","1.868","3"],
            ["ctg75","chr11","SNP","27790354","27790355","G","19","20","C","False","CAGGCGTATGG","1.859","2"],
            ["ctg75","chrY","SNP","27790357","27790358","T","22","23","A","False","GCGTATGGTGA","1.79","2"],
            ["ctg75","chr1","SNP","27790375","27790376","C","40","41","T","False","AATATCTTCCC","1.573","3"],
            ["ctg75","chr20","INS","27790378","27790378","","43","44","C","False","ATCTTCCCGTA","1.823","3"],
            ["ctg75","chrX","SNP","27790381","27790382","G","47","48","A","False","TTCCCGTAAAA","1.868","4"],
            ["ctg75","chr10","SNP","27790382","27790383","T","48","49","A","False","TCCCGTAAAAA","1.79","5"]]

    sorted_data = sort_summary_data(data=data, chromosome_name_index=1)

    for item in data:
        print(item)
    print()

    for item in sorted_data:
        print(item)


if __name__ == "__main__":
    test_sort()
