from collections import Counter
import math


def calculate_shannon_entropy(sequence):
    counter = Counter(sequence)
    entropy = 0

    for character in counter:
        frequency = counter[character]
        probability = frequency/len(sequence)

        entropy -= probability*math.log(probability, 2)

    return entropy


def find_longest_repeat(sequence):
    prev_character = None
    max_length = 0

    c = None
    for character in sequence:
        if character == prev_character:
            c += 1
        else:
            if c is not None and c > max_length:
                max_length = c
            c = 1

        prev_character = character

    if c > max_length:
        max_length = c

    return max_length


def test_entropy():
    sequence = "AGGTTGGGGGAAAAAA"

    entropy = calculate_shannon_entropy(sequence)
    max_length = find_longest_repeat(sequence)
    print(entropy, max_length)


if __name__ == "__main__":
    test_entropy()
