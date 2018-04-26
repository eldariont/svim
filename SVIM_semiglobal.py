"""Module implementing functions for sequence alignment using dynamic programming"""

import sys
import numpy as np


def nw_compute_matrix(a, b, costs=(2, -3, -2), backwards=False):
    """Compute alignment matrix for alignment. No penalty for gaps at beginning of sequences."""
    if backwards:
        a = a[::-1]
        b = b[::-1]
    matrix = np.zeros((len(a)+1, len(b)+1), dtype=int)

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i-1] == b[j-1]:
                match = matrix[i-1, j-1] + costs[0]
            else:
                match = matrix[i-1, j-1] + costs[1]
            delete = matrix[i-1, j] + costs[2]
            insert = matrix[i, j-1] + costs[2]
            matrix[i, j] = max(match, insert, delete)
    return matrix


def nw_get_alignment(a, b, matrix, costs=(2, -3, -2), backwards=False):
    """Produce alignment strings from filled alignment matrix."""
    if backwards:
        a = a[::-1]
        b = b[::-1]
    i, j = np.unravel_index(np.argmax(matrix), (len(a)+1, len(b)+1))

    alin_a = a[i:].lower()
    alin_b = b[j:].lower()

    while i > 0 and j > 0:
        if a[i-1] == b[j-1]:
            match_cost = costs[0]
        else:
            match_cost = costs[1]
        if matrix[i, j] == matrix[i-1, j-1] + match_cost:
            alin_a = a[i-1].upper() + alin_a
            alin_b = b[j-1].upper() + alin_b
            i -= 1
            j -= 1
        elif matrix[i, j] == matrix[i-1, j] + costs[2]:
            alin_a = a[i-1].upper() + alin_a
            alin_b = "-" + alin_b
            i -= 1
        elif matrix[i, j] == matrix[i, j-1] + costs[2]:
            alin_a = "-" + alin_a
            alin_b = b[j-1].upper() + alin_b
            j -= 1

    if backwards:
        return alin_a[::-1], alin_b[::-1]
    else:
        return alin_a, alin_b


def print_alignment(alin_a, alin_b, line_length=50, backwards=False):
    """Print alignment strings nicely. Break lines at given line length."""
    if backwards:
        for i in range(((max(len(alin_a), len(alin_b)) + (line_length - 1)) // line_length) * -line_length,
                       -line_length, line_length):
            print(i)
            print("|")
            print(alin_a[i:i+line_length].rjust(line_length))
            print(alin_b[i:i+line_length].rjust(line_length))
            print("")
        print(-50)
        print("|")
        print(alin_a[-line_length:].rjust(line_length))
        print(alin_b[-line_length:].rjust(line_length))
        print("")
    else:
        for i in range(0, max(len(alin_a), len(alin_b)), line_length):
            print(i)
            print("|")
            print(alin_a[i:i+line_length])
            print(alin_b[i:i+line_length])
            print("")


def get_end_of_alignment(matrix, a_len, b_len, backwards=False):
    """Get row and column indices of highest score from filled alignment matrix."""
    i, j = np.unravel_index(np.argmax(matrix), (a_len+1, b_len+1))
    if backwards:
        return a_len - i, b_len - j
    else:
        return i-1, j-1


def main():
    a = raw_input("Sequence A: ")
    b = raw_input("Sequence B: ")

    print("Computing semi-global alignment of two sequences (lengths {0} and {1})..".format(len(a), len(b)))

    costs = (3, -12, -12)
    matrix = nw_compute_matrix(a, b, costs)

    print("Retrieve alignment string..")
    alin_a, alin_b = nw_get_alignment(a, b, matrix, costs)

    print(matrix)

    print_alignment(alin_a, alin_b)


if __name__ == "__main__":
    sys.exit(main())
