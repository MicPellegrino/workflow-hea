import numpy as np


def clean_lines(filename):
    """Remove empty lines and comments (# ...)"""
    with open(filename, "r") as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.lstrip().startswith("#")
        ]


def read_scalar_field(filename, value_name):
    """
    Generic reader for BOTH files.

    Returns:
        dict:
            frame -> {
                value_name: np.array([...])
            }
    """
    lines = clean_lines(filename)

    data = {}
    i = 0

    while i < len(lines):

        header = lines[i].split()

        frame = int(header[0])
        nbin = int(header[1])

        values = np.zeros(nbin)

        for j in range(nbin):
            parts = lines[i + 1 + j].split()

            bin_id = int(parts[0]) - 1
            values[bin_id] = float(parts[1])

        data[frame] = {
            value_name: values
        }

        i += nbin + 1

    return data


def merge_fields(a_dict, b_dict, key_b="temperature"):
    """
    Merge two frame dictionaries into one.
    """
    for frame in b_dict:
        if frame not in a_dict:
            a_dict[frame] = {}

        a_dict[frame][key_b] = b_dict[frame][key_b]

    return a_dict


if __name__ == "__main__":
    main()