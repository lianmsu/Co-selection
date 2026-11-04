#!/usr/bin/env python3
"""Convert the long diamond hit table into a genome x gene count table.

Usage:
    python3 turn_to_table.py <input_file>

The last two columns of the input are interpreted as (row name, column name),
i.e. gene and genome, and their co-occurrences are counted.  The output is
written next to the input as <input_basename>_table.csv (tab separated).

Original script: turn_to_table.py
"""
import os
import sys
from collections import defaultdict

import pandas as pd


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    input_filename = os.path.basename(input_file).split(".")[0]

    data = defaultdict(lambda: defaultdict(int))
    with open(input_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) >= 2:
                row_name, col_name = cols[-2], cols[-1]
                data[row_name][col_name] += 1

    df = pd.DataFrame(data).fillna(0).astype(int)
    df = df.T

    output_file = input_filename + "_table.csv"
    df.to_csv(output_file, sep="\t")
    print(f"Output saved to {output_file}")


if __name__ == "__main__":
    main()
