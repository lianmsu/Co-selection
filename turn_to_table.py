#!/home/test/miniconda3/bin/python3
import sys
import pandas as pd
from collections import defaultdict

if len(sys.argv) != 2:
    print("Usage: python script.py <input_file>")
    sys.exit(1)
input_file = sys.argv[1]
input_filename = input_file.split('/')[-1].split('.')[0]
data = defaultdict(lambda: defaultdict(int))
with open(input_file, 'r') as f:
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) >= 2:
            row_name, col_name = cols[-2], cols[-1]
            data[row_name][col_name] += 1  

df = pd.DataFrame(data).fillna(0).astype(int)
df = df.T  
output_file = input_filename + '_table.csv'
df.to_csv(output_file, sep='\t')
print(f"Output saved to {output_file}")
