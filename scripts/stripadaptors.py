import gzip
import sys

# Ensure the correct usage of command-line arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <input_file> <output_file>")
    sys.exit(1)

input_file, output_file = sys.argv[1], sys.argv[2]

# Define the adapter sequences
adapter_fa = "XXX"
adapter_sa = "XXX"

def process_chunk(chunk, out):
    for i in range(0, len(chunk), 4):
        out.write(chunk[i])
        
        end_of_sequence = len(chunk[i+1].strip().split(adapter_fa)[0])
        out.write(chunk[i+1].strip()[:end_of_sequence] + '\n')
        
        out.write(chunk[i+2])
        
        out.write(chunk[i+3].strip()[:end_of_sequence] + '\n')

def process_file(input_file, output_file):
    chunk_size = 100000  # Number of lines to read at a time (adjust as needed)
    with gzip.open(input_file, 'rt') as data, gzip.open(output_file, 'wt') as out:
        chunk = []
        for line in data:
            chunk.append(line)
            if len(chunk) == chunk_size:
                process_chunk(chunk, out)
                chunk = []
        if chunk:
            process_chunk(chunk, out)

process_file(input_file, output_file)

print("successfully removed adaptors")