import re
import sys

def process_fastq(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            header = infile.readline().strip()
            if not header:
                break
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            # Define the pattern to match the sequences to be removed
            pattern = re.compile(r'^(?:none|A|GT|TCA)?([NATCG]{9})GTGA([NATCG]{9})GACA(.*)')
            match = pattern.search(sequence)

            if match:
                # Extract the remaining part of the sequence after the matched pattern
                sequence = match.group(1) + match.group(2) + match.group(3)

                # Calculate the number of characters to remove from the quality string
                removed_length = len(match.group(0)) - len(sequence)
                quality = quality[:len(sequence)] + quality[len(sequence) + removed_length:]

            # Write the processed lines to the output file
            outfile.write("{}\n".format(header))
            outfile.write("{}\n".format(sequence))
            outfile.write("{}\n".format(plus))
            outfile.write("{}\n".format(quality))

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_fastq(input_file, output_file)