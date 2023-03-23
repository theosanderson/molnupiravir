print("Loading modules")
import pysam
import argparse
# make argparse, with input as a bam file and output of a TSV file

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input bam file")
parser.add_argument("-o", "--output", help="output tsv file")
args = parser.parse_args()
output = open(args.output, "wt")
bam = pysam.AlignmentFile( args.input, "rb")
print("Loaded bam file", args.input)

def write_line(pos, base_counts):
    # write the position and the counts for each base
    output.write(str(pos))
    for base in "ACGT":
        output.write("\t" + str(base_counts[base]))
    output.write("\n")

# get the header
header = bam.header
# list the number of As, Cs, Gs, Ts for each position
# in the reference

# get the reference sequence
ref = header["SQ"][0]["SN"]
# get the length of the reference sequence
ref_len = header["SQ"][0]["LN"]
# initialize the counts
from collections import defaultdict, Counter
# loop through the reads in a pileup
import tqdm
for pileupcolumn in tqdm.tqdm(bam.pileup(ref)):
    # get the position
    pos = pileupcolumn.pos
    base_counts = Counter()
    # loop through the reads at that position
    for pileupread in pileupcolumn.pileups:
        # get the base
        if pileupread.is_del or pileupread.is_refskip:
            continue
        base = pileupread.alignment.query_sequence[pileupread.query_position]
        # increment the count for that base
        base_counts[base] += 1
    # if total is >50, and second makes up >10% of the total then print
    total = sum(base_counts.values())
    # if total > 50 and len(base_counts) > 1:
    #     second = sorted(base_counts.values())[-2]
    #     if second/total > 0.1:
    #         write_line(pos, base_counts)
    write_line(pos, base_counts)