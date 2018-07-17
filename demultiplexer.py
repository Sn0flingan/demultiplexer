# demultiplexer.py
#
# Runs on Python3
#
# Created by Alice (Sn0flingan) on 2018-07-16
#

import argparse
from Levenshtein import distance

def main():
    args = get_arguments()
    primer_f = "TTGATTACGTCCCTGCCCTTT"
    primer_r = "CCTTAGTAACGGCGAGTGAAA" #reverse compliment of reverse primer
    barcodes = ['AACCACTGGATGGAAA',
                'AAGTAGGGGTCAGCTC',
                'AATCGCATCAAGCGGG',
                'ACCCACATGATATTCC',
                rev_comp('AACCACTGGATGGAAA'),
                rev_comp('AAGTAGGGGTCAGCTC'),
                rev_comp('AATCGCATCAAGCGGG'),
                rev_comp('ACCCACATGATATTCC')] #Last four rev comp for lagging strand
    Barcode_matches = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 100:0}
    with open(args.input) as file:
        read_next_line = False
        for line in file:
            if line[0]=='@':
                read_next_line = True
            elif read_next_line:
                read_start = line[:150]
                read_end = line[-150:-1] #-1 due to newline character as last character of read
                #Forward strand
                print("--- Read ---")
                print("Sequence len: {}".format(len(line)))
                print("-- Start of read")
                barcode_idx = 100
                (start_pos, end_pos, primer_idx) = get_primer_pos(read_start, [primer_f, primer_r], args.verbosity)
                if start_pos is not None and (start_pos-17)>=0:
                    cand_barcode = read_start[start_pos-17:start_pos+5]
                    (start_pos, end_pos, barcode_idx) = get_primer_pos(cand_barcode, barcodes, 2)
                else:
                    cand_barcode = None
                print("Barcode idx {}".format(barcode_idx))
                Barcode_matches[barcode_idx] += 1
                print("-- End of read")
                barcode_idx = 100
                (start_pos, end_pos, primer_idx) = get_primer_pos(read_end, [primer_r, rev_comp(primer_f)], args.verbosity)
                if end_pos is not None:
                    cand_barcode = read_end[end_pos-5:end_pos+17]
                    (start_pos, end_pos, barcode_idx) = get_primer_pos(cand_barcode, barcodes, 2)
                else:
                    cand_barcode = None
                print("Barcode idx {}".format(barcode_idx))
                Barcode_matches[barcode_idx] += 1
                read_next_line = False
    print("--- Demultiplex summary ---")
    print(Barcode_matches)


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input read files in fastq format")
    parser.add_argument("-v", "--verbosity", help="Level of output while running (0 - no output)",
                        default=1, type=int)
    args = parser.parse_args()
    if args.verbosity >=1:
        print("\n--- Input arguments parsed ---")
        print("Input file(s): {}".format(args.input))
        print("Verbosity: {}".format(args.verbosity))
    return args


def get_primer_pos(seq, primers, verbosity):
    min_dist = len(primers[0])
    match_idx = (None, None, 100)
    for i in range(len(primers)):
        (dist, idx) = match_primer(seq, primers[i], verbosity)
        if dist<=7 and dist<min_dist:
            match_idx = (idx, idx+len(primers[i]), i+1)
            min_dist=dist
    return match_idx

def match_primer(sequence, primer, verbosity):
    min_dist = len(sequence)
    best_match_idx = 0
    for i in range(0,len(sequence)-len(primer)):
        seq_dist = distance(sequence[i:i+len(primer)],primer)
        if seq_dist<min_dist:
            best_match_idx = i
            min_dist = seq_dist
    coloured_match = '\x1b[6;31;48m' + sequence[best_match_idx:best_match_idx+len(primer)] + '\x1b[0m'
    if verbosity>=2:
        print("Distance : {}".format(min_dist))
        print(sequence[:best_match_idx] + coloured_match + sequence[best_match_idx+len(primer):] )
    return min_dist, best_match_idx,

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_comp_seq = "".join(complement.get(base,base) for base in reversed(seq))
    return rev_comp_seq

main()
