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
                print("-- Forward strand")
                f_dist_s = check_barcode(read_start, primer_f)
                if f_dist_s>=5:
                    print("- Reverse-comp")
                    check_barcode(read_start, rev_comp(primer_f))
                f_dist_e = check_barcode(read_end, primer_r)
                if f_dist_e>=5:
                    print("- Reverse-comp")
                    check_barcode(read_end, rev_comp(primer_r))
                print("-- Reverse strand")
                r_dist_s = check_barcode(read_start, rev_comp(primer_r))
                if r_dist_s>=5:
                    print("- Reverse-comp")
                    check_barcode(read_start, primer_r)
                r_dist_e = check_barcode(read_end, rev_comp(primer_f))
                if r_dist_e>=5:
                    print("- Reverse-comp")
                    check_barcode(read_end, primer_f)
                read_next_line = False


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

def check_barcode(sequence, primer):
    min_dist = len(sequence)
    best_match_idx = 0
    for i in range(0,len(sequence)-len(primer)):
        seq_dist = distance(sequence[i:i+len(primer)],primer)
        if seq_dist<min_dist:
            best_match_idx = i
            min_dist = seq_dist
    coloured_match = '\x1b[6;31;48m' + sequence[best_match_idx:best_match_idx+len(primer)] + '\x1b[0m'
    #print(sequence[:best_match_idx] + coloured_match + sequence[best_match_idx+len(primer):] )
    print("Distance : {}".format(min_dist))
    return min_dist

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_comp_seq = "".join(complement.get(base,base) for base in reversed(seq))
    return rev_comp_seq

main()
