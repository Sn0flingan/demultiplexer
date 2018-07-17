# demultiplexer.py
#
# Runs on Python3
#
# Created by Alice (Sn0flingan) on 2018-07-16
#

import argparse
from Levenshtein import distance
import seaborn as sns

def main():
    args = get_arguments()
    primer_f = "TTGATTACGTCCCTGCCCTTT"
    primer_r = "CCTTAGTAACGGCGAGTGAAA" #reverse compliment of reverse primer
    matches = {'Leading':{'Both':0, 'Start':0, 'End':0}, 'Lagging': {'Both':0, 'Start':0, 'End':0}, 'None': 0}
    index = {'Matching-start': [], 'Matching-end': [] }
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
                print("-- Leading strand")
                (f_dist_s, f_idx_s) = check_barcode(read_start, primer_f)
                (f_dist_e, f_idx_e) = check_barcode(read_end, primer_r)
                print("-- Lagging strand")
                (r_dist_s, r_idx_s) = check_barcode(read_start, rev_comp(primer_r))
                (r_dist_e, r_idx_e) = check_barcode(read_end, rev_comp(primer_f))
                if f_dist_s<=7 and f_dist_e<=7:
                    matches['Leading']['Both'] +=1
                    index['Matching-start'].append(f_idx_s)
                    index['Matching-end'].append(f_idx_e)
                elif r_dist_s<=7 and r_dist_e<=7:
                    matches['Lagging']['Both'] +=1
                    index['Matching-start'].append(r_idx_s)
                    index['Matching-end'].append(r_idx_e)
                elif f_dist_s<=7 and f_dist_e>7:
                    matches['Leading']['Start'] +=1
                    index['Matching-start'].append(f_idx_s)
                elif f_dist_s>7 and f_dist_e<=7:
                    matches['Leading']['End'] +=1
                    index['Matching-end'].append(f_idx_e)
                elif r_dist_s<=7 and r_dist_e>7:
                    matches['Lagging']['Start'] +=1
                    index['Matching-start'].append(r_idx_s)
                elif r_dist_s>7 and r_dist_e<=7:
                    matches['Lagging']['End'] +=1
                    index['Matching-end'].append(r_idx_e)
                else:
                    matches['None'] +=1
                read_next_line = False
    print(matches)
    histogram_plot(index['Matching-start'], "Start")
    histogram_plot(index['Matching-end'], "End")
    #print(sorted(index['Matching']))



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
    return (min_dist, best_match_idx)

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_comp_seq = "".join(complement.get(base,base) for base in reversed(seq))
    return rev_comp_seq

def histogram_plot(list, name):
    sns.set(style="whitegrid")
    plot = sns.distplot(list, kde=False, color="b")
    fig = plot.get_figure()
    fig.savefig(name + "_histogram.png")
main()
