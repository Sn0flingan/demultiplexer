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
    primer_r = "TTTCACTCGCCGTTACTAAGG" #both 5' to 3' sequences 'as ordered'
    barcodes = get_barcodes("barcodes.csv")
    for bc in barcodes:
        print(bc)
    return
    barcodes = ['AACCACTGGATGGAAA',
                'AAGTAGGGGTCAGCTC',
                'AATCGCATCAAGCGGG',
                'ACCCACATGATATTCC',
                rev_comp('AACCACTGGATGGAAA'),
                rev_comp('AAGTAGGGGTCAGCTC'),
                rev_comp('AATCGCATCAAGCGGG'),
                rev_comp('ACCCACATGATATTCC')] #Last four rev comp for lagging strand
    barcode_matches = {}
    with open(args.input) as file:
        read_next_line = False
        for line in file:
            if line[0]=='@':
                read_next_line = True
            elif read_next_line:
                read_start = line[:150]
                read_end = line[-150:-1] #-1 due to newline character as last character of read
                #Forward strand
                #print("--- Read ---")
                #print("Sequence len: {}".format(len(line)))

                #print("-- Start of read")
                barcode_idx_s = 100
                (start_pos, end_pos, primer_idx) = get_primer_pos(read_start, [primer_f, rev_comp(primer_r)], 9, args.verbosity)
                if start_pos is not None and (start_pos-21)>=0:
                    cand_barcode = read_start[start_pos-21:start_pos+5]
                    (start_pos, end_pos, barcode_idx_s) = get_primer_pos(cand_barcode, barcodes, 4, args.verbosity)
                else:
                    cand_barcode = None
                #print("Barcode idx {}".format(barcode_idx_s))

                #print("-- End of read")
                barcode_idx_e = 100
                (start_pos, end_pos, primer_idx) = get_primer_pos(read_end, [primer_r, rev_comp(primer_f)], 9, args.verbosity)
                if end_pos is not None:
                    cand_barcode = read_end[end_pos-5:end_pos+21]
                    (start_pos, end_pos, barcode_idx_e) = get_primer_pos(cand_barcode, barcodes, 4, args.verbosity)
                else:
                    cand_barcode = None
                #print("Barcode idx {}".format(barcode_idx_e))

                barcode_name = None
                if barcode_idx_s==100 and barcode_idx_e==100:
                    barcode_name = 'None'
                elif barcode_idx_s==100:
                    if barcode_idx_e>4:
                        barcode_name ='BC_s' + str(barcode_idx_e-4)
                    else:
                        barcode_name = 'BC_e' + str(barcode_idx_e)
                elif barcode_idx_e==100:
                    if barcode_idx_s>4:
                        barcode_name = 'BC_e' + str(barcode_idx_s-4)
                    else:
                        barcode_name = 'BC_s' + str(barcode_idx_s)
                else:
                    if barcode_idx_s>4 and barcode_idx_e>4:
                        barcode_name = 'BC_' + str(barcode_idx_e-4) + '-' + 'BC_' + str(barcode_idx_s-4)
                    elif barcode_idx_s<=4 and barcode_idx_e<=4:
                        barcode_name = 'BC_' + str(barcode_idx_s) + '-' + 'BC_' + str(barcode_idx_e)
                    else:
                        barcode_name = 'Conflicting'
                    
                if barcode_name in barcode_matches:
                    barcode_matches[barcode_name] +=1
                else:
                    barcode_matches[barcode_name] = 1


                read_next_line = False
    print("--- Demultiplex summary ---")
    for b, c in barcode_matches.items():
        print("{}: {}".format(b, c))


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


def get_primer_pos(seq, primers, dist_thresh, verbosity):
    min_dist = len(primers[0])
    match_idx = (None, None, 100)
    for i in range(len(primers)):
        (dist, idx) = match_primer(seq, primers[i], verbosity)
        if dist<=dist_thresh and dist<min_dist:
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

def get_barcodes(file):
    barcodes = []
    with open(file) as file:
        for row in file:
            columns = row.split(';')
            if len(columns)==2:
                barcode = Barcode(name=columns[0], start_name=columns[0],
                                  start_seq=columns[1].rstrip())
            else:
                barcode = Barcode(name=columns[0],
                                  start_name=columns[1], start_seq=columns[2],
                                  end_name=columns[3], end_seq=columns[4].rstrip())
            #print("Name: {}".format(columns[0]))
            barcodes.append(barcode)
    return barcodes

class Barcode():
    def __init__(self, name, start_name="", start_seq="", end_name="", end_seq=""):
        self.name = name
        self.start_name = start_name
        self.start_seq = start_seq.upper()
        if not end_name:
            self.isDual = False
            self.end_name = start_name
            self.end_seq = start_seq.upper()
        else:
            self.isDual = True
            self.end_name = end_name
            self.end_seq = end_seq.upper()

    def __str__(self):
        n = "Name: {}\n".format(self.name)
        sn = "Start name: {}\n".format(self.start_name)
        ss = "Start seq: {}\n".format(self.start_seq)
        en = "End name: {}\n".format(self.end_name)
        es = "End seq: {}\n".format(self.end_seq)
        return n + sn + ss + en + es

main()
