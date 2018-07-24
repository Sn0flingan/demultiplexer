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
    barcodes = get_barcodes(args.barcodes)
    dist_thresh = 15

    barcode_matches = {}
    with open(args.input) as file:
        read_next_line = False
        for line in file:
            if line[0]=='@':
                read_next_line = True
            elif read_next_line:
                tags = get_tags(line, [primer_f, rev_comp(primer_r)], dist_thresh, args.verbosity)
                
                matching_barcode = match_barcode(tags, barcodes)
                    
                if matching_barcode in barcode_matches:
                    barcode_matches[matching_barcode] +=1
                else:
                    barcode_matches[matching_barcode] = 1
                
                read_next_line = False

    print("--- Demultiplex summary ---")
    for b, c in barcode_matches.items():
        print("{}: {}".format(b, c))


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input read files in fastq format")
    parser.add_argument("-b", "--barcodes", help="Barcode file in .csv format, each line one barcode pair. \
                        Singleplex: Name;Seq   Dualplex: Name;Name-for;Seq-for;Name-rev;Seq-rev")
    parser.add_argument("-v", "--verbosity", help="Level of output while running (0 - no output)",
                        default=1, type=int)
    args = parser.parse_args()
    if args.verbosity >=1:
        print("\n--- Input arguments parsed ---")
        print("Input file(s): {}".format(args.input))
        print("Verbosity: {}".format(args.verbosity))
    return args

def get_tags(read, primers, dist_thresh, verbosity):
    read_start = read[:150]
    read_end = read[-151:-1]
    #leading strand
    (s_dist, s_idx) = match_primer(read_start, primers[0], dist_thresh, verbosity) #start of read
    (e_dist, e_idx) = match_primer(read_end, primers[1], dist_thresh, verbosity) #end of read
    #lagging strand
    (rc_s_dist, rc_s_idx) = match_primer(read_start, rev_comp(primers[1]), dist_thresh, verbosity) #start of read
    (rc_e_dist, rc_e_idx) = match_primer(read_end, rev_comp(primers[0]), dist_thresh, verbosity) #end of read
    
    read_is_lagging = (rc_s_dist + rc_e_dist)<(s_dist + e_dist) #might need more sophisticated method
    tags = ["Primer not found"]*2
    if read_is_lagging:
        if rc_s_idx and (rc_s_idx-21)>=0:
            tags[1] = rev_comp(read_start[rc_s_idx-21:rc_s_idx+5]) #reverse order and rev_comp to adjust lagging to leading
        if rc_e_idx:
            tags[0] = rev_comp(read_end[rc_e_idx-5:rc_e_idx+21])
    else:
        if s_idx and (rc_s_idx-21)>=0:
            tags[0] = read_start[s_idx-21:s_idx+5]
        if e_idx:
            tags[1] = read_end[e_idx-5:e_idx+21]
    return tags

def match_primer(sequence, primer, dist_thresh, verbosity):
    min_dist = len(sequence)
    best_match_idx = 0
    for i in range(0,len(sequence)-len(primer)):
        seq_dist = distance(sequence[i:i+len(primer)],primer)
        if seq_dist<min_dist and seq_dist<dist_thresh:
            best_match_idx = i
            min_dist = seq_dist
    coloured_match = '\x1b[6;31;48m' + sequence[best_match_idx:best_match_idx+len(primer)] + '\x1b[0m'
    if verbosity>=2:
        print("Distance : {}".format(min_dist))
        print(sequence[:best_match_idx] + coloured_match + sequence[best_match_idx+len(primer):] )
    return min_dist, best_match_idx,

def match_barcode(tags, barcodes):
    min_dist = 1000
    matching_barcode = "None"
    dist_thresh = 7
    verbosity = 1

    for i in range(0, len(barcodes)):
        if tags[0]=="Primer not found":
            s_dist = 100
        else:
            (s_dist, s_pos) = match_primer(tags[0], barcodes[i].start_seq, dist_thresh, verbosity)
        if tags[1]=="Primer not found":
            e_dist = 100
        else:
            (e_dist, e_pos) = match_primer(tags[1], barcodes[i].end_seq, dist_thresh, verbosity)

        if (s_dist+e_dist)==200:
            matching_barcode = "Primers not found"
            break;
        elif (s_dist+e_dist)<min_dist:
            min_dist = (s_dist+e_dist)
            matching_barcode = barcodes[i].name
                
    if 52>min_dist>=26:
        matching_barcode = matching_barcode + " single(b)"
    elif min_dist==52:
        matching_barcode = "Barcodes not found"
    elif 200>min_dist>=100:
        matching_barcode = matching_barcode + " single(p)"

    return matching_barcode



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
