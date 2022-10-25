#!/usr/bin/env python

#Imports:
import argparse
import re

# Functions:

#Accept command line input for fields mentioned in assignment README.md
def get_args():
#     * -f, --file: input sorted SAM file (absolute path)
#     * -o, --outfile: output sorted, duplicate-free SAM file (absolute path)
#     * -u, --umi: input file w/ list of UMIs (absolute path)
#     * -h, --help: help message for deduper function
    parser = argparse.ArgumentParser(description="A program to remove duplicate reads from a given SAM file.")
    parser.add_argument("-f","--file",help="Input sorted SAM file (absolute path)",type=str)
    parser.add_argument("-o","--outfile",help="Output file name (outfile.sam)",type=str)
    parser.add_argument("-u","--umi",help="Input file with known UMIs (absolute path)",type=str)
    parser.add_argument("-h","--help",help="HELP",type=str)
    return parser.parse_args()

#Get current record and return necessary variables using file handle input (input_fh) - see PS8/Demux.py for code
def getRecord(input_fh):
#     * Save current record to variable
#     * Split record for each required variable. Index positions per variable after record split: QNAME = 0, POS = 3, RNAME = 2, FLAG = 1, CIGAR = 5)
#     * Split QNAME again by ":" and get UMI (index 7)
#     * Return list of variables (current record, UMI, POS, RNAME, FLAG, CIGAR)
#     * i.e. Input = fh for test.sam, Output = ["NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	71M...", CTGTTCAC, 76814284, 2, 0, 71M]
    pass

#Open and write current record to output SAM file specified by argparse -o (output_fh) - see Demux.py for code
def writeRecord(output_fh,current_record):
#     * Prints current (complete) record to out.sam file
#     * Returns NULL
#     * i.e. Input = fh for out.sam, "NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	71M...", Output = NULL (but creates/prints to output file out.sam)
    pass

#Takes input bitwise flag of current record and returns strandedness.
def getStrand(flag):
#     * Evaluates bit 16 -> (flag & 16) == 16. TRUE = "reverse", FALSE = "forward"
#     * Returns "forward" or "reverse"
#     * i.e. Input = 16, Output = "reverse"
    pass

#Takes input CIGAR string and original pos and returns corrected position (adjusted for soft-clipping).
def processCIGAR(pos,cigar):
#     * Parses cigar string to identify location (beginning or end of read) and amount (#bases) of soft-clipping
#     * Adjusts pos based on cigar info
#     * Returns adj_pos
# need to add deletions and ns for rev (-) read strands - split into two functions
# Rev strand reads: (stpos + ds + ns + ending ss)
# For strand reads: (beg ss + stpos)
# Ms, Is, Ss, give length of read
#     * i.e. Input = 100,2S50M, Output = 102
    pass

#Main Pseudocode:

#Retrieve argparse inputs (see Functions)
args = get_args()
in_file = args.file #i.e. 
out_file = args.outfile #i.e. "outfile.sam"
umi_file = args.umi #i.e. "STL96.txt"


# * Initialize variables for UMI/QNAME, position/POS, chromosome/RNAME, bitwise flag, and CIGAR string (will be assigned using getRecord())
# * Initialize set to hold known UMIs and populate set with each line in STL96.txt (-u argparse)
#     * Use while True loop to go through each line in file, remove \n, save UMI string to set
# * Initialize dictionary to track written records. Key = (UMI,RNAME,strand), Value = (POS,CIGAR).
# * Use "while True" loop to go through each line of input sorted (by QNAME) SAM file (refer to PS8/Demux.py for code). For each line in SAM file (until end of file):
#     * Check if line begins w/ "@"
#         * Open output file (-o argparse)
#         * Write current line
#     * Else: getRecord(input_fh) - provides list of [current record, UMI, POS, RNAME, FLAG, CIGAR] (replaced w/ each loop iteration). Referred to as recordList below.
#     * Check if current UMI (recordList[1]) is in set of known UMIs. If yes:
#         * Check if UMI (recordList[1]) is in "written" dict as key[0]. If yes:
#             * Check if RNAME (recordList[3]) is in "written" dict as key[1]. If yes:
#                 * Check if getStrand(recordList[4]) equals strandeness value in "written" dict as key[2]. If yes:
#                     * Check if CIGAR (recordList[5]) contains "S" (for soft-clipping). If yes:
#                         * Get adjusted position from processCIGAR(recordList[2],recordList[5])
#                         * Check if adj_POS equals "written" dict value[0] for current (UMI,RNAME,strand) key. If yes:
#                             * THIS IS A DUPLICATE - DO NOT WRITE TO FILE
#                         * If no (POS is different):
#                             * Open output file (-o argparse)
#                             * writeRecord(output_fh,recordList[0])
#                             * Add necessary variables of current record to "written" dict (key=(UMI,RNAME,getStrand(recordList[4])), value=(adj_POS,recordList[5]))
#                     * If no (no soft-clipping was applied for this record):
#                         * Check if POS (recordList[2]) equals "written" dict value[0] for current (UMI,RNAME,strand) key. If yes:
#                             * THIS IS A DUPLICATE - DO NOT WRITE TO FILE
#                         * If no (POS is different):
#                             * Open output file (-o argparse)
#                             * writeRecord(output_fh,recordList[0])
#                             * Add necessary variables of current record to "written" dict (key=(UMI,RNAME,getStrand(recordList[4])), value=(recordList[2],recordList[5])) 
#                 * If no (strandedness is different):
#                     * Open output file (-o argparse)
#                     * writeRecord(output_fh,recordList[0])
#                     * Add necessary variables of current record to "written" dict (key=(UMI,RNAME,getStrand(recordList[4])), value=(recordList[2],recordList[5])) 
#             * If no (current record has new RNAME):
#                 * Open output file (-o argparse)
#                 * writeRecord(output_fh,recordList[0])
#                 * Add necessary variables of current record to "written" dict (key=(UMI,RNAME,getStrand(recordList[4])), value=(recordList[2],recordList[5])) 
#         * If no (current record has new, known UMI):
#             * Open output file (-o argparse)
#             * writeRecord(output_fh,recordList[0])
#             * Add necessary variables of current record to "written" dict (key=(UMI,RNAME,getStrand(recordList[4])), value=(recordList[2],recordList[5]))