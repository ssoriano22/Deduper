#!/usr/bin/env python

#Imports:
import argparse
import re

# Functions:

#Accept command line input for fields mentioned in assignment README.md - *ADD HELP LINE*
def get_args():
    parser = argparse.ArgumentParser(description="A program to remove duplicate reads from a given SAM file.")
    parser.add_argument("-f","--file",help="Input sorted SAM file (absolute path)",type=str)
    parser.add_argument("-o","--outfile",help="Output file name (outfile.sam)",type=str)
    parser.add_argument("-u","--umi",help="Input file with known UMIs (absolute path)",type=str)
    #parser.add_argument("-h","--help",help="HELP",type=str)
    return parser.parse_args()

#Get current record and return necessary variables using file handle input (input_fh)
def getRecord(input_fh):
    record_stripped: str = input_fh.readline().strip()
    # Split record for each required variable. Index positions per variable after split: QNAME = 0, POS = 3, RNAME = 2, FLAG = 1, CIGAR = 5.
    record_split: list = record_stripped.split("\t")
    if record_split[0].startswith("@") == True:
        #Header case
        record_list: list = ["HEADER",record_stripped]
    else:
        #All other lines
        #Split QNAME again by ":" and get UMI (index 7)
        umi: str = record_split[0].split(":")[7]
        bitflag: int = record_split[1]
        rname: str = record_split[2]
        pos: int = record_split[3]
        cigar: str = record_split[5]
        record_list: list = [umi,bitflag,rname,pos,cigar,record_stripped]
    return record_list

#Open and write current record to output SAM file specified by argparse -o (output_fh)
def writeRecord(output_fh,current_line):
    output_fh.write(current_line +"\n")
    return None

#Takes input bitwise flag of current record and returns strandedness (forward or reverse).
def getStrand(flag):
    int_flag: int = int(flag)
    strandedness: str = ""
    if((int_flag & 16) == 16):
        strandedness = "reverse"
    else:
        strandedness = "forward"
    return strandedness

#Takes input CIGAR string and original pos and returns corrected position (adjusted for soft-clipping).
def processCIGAR(pos,cigar,strand):
    adj_pos: int = pos #Initialize adj_pos to initial position from current record
    #Initialize relevant cigar options to 0
    dels: int = 0
    nns: int = 0
    end_s: int = 0
    beg_s: int = 0
    #Split cigar w/ regex into list of letters and list of numbers (indexes should be compatible btwn lists)
    cigar_letters: list = re.split("\d+",cigar)[1:]
    #print(cigar_letters)
    cigar_numbers: list = re.split("[A-Z]+",cigar)[:-1]
    #print(cigar_numbers)
    if strand == "reverse":
        #Rev strand reads: (stpos + ds + ns + ending ss)
        if "D" in cigar_letters:
            #Yes - cigar contains Ds
            D_index: int = cigar_letters.index("D")
            dels = int(cigar_numbers[D_index])
        if "N" in cigar_letters:
            #Yes - cigar contains Ns
            N_index: int = cigar_letters.index("N")
            nns = int(cigar_numbers[N_index])
        if cigar_letters[-1] == "S":
            #Yes - cigar ends with S
            end_s = int(cigar_numbers[-1])
        adj_pos = pos + dels + nns + end_s
    else:
        #For strand reads: (beg ss + stpos)
        if cigar_letters[0] == "S":
            beg_s = int(cigar_numbers[0])
        adj_pos = beg_s + pos
    return adj_pos
#print(processCIGAR(5,"3S2I3D71M3S","reverse")) #FOR TESTING CIGAR FUNCTION

#MAIN CODE:

#Retrieve argparse inputs (see Functions)
args = get_args()
in_file = args.file #i.e. test.sam
out_file = args.outfile #i.e. "outfile.sam"
umi_file = args.umi #i.e. "STL96.txt"

#Initialize set to hold known UMIs
known_umis = set()
with open(umi_file,"r") as fh_umi:
    for line in fh_umi: #Populate set with each line in umi file
        k_umi: str = line.strip()
        known_umis.add(k_umi)

#Initialize list to track written records. Format = (UMI,RNAME,strand,adj POS). *WHAT STRUCTURE SHOULD THIS BE*
written_reads = []
writ_umis = set()
writ_rns = set()

#Open output file for writing
fh_out = open(out_file,"w")

i = 0 #Sam line/record counter
with open(in_file,"r") as fh_input:
    while True:
        i+=1 #Increment record counter (start at i=1)
        sam_rec_list: list = getRecord(fh_input)
        if sam_rec_list[0] == "HEADER":
            #Header line cases - print all headers to output file
            writeRecord(fh_out,sam_rec_list[1])
        elif sam_rec_list[0] == "":
            #EOF case - break while True loop for reading input file
            break 
        else:
            #Read line - evaluate for duplicate characteristics
            if sam_rec_list[0] in known_umis:
                #Yes - current UMI is in list of known UMIs
                #Determine strandedness from bitwise flag for written set - needed to add written_reads to list going forward
                current_strand: str = getStrand(sam_rec_list[1])
                #Determine adjusted position based on cigar - needed to add written_reads to list going forward
                pos_adj: int = processCIGAR(int(sam_rec_list[3]),sam_rec_list[4],current_strand)
                #Format current read to a list w/ same parameters as the lists in written_reads list
                rec_list: list = [sam_rec_list[0],sam_rec_list[2],current_strand,pos_adj] 
                #print(rec_list)
                if rec_list[0] not in writ_umis:
                    #Yes - current UMI is NOT in written_reads list
                    writeRecord(fh_out,sam_rec_list[5])
                    #Add current read (rec_list) to written_reads list
                    written_reads.append(rec_list)
                    #Add current umi/rname to the written sets for each
                    writ_umis.add(rec_list[0])
                    writ_rns.add(rec_list[1])
                else:
                    #No - current UMI has been written before
                    writ_rns = set() #Create set containing currently written rnames
                    for written_rec in written_reads:
                        writ_rns.add(written_rec[1])
                    if rec_list[1] not in writ_rns:
                        #Yes - current UMI is NOT in written_reads list
                        writeRecord(fh_out,sam_rec_list[5])
                        #Add current read (rec_list) to written_reads list
                        written_reads.append(rec_list)
                    else:
                        pass
            else:
                #No - current UMI is unknown/error
                pass #DO NOT PRINT RECORD
        if i>27: #FOR TESTING
            break #FOR TESTING
#print(written_reads)
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

fh_out.close()