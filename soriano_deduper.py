#!/usr/bin/env python

#Imports:
import argparse
import re

# Functions:

#Accept command line input for fields mentioned in assignment README.md - *ADD HELP LINE*
def get_args():
    '''Argparse Function: Retrieves user-input arguments from command line.'''
    parser = argparse.ArgumentParser(description="A program to remove duplicate reads from a given SAM file. Please use the following format to use this program to deduplicate SAM files: ./soriano_deduper.py -f INPUT_FILENAME.sam -o OUTPUT_FILENAME.sam -u KNOWN_UMI_FILE.txt")

    parser.add_argument("-f","--file",help="Input sorted SAM file (absolute path)",type=str)
    parser.add_argument("-o","--outfile",help="Output file name (outfile.sam)",type=str)
    parser.add_argument("-u","--umi",help="Input file with known UMIs (absolute path)",type=str)
    return parser.parse_args()

#Get current record and return necessary variables using file handle input (input_fh)
def getRecord(input_fh):
    '''Function that takes the file handle of the SAM input file and returns a list of the relevant sections of the current record/line in SAM file. Line sections returned as record_list = [UMI,BITFLAG,RNAME,POS,CIGAR,complete_record(w/o newline)].'''
    record_stripped: str = input_fh.readline().strip()
    # Split record for each required variable. Index positions per variable after split: QNAME = 0, POS = 3, RNAME = 2, FLAG = 1, CIGAR = 5.
    record_split: list = record_stripped.split("\t")
    if record_split[0] == "":
        record_list: list = ["EOF"]
    elif record_split[0].startswith("@") == True:
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
    '''Takes output file handle and current, complete record line and writes current record to output file.'''
    output_fh.write(current_line +"\n")
    return None

#Takes input bitwise flag of current record and returns strandedness (forward or reverse).
def getStrand(flag):
    '''Takes bitwise flag (from current record) and returns read strandedness ("reverse" or "forward") by evaluating bit 16.'''
    int_flag: int = int(flag)
    strandedness: str = ""
    if((int_flag & 16) == 16):
        strandedness = "reverse"
    else:
        strandedness = "forward"
    return strandedness

#Takes input CIGAR string and original pos and returns corrected position (adjusted for soft-clipping).
def processCIGAR(pos,cigar,strand):
    '''Takes read start position (POS), cigar string (CIGAR), and strandedness (from getStrand(bitflag)) of current record and returns adjusted read start position. Accounts for how cigar and strandeness affect actual read start position when checking for duplicate reads.'''
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
        adj_pos = pos - beg_s
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
written_reads = set()
# writ_umis = set()
# writ_rns = set()

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
        elif sam_rec_list[0] == "EOF":
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
                rec_list = (sam_rec_list[0],sam_rec_list[2],current_strand,pos_adj)
                #print(rec_list)
                if rec_list not in written_reads:
                    #Yes - current record is NOT in written_reads set
                    writeRecord(fh_out,sam_rec_list[5])
                    #Add current read (rec_list) to written_reads set
                    written_reads.add(rec_list)

        # if i>99: #FOR TESTING
        #     break #FOR TESTING
#print(written_reads) #FOR TESTING

fh_out.close()