#! /usr/bin/env python

## Marie Lisandra Zepeda Mendoza
## v 1
## 24.10.14

####################
###Get the options
import argparse
parser = argparse.ArgumentParser(description='Filter the multiplexed sequences by the presence of the sequence in the different PCR reactions at a given minimum abundance, and its length')
parser.add_argument("-psInfo", required=True, help='Text file with the information on the tag combination in each PCR reaction for every sample [Format: sampleName\tTagNameForwardFromPCR1\tTagNameReverseFromPCR1\tPool#\nsampleName\tTagNameForwardFromPCR2\tTagNameReverseFromPCR2\tPool#\n...]')
parser.add_argument("-x", type=int,default=2,  help='Number of PCR rxns performed per sample')
parser.add_argument("-y",type=int, default=1, help='Number of PCR rxns in which the sequence has to be present')
parser.add_argument("-p", type=int, default=1, help='The number of pools in which the samples were divided for sequencing (in case of tag combinations repetition due to processing so many samples) [default 1]\nNOTE: If using pools, each fastq must be in a folder called pool#, in which the sort.py was run for each pool inside the corresponding folder, and this program chimeraCheck.py is run in the parent directory of the pools directories')
parser.add_argument("-t",type=int, default=1, help='Number of times a unique sequence has to be present')
parser.add_argument("-l",type=int, default=100, help='Minimum sequence length')
parser.add_argument("--chimeraChecked",  help='Use this parameter if you have performed a chimera check on the sorted collapsed sequence files [default not set]', action="store_true")

args = parser.parse_args()

####################
import sys
from modules_filter import *

######################
##Initialize variables
PSinfo= args.psInfo 
X= args.x
Y= args.y
P= args.p
T= args.t
L= args.l

chimeraChecked=args.chimeraChecked

OUT=open("Comparisons_%sPCRs.txt"%(X), "w")  
OUTYX=open("Comparisons_%soutOf%sPCRs.txt"%( Y, X), "w")  
OUTthresh=open("Comparisons_%soutOf%sPCRs.countsThreshold%s.txt"%( Y, X, T), "w")  
OUT_fas=open("Comparisons_%sPCRs.fasta"%( X), "w")  
OUTYX_fas=open("FilteredReads_atLeast%s.fasta"%(Y), "w")  
OUTthresh_fas=open("FilteredReads_atLeast%s.threshold.fasta"%( Y), "w")  
OUTthreshLen_fas=open("FilteredReads.fna", "w") # The one where the threshold reads were further filtered by length. The final one to keep using 


#################################################################### Main program

#### Make the PS#.txt files OUTs of the PSinfo.txt file
makePSnumFiles(PSinfo, X, P, chimeraChecked)

##### Make the comparison of how many counts are per PSset per seq 

#Open the PS#_files.txt IN handlers and read them
PSinsLines=ReadPSnumFiles(X)

#Get the names of the samples
sampleName=MakeSampleNameArray(PSinfo)

# Compare all the samples
for i in range(len(PSinsLines["1"])):
	#Read the haps files that are on PSinsLines paths for sample A
	haps=ReadHapsForASample(X, PSinsLines, i )
	#Get all the seqs (as well as its F, R, and counts info) for each PSnum files and get a nr set of them (seqsALL)
	seqsALL, F, R, counts, seqs = getSeqsSetsAndFRcounts(X, haps)
	#Compare how many times each seq from seqsALL are in each PSnumber and get the PSs seqs passing the thresholds of Y out of X 
	MakeComparisonFile(X, seqsALL, haps, F, R, counts, seqs, OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas, OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i)

# Close the output files
OUTthresh.close()
OUT.close()
OUTYX.close()
OUTthresh_fas.close()
OUT_fas.close()
OUTYX_fas.close()
OUTthreshLen_fas.close()

