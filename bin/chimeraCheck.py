#! /usr/bin/env python

## Marie Lisandra Zepeda Mendoza
## v 1
## 24.10.14

####################
###Get the options
import argparse

parser = argparse.ArgumentParser(description='Create necessary files to operate on sequences per PCR reaction')
parser.add_argument("-psInfo", required=True, help='Text file with the information on the tag combination in each PCR reaction for every sample [Format: sampleName\tTagNameForwardFromPCR1\tTagNameReverseFromPCR1\tPool#\nsampleName\tTagNameForwardFromPCR2\tTagNameReverseFromPCR2\tPool#\n...]')
parser.add_argument("-x", required=True, type=int, help='Number of PCR rxns performed per sample')
parser.add_argument("-p", type=int, default=1, help='The number of pools in which the samples were divided for sequencing (in case of tag combinations repetition due to processing so many samples) [default 1]\nNOTE: If using pools, each fastq must be in a folder called pool#, in which the sort.py was run for each pool inside the corresponding folder, and this program chimeraCheck.py is run in the parent directory of the pools directories')

args = parser.parse_args()

####################
import sys
from modules_chimeraCheck import *

######################
##Initialize variables
PSinfo= args.psInfo 
X= args.x
P= args.p


######################
## Main program
######################

#Create tag files

if P==1:
	makeTagFiles(PSinfo, X)
else:
	makeTagFilesWithPools(PSinfo, X)

#Make the sizeout fasta for each PCR rxn (which is each tag combination from each sample)
MakeSizeOutFastas(P, X)

#Sort the sizeout fastas
SortFasta(P)

## Bring the noChim fasta back to the hap format fileS... sort it into the tagCombs
#Make the fastas one liners
MakeFasSeqOneLine(P)
#sort it into the tagCombs  
MakeNoChimHaps(P)


