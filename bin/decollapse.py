#! /usr/bin/env python

## Marie Lisandra Zepeda Mendoza
## v 1
## 24.10.14

####################
###Get the options
import argparse
parser = argparse.ArgumentParser(description='This program puts X number of times a sequence reported like this: Primer	Tag1	Tag2	Freq	Seq	Qual	Header. The utput is a fasta withe >Tag1.Tag2.Freq_RN , where RN is the writen line number in the output.')
parser.add_argument("-input", required=True, help='Text file with the information on the tag combination and freq of each unique seq')
parser.add_argument("-outFas", default="Decollapsed.fasta", help='Output fasta file name with the unique sequences being repeated as many times as their reported freq [default "Decollapsed.fasta"]')

args = parser.parse_args()

####################
import sys
from modules_filter import *

######################
##Initialize variables
IN=args.input
OUT=args.outFas

##################### Main program

IN=open(IN)
OUT=open(OUT, "w")

seq_id=0
line=IN.readline()
while line:
	line=line.rstrip().split()
	count=1
	while count<=int(line[3]):
		seq_id+=1
		OUT.write( ">" + line[1] + "." + line[2] + "." + line[3] + "_" + str(seq_id) + "\n" + line[4] + "\n")
		count+=1
	line=IN.readline()

IN.close()
OUT.close()

