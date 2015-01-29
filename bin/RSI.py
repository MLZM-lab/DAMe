#!/usr/bin/python

import sys
import numpy
from optparse import OptionParser

## Get options and display help
usage = "usage: %prog [-h help] [-e explicit] [-o output FILE] <in.txt>"
parser = OptionParser(usage)
        
parser.add_option("-e", "--explicit", action="store_true", dest="explicit",
                  help="The output file contains the explicit RSI value for every pairwise comparison.")
parser.add_option("-o", "--output", dest="outfile", metavar="FILE", help="Write output to FILE.")
                  
(options, args) = parser.parse_args()
if not len(args):
	parser.error('You have to provide an input file.')
	

## Defining 'compare' function
def compare(matrix,j,a,b):
	print 'Comparing replicates '+repr(a)+' and '+repr(b) + ' from sample '+j 
	CP = []
	sum = matrix.sum(axis=0)
	if sum[0] == 0:
		sum[0] = 1
		print 'This sample gave zero in replicate: '+repr(a)
	if sum[1] == 0:
		sum[1] = 1
		print 'This sample gave zero in replicate: '+repr(b)
	a = numpy.array([row[0]/float(sum[0]) for row in matrix])
	b = numpy.array([row[1]/float(sum[1]) for row in matrix])
	percent = numpy.column_stack((a,b))
	return(1-numpy.sum([row.min() for row in percent]))
	

## Main
data = []
output = []
rkn = []

## Reading data
with open(args[0]) as file:
	for line in file:
		parts = line.split()
		data.append([q for q in parts])
data = numpy.array(data)
names = set([row[0] for row in data])
no_rep = (len(data[0])-2)/2
if no_rep < 2:
	print('There are no replicates in the file.')
	sys.exit(0)

## Calculate pairwise RSI values
for i in names:
	rep = 0
	output = 0
	subset = numpy.array([row for row in data if row[0]==i])
	## one by one
	if options.explicit:
		for A in range(1,no_rep):
			for B in range(A+1,no_rep+1):
				sample = subset[:, [A*2,B*2]]
				sample = sample.astype(int)
				output = compare(sample,i,A,B)
				rkn.append([i,A,B,output])
	## by sample
	else:
		for A in range(1,no_rep):
			for B in range(A+1,no_rep+1):
				sample = subset[:, [A*2,B*2]]
				sample = sample.astype(int)
				output += compare(sample,i,A,B)
				rep=rep+1
		rkn.append([i,output/rep])

## Write output file
if options.outfile:
	outfile = options.outfile
else:
	outfile='RSI_output.txt'

with open(outfile,'w') as out_file:
	## one by one
	if options.explicit:
		out_file.write('Sample\tReplicateA\tReplicateB\tRSI\n')
		for row in rkn:
			out_file.write(str(row[0])+'\t'+str(row[1])+'\t'+str(row[2])+'\t'+str(row[3])+'\n')
	## by sample
	else:		
		out_file.write('Sample\tRSI\n')
		for row in rkn:
			out_file.write(str(row[0])+'\t'+str(row[1])+'\n')

