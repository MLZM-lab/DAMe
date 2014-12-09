## Marie Lisandra Zepeda Mendoza
## v 1
## 24.10.14

################Import libraries

import re
import os
import sys
import subprocess

################Define the functions


def makeTagFiles(PSinfo, X):
	#Open the PS#.txt out handlers
	PSouts=[]
	for i in range(X):
		PSouts.append(open("PS%s.tags.txt"%(i+1), "w"))
	#Fill the opened handlers. Parse the PSinfo.txt
	PSinfo=open(PSinfo)
	PS=PSinfo.readlines()
	PSinfo.close()
	for NR, psinfo in enumerate(PS):
		NR=NR+1
		psinfo=psinfo.rstrip() #The chomp equivalent
		psinfo= psinfo.split()
		residue=NR%X
		if residue != 0:
			PSouts[residue-1].write("%s\t%s\n"%(psinfo[1], psinfo[2]))
		else:
			PSouts[X-1].write("%s\t%s\n"%(psinfo[1], psinfo[2]))
	#Close handlers
	for i in range(len(PSouts)):
		PSouts[i].close()


def makeTagFilesWithPools(PSinfo, X):
	#Open the PS#.txt out handlers
	PSouts=[]
	for i in range(X):
		PSouts.append(open("PS%s.tags.txt"%(i+1), "w"))
	#Fill the opened handlers. Parse the PSinfo.txt
	PSinfo=open(PSinfo)
	PS=PSinfo.readlines()
	PSinfo.close()
	for NR, psinfo in enumerate(PS):
		NR=NR+1
		psinfo=psinfo.rstrip() #The chomp equivalent
		psinfo= psinfo.split()
		residue=NR%X
		if residue != 0:
			PSouts[residue-1].write("%s\t%s\t%s\n"%(psinfo[1], psinfo[2], psinfo[3]))
		else:
			PSouts[X-1].write("%s\t%s\t%s\n"%(psinfo[1], psinfo[2], psinfo[3]))
	#Close handlers
	for i in range(len(PSouts)):
		PSouts[i].close()




def MakeSizeOutFastas(P, X):    
	#Open a fasta for each pool (the samples pooled are all amplified togeter in a PCR)
	OUTS=[]
	for pool in range(P):
		OUTS.append(open("Pool%s.fasta"%(pool+1), "w"))
	for num in range(X):
		#open PS${num}.tags.txt
		file=open("PS%s.tags.txt"%(num+1)) 
		line= file.readline()	
		while line:
			line=line.rstrip().split()
			hap="_".join([line[0], line[1]])
			hap=".".join([hap,"txt"])
			if P>1:
				pool= "./pool" + str(line[2]) + "/"
				hap= pool + hap
			#Open the sorted collapsed countes seqs files and fill the pool fasta
			if os.path.exists(hap):
				hap=open(hap)
			else:	
				line=file.readline()
				continue
			seq= hap.readline()
			idNum=1 
			while seq:
				seq=seq.rstrip().split()
				a= ">"+ "_".join([ seq[0], line[0], line[1], str(idNum)]) + ";size=" + str(seq[3]) + "\n" +  seq[4] ## >Gene_Tag1_Tag2_idNum;size=Freq\nSeq
				if P>1:
					OUTS[int(line[2])-1].write("%s\n"%(str(a)))
				else:
					OUTS[0].write("%s\n"%(a))
				idNum+=1
				seq= hap.readline()
			hap.close()
			line=file.readline()  
		file.close()
	#Now close the pools fastas
	for pool in range(P):
		OUTS[pool].close()



def SortFasta(P):
	for pool in range(P):
		input="Pool%s.fasta"%(pool+1)
		output="Pool%s.sort.fasta"%(pool+1)
		cmd= "usearch --sortsize " + input + " --output " + output  #First I need to sort the file
		p_core= subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p_core.communicate()
		fh= open("sort%s.out"%(pool+1), 'w')
		fh.write(stdout)
		fh.close()
		fh= open("sort%s.err"%(pool+1), 'w')
		fh.write(stderr)
		fh.close()
		# Now check for chimeras
		input=output
		output1="Pool%s.Chim.fasta"%(pool+1)
		output2="Pool%s.noChim.fasta"%(pool+1)
		cmd=  "usearch -uchime " + input + " -chimeras " + output1 + " -nonchimeras " + output2   #Identify and remove chimeras de novo
		p_core= subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p_core.communicate()
		fh= open("chimeraCheck%s.out"%(pool+1), 'w')
		fh.write(stdout)
		fh.close()
		fh= open("chimeraCheck%s.err"%(pool+1), 'w')
		fh.write(stderr)
		fh.close()



def MakeFasSeqOneLine(P):
	for pool in range(P):
		fasta=open("Pool%s.noChim.fasta"%(pool+1))
		fastaOne=open("Pool%s.noChim.oneLiner.fasta"%(pool+1), "w")
		line= fasta.readline()	
		seq=""
		while line:
			line=line.rstrip()
			if line.find(">")==0:  #If it's a header
				if len(seq)>0: ## Print the seq if this is not the first entry
					fastaOne.write("%s\n"%(seq))
				fastaOne.write("%s\n"%(line))
				seq=""
			else:
				seq+=line
			line= fasta.readline()
		fastaOne.write("%s\n"%(seq)) #Print the last seq
		fastaOne.close()



def MakeNoChimHaps(P):
	HAP={}
	for pool in range(P):
		fasta=open("Pool%s.noChim.oneLiner.fasta"%(pool+1))
		line= fasta.readline()	
		while line:
			line=line.rstrip()
			if line.find(">")==0:  #If it's a header
				line=re.sub(">", "", line)
				primerName=line.split("_")[0]
				tagName1=line.split("_")[1]
				tagName2=line.split("_")[2]
				freq=line.split("=")[1]
				#Add tags to HAP if it doesn't have it already.  
				tagHapKey="_".join([tagName1, tagName2])
				if not HAP.has_key(tagHapKey):   ##key is the Tag1_Tag2  # When they are the first one
					HAP[tagHapKey]=[[],[],[],[],[]] # primerName, tag1, tag2, Freq, Seq  
					HAP[tagHapKey][0]=[primerName]
					HAP[tagHapKey][1]=[tagName1]
					HAP[tagHapKey][2]=[tagName2]
					HAP[tagHapKey][3]=[freq]
				else: # If they are not the first one
					HAP[tagHapKey][0].append(primerName)
					HAP[tagHapKey][1].append(tagName1)
					HAP[tagHapKey][2].append(tagName2)
					HAP[tagHapKey][3].append(freq)
			else: ##If its the seq
				if len(HAP[tagHapKey][4])==0: ## If it's the first one
					HAP[tagHapKey][4]=[line+"\n"]
				else:
					HAP[tagHapKey][4].append(line+"\n")
			line= fasta.readline()
		fasta.close()
	#Now print the HAP		
	for TagComb in HAP:
		out=open("%s.noChim.txt"%(TagComb), "w")
		for i in range(len(HAP[TagComb][0])):
			a="\t".join([ HAP[TagComb][0][i], HAP[TagComb][1][i], HAP[TagComb][2][i], HAP[TagComb][3][i], HAP[TagComb][4][i] ])
			out.write(a)
		out.close()


