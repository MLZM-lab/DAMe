
def makePSnumFiles(PSinfo, X, P, chimeraChecked):
	#Open the PS#.txt out handlers
	PSouts=[]
	for i in range(X):
		PSouts.append(open("PS%s_files.txt"%(i+1), "w"))
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
			if not chimeraChecked:
				if P==1:
					PSouts[residue-1].write("%s_%s.txt\n"%(psinfo[1], psinfo[2]))
				else:  
					PSouts[residue-1].write("./pool%s/%s_%s.txt\n"%(psinfo[3], psinfo[1], psinfo[2]))
			else:
				if P==1:
					PSouts[residue-1].write("%s_%s.noChim.txt\n"%(psinfo[1], psinfo[2]))
				else:  
					PSouts[residue-1].write("./pool%s/%s_%s.noChim.txt\n"%(psinfo[3], psinfo[1], psinfo[2]))
		else:
			if not chimeraChecked:
				if P==1:
					PSouts[X-1].write("%s_%s.txt\n"%(psinfo[1], psinfo[2]))
				else:  
					PSouts[X-1].write("./pool%s/%s_%s.txt\n"%(psinfo[3], psinfo[1], psinfo[2]))
			else:
				if P==1:
					PSouts[X-1].write("%s_%s.noChim.txt\n"%(psinfo[1], psinfo[2]))
				else:  
					PSouts[X-1].write("./pool%s/%s_%s.noChim.txt\n"%(psinfo[3], psinfo[1], psinfo[2]))
	#Close handlers
	for i in range(len(PSouts)):
		PSouts[i].close()



def ReadPSnumFiles(X):
	PSins=[]
	for i in range(X):
		PSins.append(open("PS%s_files.txt"%(i+1)))
	PSinsLines={}
	for i in range(X): 
		PSinsLines[str(i)]=PSins[i].readlines()
		PSins[i].close()		
	return(PSinsLines)


def MakeSampleNameArray(PSinfo):
	sampleName=[]
	file=open(PSinfo)
	line=file.readline()
	while line:
		if not line.split()[0] in sampleName:
			sampleName.append(line.split()[0])
		line=file.readline()
	file.close()
	return(sampleName)



def ReadHapsForASample(X, PSinsLines, i ):
	haps={}
	for j in range(X):
		haps[str(j)] = list()
		if PSinsLines[str(j)][i].rstrip() != "empty" :
			f = open(PSinsLines[str(j)][i].rstrip())
			line= f.readline()
			while line:
				haps[str(j)].append(line.split())
				line= f.readline()
	return(haps)


def getSeqsSetsAndFRcounts(X, haps):
	F={}
	R={}
	counts={}
	seqs={}
	seqsALL=[]
	for j in range(X):
		if len(haps[str(j)])!=0: #if it's not empty
			seqs[str(j)]=[]
			F[str(j)]=haps[str(j)][0][1]
			R[str(j)]=haps[str(j)][0][2]
			counts[str(j)]=[]
			for k in range(len(haps[str(j)])):
				counts[str(j)].append(haps[str(j)][k][3])
				seqs[str(j)].append(haps[str(j)][k][4])
				seqsALL.append(haps[str(j)][k][4])
	seqsALL=set(seqsALL) 
	return (seqsALL, F, R, counts, seqs)


def MakeComparisonFile(X, seqsALL, haps, F, R, counts, seqs, OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas, OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i):
	idnum=1
	for seq in seqsALL:
		line=sampleName[i] + "\t"
		lineFasIDs=">" + sampleName[i] + "\t"
		lineFasCounts= "\t"
		y=0
		t=0
		for j in range(X):
			if len(haps[str(j)])!=0: #if it's not empty
				pos = [pos for pos,s in enumerate(seqs[str(j)]) if seq == s] ##Get the pos where it overlaps
				if len(pos)==0:
					count=0
				else:
					y=y+1
					count=counts[str(j)][pos[0]]
					if int(count) < T :
						t=t+1
				line = line + F[str(j)] + "-" + R[str(j)] + "\t" + str(count) + "\t"
				if j < (X-1):
					lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "." 
					lineFasCounts= lineFasCounts + str(count) + "_"
				else:
					lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "_" + str(idnum) + "\t"
					lineFasCounts= lineFasCounts + str(count) + "\n" + seq
			if len(haps[str(j)])==0: #if it is empty
				line = line + "empty\t0\t"
				if j < (X-1):
					lineFasIDs = lineFasIDs + "empty-empty." 
					lineFasCounts= lineFasCounts + "0_"
				else:
					lineFasIDs = lineFasIDs + "empty-empty_" + str(idnum) + "\t"
					lineFasCounts= lineFasCounts + "0\n" + seq
		line= line + seq + "\n"
		lineFas= lineFasIDs + lineFasCounts + "\n"
		OUT.write(line)
		OUT_fas.write(lineFas)
		if y >= Y :
			OUTYX.write(line)
			OUTYX_fas.write(lineFas)
		if (y-t) >= Y :
			OUTthresh.write(line)
			OUTthresh_fas.write(lineFas)
			if len(seq) >= L:
				 OUTthreshLen_fas.write(lineFas) 		
		idnum=idnum+1


