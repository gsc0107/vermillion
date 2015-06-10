#! /usr/bin/env python
import sys
from Bio import SeqIO


#script takes blast output and fasta files for reads hitting to a particular genome (but not hitting only to that, i.e. likely a viral element in a larger genome)
#script collates pairs of reads that have hits to the genome of interest and then trims off those portions, leaving only read sections that hit to the external genome
#usage
#python trimInternalGenome.py read1BlastFile read2BlastFile read1FastaFile read2FastaFile outputName

#get the list of files from the command line and the name for the output file
try:
	read1=sys.argv[1]
except IndexError:
	print "\n read 1 blast output file not supplied."
	sys.exit()
except IOError:
	print "\n read 1 blast output file not found in directory."
	sys.exit()
try:
	read2=sys.argv[2]
except IndexError:
	print "\n read 2 blast output file not supplied."
	sys.exit()
except IOError:
	print "\n read 2 blast output file not found in directory."
	sys.exit()
try:
	fasta1=sys.argv[3]
except IndexError:
	print "\n read 1 fasta file not supplied."
	sys.exit()
except IOError:
	print "\n read 1 fasta file not found in directory."
	sys.exit()
try:
	fasta2=sys.argv[4]
except IndexError:
	print "\n read 2 fasta  file not supplied."
	sys.exit()
except IOError:
	print "\n read 2 fasta  file not found in directory."
	sys.exit()
try:
	out=sys.argv[5]
except IndexError:
	print "\n output filename file not supplied."
	sys.exit()

#populate a dictionary with the linked blast results between read 1 and read 2
observed = {} # make an empty dictionary observed
file = open(read1, "r")# open a blast file for read one - this has been formatted to contain only one hit per query and show only hits
for line in file:
	line=line.split("\t")
	observed[line[0]]=[int(line[6]),int(line[7])],[0,0] #make the header the key in the dictionary and set the query start and stop nucleotides, set a place holder for the second reads start and stop coordinates
file.close()

file = open(read2, "r") # open the blast file for pair 2
for line in file:
	line=line.split("\t")
	if observed.has_key(line [0]): # if the header is already in the dictionary add start and stop values
		observed[line[0]][1][0] = int(line[6])
		observed[line[0]][1][1] = int(line[7])
	else:
		observed[line[0]]=[0,0],[ int(line[6]), int(line[7])]# now make the dictionary of the fasta reads # if the header is not in the dictionary add 00 for read one and start stop for read 2
file.close()



#parse the results to find only read pairs that are related to the genome of interest (ALVE)
reads = {} # create an empty dictionary reads
for record in SeqIO.parse(fasta1, 'fasta'): # open pair one fasta file
	if observed.has_key(record.id): # if the header from the fasta file is in the dictionary observed do the following:
		header = record.id
		sequence = str(record.seq)
		if observed[header][0][0]==0 and observed[header][0][1]==0: # if the start and stop coordinates are zero ie. no hit to ALVE therefore all chicken 
			reads[record.id]=[sequence,0] # take the entire sequence and put in reads dictionary
		elif observed[header][0][0] > 0 and observed[header][0][1] > 0: # if start and stop were not zero - ie. there was a hit to ALVE
			substring=sequence[:observed[header][0][0]-1]+sequence[observed[header][0][1]:] # trim the ALVE portions and set it as substring
			if len(substring)>=20: # if the substring is > or = 20 put the substring in the reads dictionary
				reads[record.id]=[substring,0] # hold a place for the second read
			else:
				reads[record.id]=[0,0] 
for record in SeqIO.parse(fasta2, 'fasta'): #open pair two fasta file
	if observed.has_key(record.id): #if the header is in the observed dictionary do the following
		header = record.id
		sequence = str(record.seq)
		if observed[header][1][0]==0 and observed[header][1][1]==0: # no hit ie. all  chicken
			reads[record.id][1]=sequence 
		elif observed[header][1][0] > 0 and observed[header][1][1] > 0: # partial or complete ALVE hit
			substring=sequence[:observed[header][1][0]-1]+sequence[observed[header][1][1]:] #trim away ALVE
			if len(substring)>=20:
				reads[record.id][1]=substring

#output the reads of interest to file				
outfile= open(out, 'w')#this code writes the output fasta file of dictionary reads
for header in reads:
	if reads[header][0]!=0 and reads[header][1]!=0:
		outfile.write('>'+header+'_1\n') # print > to recognize as fasta plus the header with a "_1" to show it is read 1 plus end of line
		outfile.write(reads[header][0]+'\n') #print the first sequence and end of line
		outfile.write('>'+header+'_2\n') # print header in fasta format with "_2" to show read 2 and end of line
		outfile.write(reads[header][1]+'\n')
	elif reads[header][0]==0 and reads[header][1]!=0:
		outfile.write('>'+header+'_2\n') # print header in fasta format with "_2" to show read 2 and end of line
		outfile.write(reads[header][1]+'\n')
	elif reads[header][0]!=0 and reads[header][1]==0:
		outfile.write('>'+header+'_1\n') # print > to recognize as fasta plus the header with a "_1" to show it is read 1 plus end of line
		outfile.write(reads[header][0]+'\n') #print the first sequence and end of line
	elif reads[header][0]==0 and reads[header][1]==0:
		continue
outfile.close()
sys.exit()
