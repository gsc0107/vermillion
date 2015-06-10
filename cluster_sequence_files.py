#! /usr/bin/env python
from Bio import SeqIO
from collections import defaultdict 
import sys

#script takes a file of clustered read headers (output from cluster_formation.py or similar) and creates a fasta file per cluster
#usage
#python cluster_files.py clusterFormationOutputFile filePrefix read1FastaFile read2FastaFile

#get the list of files from the command line and the name for the output file
try:
	cluster_file=sys.argv[1]
except IndexError:
	print "\n cluster file not supplied."
	sys.exit()
except IOError:
	print "\n cluster file not found in directory."
	sys.exit()
try:
	string=sys.argv[2]
except IndexError:
	print "\n output prefix not supplied."
	sys.exit()
try:
	fasta_p1=sys.argv[3]
except IndexError:
	print "\n read 1 fasta file not supplied."
	sys.exit()
except IOError:
	print "\n read 1 fasta file not found in directory."
	sys.exit()
try:
	fasta_p2=sys.argv[4]
except IndexError:
	print "\n read 2 fasta file not supplied."
	sys.exit()
except IOError:
	print "\n read 2 fasta file not found in directory."
	sys.exit()

#read in the cluster information
cluster=defaultdict(dict)
file = open(cluster_file, "r")
for line in file:
	line=line.rstrip()
	#start and chromosome are split by tabs so extract them
	#the headerStr is the headers separated by a space each
	start, chr, headerStr=line.split("\t")
	#take the header string and place into a list by separating by spaces
	headerList=headerStr.split(" ")
	cluster[str(chr)][int(start)] = headerList
file.close()

#for each cluster, extract the relevant fasta reads and output to file
for key in cluster:
	subcluster=cluster[key]
	for start in subcluster:
	#need to write my output file in here
		outfile=open(string + "_" + str(key) + "_" + str(start) +".fasta", 'w')
		headers=subcluster[start]
		# for h in headers:
		for record in SeqIO.parse(fasta_p1, 'fasta'): # open pair one fasta file
			if record.id in headers: # if the header from the fasta file is in the dictionary observed do the following:
				#print record.id + "!\n"
				header = record.id
				sequence = str(record.seq)
				outfile.write( ">" + header + "_f\n" + sequence + "\n")
		for record in SeqIO.parse(fasta_p2, 'fasta'): # open pair one fasta file
			if record.id in headers: # if the header from the fasta file is in the dictionary observed do the following:
				header = record.id
				sequence = str(record.seq)
				outfile.write( ">" + header + "_r\n" + sequence + "\n")	
		outfile.close()
sys.exit()		
