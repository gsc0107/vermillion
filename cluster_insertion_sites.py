#! /usr/bin/env python
import sys

#script takes a blast file and clusters the hits based on start location in the database genome
#Note: only keeps clusters with more than one member
#usage
#python cluster_formation.py blastFile outputFileName

#get the blast output file and the name for the output file from this script
try:
	gallus_hits=sys.argv[1]
except IndexError:
	print "\n blast output file not supplied."
	sys.exit()
except IOError:
	print "\n blast output file not found in directory."
	sys.exit()
try:
	outfile=sys.argv[2]
except IndexError:
	print "\n output filename file not supplied."
	sys.exit()

#cluster the hits based on where they are in the genome.
cluster = {} # make an empty dictionary clusters
file = open(gallus_hits, "r")
for line in file:
	line=line.rstrip()
	line=line.split("\t")
	found=0
	for key in cluster:
		if abs(int(line[8])-int(key))<=1000: #first check and see if the start position is within 1000 nt of the key. The key is the start position it is a string for a mathematical formula you need to use integers only
			if line[1]==cluster[key][0]: # if the start is within 1000 check and see if the chr is the same. 
				cluster[key][1].append(line[0])#append the header to the dictionary
				found=1 # if key is found 
				break # stop the loop
	if found==0: #if the key has not been found then found is still set to zero
		cluster[line[8]]=[line[1],[line[0]]] # make a new cluster with the key in the start and populate it with chr and header
file.close()

for k, v in cluster.items():
	if len(v[1]) == 1:
		del cluster[k]

#output the clusters to file
outfile=open(outfile, 'w')		
for key, value in cluster.items():
	#key is the start number so write this to file
	outfile.write(key + '\t')
	#value has the chromosome number in element 0 and a list in element 1
	#write the chromosome number out
	outfile.write(value[0] + '\t')
	#for each header in the list (value[1]) write to file with a space in between
	outfile.write(' '.join(value[1]))
	outfile.write('\n')
outfile.close()	
sys.exit()