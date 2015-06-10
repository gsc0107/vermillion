# vermillion
A bioinformatics package for the processing of targeted DNA sequencing to discover novel endogenous retroviruses.
#Sequencing files are received in fastq format from genome quebec. They do the QC (I think trimming to phred score of 30)

#unzip files gunzip filename
#cat files if from more than one lane 

#convert from fastq to fasta
python
from Bio import SeqIO
SeqIO.convert('filename', "fastq", 'outfile name', "fasta")

#make a BLAST db for ev-1
makeblastdb -in file.fasta -title name of db -dbtype nucl -out dbname -parse_seqids
#BLAST pair one against ev-1
#BLAST pair two against ev-1


#note use culling_limit 1 to get only one hit or else you get hits to both LTRs
#agrument -num_threads didn't work on moa, should ask Karl if resolved
#outfmt 6 = tab delimited

blastn -query -db ALVE/ALVE1 -culling_limit 1 -out -outfmt 6

#Blast output looks like this:
HWI-ST913:232:C431UACXX:4:1101:1882:2212        gi|13508434|gb|AY013303.1|	100.00  100     0	0	1	100     2454    2355    8e-51    185
HWI-ST913:232:C431UACXX:4:1101:4172:2094        gi|13508434|gb|AY013303.1|	97.59   83	2	0	1	83	6058    6140    3e-39    147
HWI-ST913:232:C431UACXX:4:1101:4964:2087        gi|13508434|gb|AY013303.1|	96.25   80	3	0	1	80	275     196     6e-36    135



#take the blast output and capture the information on what part of the sequences match ALVE-1 - the goal is to trim these sequences to remove ALVE portion and store trimmed reads in a new file
#need to capture header and start and stop coordinates, for both paired reads, of where the sequence matches ALVE-1

#use the python script trimInternalGenome.py

python trimInternalGenome.py blast_output_pair_1.txt blast_output_for_pair_2.txt fasta_file_for_pair_1.fasta fasta_file_for_pair_2.fasta outfile_name.fasta

#the fasta file generated from this script has chicken only sequence - the ALVE sequences have been trimmed.
#if the chicken portion of the sequence was less than 20 nt it was deleted.

#next step is to BLAST new, trimmed fasta file to chicken genome
#make a BLAST database for chicken genome - needed to reformat chrZ and chrW - seemed to be problems mixing integers and letters later on. 
#replaced chrZ = chr100, chrW =chr200 using sed command
sed 's/chrZ/chr100/g' galgal4.fa > galgal4_edit.fa
#also removed the unknown chr from the database this was done by searching for headers with chr\d+, one or more digit
#I think I used Morgan's code to do this!

makeblastdb -in file.fasta -tile -dbtype nucl -out db name -parse_seqids


blastn -query -db gallus/gallus -culling_limit 3 -evalue 1e-30 -out -outfmt 6

#blast output file looks like this:
HWI-ST913:232:C431UACXX:4:2113:6249:9539_2	chr100  100.00  100     0	0	1	100     10540416        10540515        9e-46    185
HWI-ST913:232:C431UACXX:5:2108:13539:96064_1    chr9    100.00  100     0	0	1	100     11332287        11332188        9e-46    185
HWI-ST913:232:C431UACXX:5:2115:12140:44310_2    chr100  100.00  96	0	0	1	96	10540446        10540541        1e-43    178

#there are some formatting issues with this when doing the next step, which is reading the file and capturing the chr and start position and reading into a dictionary
#three modifications - remove "_1" and "_2" in the header name, replace with nothing

sed 's/_1//g' file_name.txt > new_file_name.txt 

#replace "chr1" with "1" - problems in next step because not an integer
sed 's/chr//g' file_name.txt > new_file_name.txt 

#next step is to cluster the reads which are aligning to the same region of the genome - all these hits are providing information for the same insertion
#python script to do this is cluster_insertion_sites.py

python cluster_insertion_sites.py blast_output_gallus_hits.txt outfile.txt

#sample of the outfile is below start /t chr /t headers sep by space
16717828        12	HWI-ST748:216:C1NV4ACXX:1:2305:18740:139030 HWI-ST748:216:C1NV4ACXX:1:1104:18030:136722
52624682        5	HWI-ST748:216:C1NV4ACXX:1:2201:12187:72988 HWI-ST748:216:C1NV4ACXX:1:2108:18063:104650


#now take this file and collect all the sequences (both pairs) for each header in the cluster. Write each cluster to a new file so you can look for insertion sites

python cluster_sequence_files.py cluster_file.txt file_name_information original_fasta_file_pair_1.fasta original_fasta_file_pair_2.fasta

#now each cluster has its own file COTW3(string)_chr_start.fasta that when you open you can look for junctions

#to find the sequences which have LTRs - to make the junction search easier can do a substring search

for i in *.fasta; do python flitering_3LTR.py $i"_3LTR.txt" $i; done;
#finds all the files with that substring and writes them to a new file
# will write the new file even if there is no information
#find files with zero bytes and move
#it would be better to do a search which allows for mis-matches - I think Scott has done this!
mkdir files_zero
find . -type f -size 0 -exec mv {} files_zero/ \; 

#what I could do: BLAST e-value cut-offs rather than culling_limit of three, I think I should leave this out the filtering process picks up on this later
#call all thses commands in one file









