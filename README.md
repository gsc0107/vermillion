#Vermillion
A bioinformatics package for the processing of targeted DNA sequencing to discover novel endogenous retrovirus insertion sites.

#Workflow

1. Blast targetted sequences against your transposable element of interest

    * Make a BLAST database for your transposable element of interest (e.g. ALVE-1)

            makeblastdb -in transposable_element.fasta -dbtype nucl -out transposable_element_db -parse_seqids

    * Notes: 
        * Ensure your input targetted sequencing files are in fasta format (not fastq)
        * Blast DB Options: culling_limit 1 (to keep top hit only) outfmt 6 (tab delimited output)

    * BLAST against pair 1 

            blastn -query -db transposable_element_db -culling_limit 1 -out blast_output_pair_1.txt -outfmt 6

    * BLAST against pair 2

            blastn -query -db transposable_element_db -culling_limit 1 -out blast_output_pair_2.txt -outfmt 6

2. Filter for informative sequences and trim away transposable element (leaving only genome sequence > 20 nt)

        python trimInternalGenome.py blast_output_pair_1.txt blast_output_for_pair_2.txt targetted_seqs_pair_1.fasta targetted_seqs_pair_2.fasta informative_seqs.fasta

3. Blast informative and trimmed sequences against the host genome

        makeblastdb -in genome.fasta -dbtype nucl -out genome_db name -parse_seqids

        blastn -query informative_seqs.fasta -db genome_db -evalue 1e-30 -out genome_hits_raw.txt -outfmt 6

    * Do some slight formatting to the output file before clustering

            sed 's/_1//g;s/chr//g' genome_hits_raw.txt > genome_hits.txt 

4. Cluster the reads which are aligning to the same region of the genome 

        python cluster_insertion_sites.py genome_hits.txt clusters.txt

5. Output original sequences for each cluster to identify genome insertions sites and for possible PCR primer design

        python cluster_sequence_files.py clusters.txt file_prefix targetted_seqs_pair_1.fasta targetted_seqs_pair_2.fasta

