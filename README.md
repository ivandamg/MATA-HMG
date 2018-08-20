# MATA-HMG



1. Find similar genes in genome. 

# get query sequence

# get target genome assemblies
https://www.ncbi.nlm.nih.gov/assembly/

 # make db
 a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f1) ;makeblastdb -dbtype nucl -in $i -parse_seqids -out db_genomes/$(echo $i | cut -d'_' -f1)_db ; done

 # blast sequences
 a=0;for i in $(ls *.fna); do echo $(echo $i | cut -d'_' -f1) ;tblastn -db db_genomes/$(echo $i | cut -d'_' -f1)_db -outfmt 6 -evalue 1e-6 -show_gis -num_threads 30 -out db_genomes/blast_hmg127_$(echo $i | cut -d'_' -f1).xml -query ~/Documents/Co-inoculation_manuscript/v7_NewPhytol/MATA_HMG_GENOMES/Mata-hmg127_n6_raw.fasta ; done

  # Make summary table of blast hits. evalue > 1.e-06 > length min 90 %
  
2. 
