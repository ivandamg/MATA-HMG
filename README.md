# MATA-HMG


# 1. Get all MATA-hmg genes in predicted nu6 

 - Download all sequences from organism of PF00505 HMG-box

 - Discard sequences that are incomplete in aliview

 - Prepare protein file from target genome gff to fasta

 - from gff to fasta bash command

       cat Predicted_prot_hint_Glomus_Nu6_genome_masked.gff | grep "#" | grep -v "CDS" | grep -v "UTR" | grep -v "Evidence" | grep -v "%" | grep -v "end" | grep -v "hint" | grep -v "###" | sed 's/start gene />/' | sed 's/protein sequence //' | tr ' = [' ' ' | tr ']' ' ' | sed 's/   //' | tr '#' ' ' | sed 's/  //' | grep -v 'E:' | grep -v "sequence number" | grep -v "none" | grep -v "command" > Predicted_prot_hint_Glomus_Nu6_genome_masked.faa

       open in text editor and erase first rows and double spaces

- blast all PFAM sequence on organism all protein file.
    
    makeblastdb -dbtype prot -in Predicted_prot_hint_Glomus_Nu6_genome_masked.faa -parse_seqids -out db_genomes/Predicted_n6_db
    
    blastp -db db_genomes/Predicted_n6_db -outfmt 6 -evalue 1e-6 -show_gis -num_threads 30 -out db_genomes/blast_PF00505-HMG-box_Predictedn6.xml -query ~/Documents/Co-inoculation_manuscript/v7_NewPhytol/MATA_HMG_GENOMES/HMG-BOX_n6/PF00505_HMG-box_complete_sequences.fa 
    
- Select only 100 % identity    
      
      cat db_genomes/blast_PF00505-HMG-box_Predictedn6.xml | awk '$3 == 100' > MATA-HMG_Predictedn6.txt

- Discard duplicate matchs on genome by discarding a second match. in excel.
   MATA-HMG_Predictedn6_unique_VF.txt

# 2. Identify gene transcription on those genes. 

- Create list of HMG-box in predicted-n6 (in excel)
   HMG-box_Predictedn6_names_vf.txt
   
- Identify the list of genes on transcriptomic data and plot results
Script RNA_double_inoc_onlyMATA-HMG_values_v7.R 
Extract all counts per sample, select only the HMG-box and plot the ones that display a confident number of reads

