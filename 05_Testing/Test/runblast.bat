@echo off

echo Running BLAST search...
blastn -query c6.fasta -db hg13 -out result_hg13.txt -evalue 1e-5 -outfmt 7
blastn -query c6.fasta -db hg15 -out result_hg15.txt -evalue 1e-5 -outfmt 7