@echo off

echo Make database
makeblastdb -in chr5_hg13.fa -dbtype nucl -out hg13
makeblastdb -in chr5_hg15.fa -dbtype nucl -out hg15

pause