
#!/bin/bash -l
cd /crex/proj/snic2021-23-14/Xue/PRJDB10113/single_kallisto
for i in $(seq 54 64);
do
cat ./DRR2353$i/abundance.tsv |awk -F "\t" '{print $5}' >DRR2353$i.tsv
ls -v DRR2353$i.tsv
done
paste ./DRR235353/abundance.tsv *.tsv >>total.tsv

# submit the bash file: bash renamefilenames.sh