#log in UPPMAX server
ssh -AX -o ServerAliveInterval=10 xung0001@rackham.uppmax.uu.se 
# put password
# go to the right directory
cd /crex/proj/snic2021-23-14/Xue
# create script
cd script
# Step1 : get raw reads from GEO dataset
nano rawread.sh
#from row 11 to row 25 are the contents of bash script
#!/bin/bash -l
#SBATCH -A snic2022-22-89
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH -J radreadsget
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
cd /crex/proj/snic2021-23-14/Xue/raw 
module load bioinfo-tools sratools/2.9.6-1 
for i in $(seq 08 09); 
do fastq-dump --gzip --split-3  SRR53676$i 
done 
# submit the tast to the server
# pre-submit and check if the script has problem
chmod +x rawread.sh
bash rawread.sh 
# run script
sbatch rawread.sh 
# check if the job is submited to uppmax and is in the waitinglist
jobinfo -u xung0001 

# Step2 (optionally) run fastqc
# write script
nano fastqc.sh
# from row 39 to row 50
#!/bin/bash -l
#SBATCH -A snic2022-22-89
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH -J fastqc
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools FastQC/0.11.9 MultiQC
mkdir fastqc 
fastqc -o   
/crex/proj/snic2021-23-14/Xue/tissue_RNA_seq/fastqc   /crex/proj/snic2021-23-14/Xue/tissue_RNA_seq/raw/*.fastq.gz 
multiqc .
# submit the script
sbatch fastqc.sh 

# step 3 remove rRNA with sortmeRNA
 # index rRNA database (only need to do once if the indexed file is always in the directory)
 Module load bioinfo-tools SortMeRNA/4.3.3 
sortmerna --index 1 --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5.8s-database-id98.fasta \
 --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5s-database-id98.fasta \
 --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-16s-id95.fasta \
  --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-23s-id98.fasta \
  --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-16s-id90.fasta \
  --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-23s-id98.fasta \
  --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-18s-id95.fasta \
  --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-28s-id98.fasta  
# creat bash script for sortmeRNA
nano sortmerna_bash.single.sh  # from row 67 to row 89 is bash script
#!/bin/bash -l 

set -eu 
if [ -d /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna/kvdb ]; then 
 rm -rf /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna/kvdb/* 
Fi 

sample=$(basename ${1/.fastq.gz}) 
echo "sample name is $sample" 
module load bioinfo-tools SortMeRNA/4.3.3 
sortmerna --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5.8s-database-id98.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5s-database-id98.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-16s-id95.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-23s-id98.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-16s-id90.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-23s-id98.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-18s-id95.fasta \ 
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-28s-id98.fasta \ 
--reads $1 \ 
--workdir /crex/proj/snic2021-23-14/Xue/sortmerna_mRNA/$sample --fastx --other  
mv /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/$sample/out/aligned.log /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/$sample/${sample}_aligned.log 
mv /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/$sample/out/other.fq.gz  /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/${sample}_sortmerna.fq.gz
# create a loop for running sortmemrna 
nano sortmerna_single.sh # from row 92 to row 102 is a script
#!/bin/bash -l 
#SBATCH -A snic2022-22-89 
#SBATCH -p core 
#SBATCH -n 5 
#SBATCH -t 12:00:00 
#SBATCH -J sortmerna 
#SBATCH --mail-user xue.zhang@slu.se 
#SBATCH --mail-type=ALL 
dir /crex/proj/snic2021-23-14/Xue/GSE96945/rawreads_singlemRNA/ 
for i in $(seq 08 09); 
do bash /crex/proj/snic2021-23-14/Xue/scripts/sortmerna_bash_single.sh /crex/proj/snic2021-23-14/Xue/GSE96945/rawreads_singlemRNA/DRR2353$i.fastq.gz 
done 
# submit the script 
sbatch sortmerna_single.sh 

# step4 trimmomatic reads
# requires: Trimmomatic_bash.sh; trimmomatic.sh 
cd /crex/proj/snic2021-23-14/Xue/GSE96945
mkdir sortmerna_mRNA/fastqc-trimmomatic 
mkdir trimmomatic_mRNA 
nano trimmomatic_bash_single.sh
# from row 113 to 122 as one single script
#!/bin/bash -l
sample=$(basename ${1/_sortmerna.fq.gz})
echo "sample name is $sample"
module load bioinfo-tools trimmomatic/0.36 FastQC/0.11.9  MultiQC/1.11
trimmomatic SE -threads 1  $1 /crex/proj/snic2021-23-14/Xue/GSE96945/trimmomatic_mRNA/${sample}_trimmomatic.fq.gz ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50
fastqc -o /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/fastqc-trimmomatic --noextract -t 1 /crex/proj/snic2021-23-14/Xue/GSE96945/trimmomatic_mRNA/${sample}_trimmomatic.fq.gz\
multiqc -f -o  /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/fastqc-trimmomatic/ .

 nano trimmomatic_single.sh # from row 124 to 135 as a script
#!/bin/bash -l
#SBATCH -A snic2022-22-89
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 96:00:00
#SBATCH -J trimmomatic
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
dir /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA
for fwd in $(find /crex/proj/snic2021-23-14/Xue/GSE96945/sortmerna_mRNA/ -name  "*_sortmerna.fq.gz");
do bash /crex/proj/snic2021-23-14/Xue/scripts/trimmomatic_bash_single.sh $fwd
done
# check and submit the work
chmod +x trimmomatic_single.sh
bash trimmomatic_single.sh
# if no problem
sbatch trimmomatic_single.sh
jobinfo -u xung0001
# check the multiqc file 
firefox multiqc_report.html
# or transform to local computer and check carefully 
scp xung0001@rackham.uppmax.uu.se:/crex/proj/snic2021-23-14/Xue/GSE96945/trimmomatic_mRNA/*.html ~/Desktop
# fill password
# found DRR *55 still have adapter, need to trim again 
# I USED trimmglore for 55

# Step 5 align-free method quantification
# Kllistoquant
nano KallistoquantSE1.sh
#!/bin/bash -l

#SBATCH -A snic2022-22-89
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 24:00:00
#SBATCH -J quantification
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools kallisto/0.46.2
index=/crex/proj/snic2021-23-14/Xue/genome/at11_kallisto.idx
cd /crex/proj/snic2021-23-14/Xue/GSE96945/single_kallisto
for f in $(find /crex/proj/snic2021-23-14/Xue/GSE96945/trimmomatic_mRNA/  -name "*_trimmomatic.fq.gz");
 do
fnam=$(basename ${f/_trimmomatic.fq.gz/})
 kallisto quant -i $index \
 -o ./$fnam \
 -b 100  --single -l 180 -s 20 $f
done
bash KallistoquantSE1.sh # it works
sbatch  KallistoquantSE1.sh


