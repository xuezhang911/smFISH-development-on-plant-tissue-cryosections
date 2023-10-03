# build up conda environment in terminal # activate conda environment
# download package blast in conda environment
Conda install -c bioconda blast
# install perl
conda install perl-digest-md5
# download wget function
# set up brew tool https://code2care.org/howto/install-homebrew-brew-on-m1-mac
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
Brew install wget

# OPTION1: blast probe sequences with genomic fasta file and annotate genes with annotation file 
## STEP1 : download arabidopsis genomic fasta file or another species such as barley: wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 
gzip -d Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz # gunzip Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa.gz 
### make a indexed blast database for Arabidopsis 
Makeblastdb -in Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -dbtype nucl -out TAIR10 -parse_seqids # change to barley:Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa , change Tair10 to Hv
### modify the path of blast 
Touch .ncbirc
BLASTDB=/User/xue/Database/blast
BLASTDB_NUCL_DATA_LOADER=/User/xue/Database/blast/nt
### use update_blastdb.pl to update and download nt file 
nohup time update_blastdb.pl nt nr > log &
## Step2 creat a local fasta file with smFISH probe sequences designed mannually or on stellaris website 
Nano q.fa
## Step3 blast probe sequences with indexed databast option1: Blastn -task blastn -db TAIR10 -query q.fa -outfmt 7 -out q.txt chose the output format you prefer 
Blastn -task blastn -db TAIR10 -query q.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore" -out q.txt
wc q.txt
less q.txt
### filter q.txt with the overlap length greater than 16bp
awk -F "\t" '$4-$5-$6 >16' q.txt>nq.txt
### prepare arabidopsis gtf file or other species: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/gff3/hordeum_vulgare/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.55.chr.gff3.gz
wget ftp://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.54.gtf.gz
mv Arabidopsis_thaliana.TAIR10.54.gtf ath.gtf
### extract the chromosome sequence
 sed 's/"/\t/g' ath.gtf |awk -v type="gene" 'BEGIN{OFS=FS="\t"}{if($3=="gene") {print $1,$4-1,$5,$10,".",$7}}' > ath.geneid.bed #for barley is intead $9
### install bedtools
 conda install -c bioconda bedtools
### prepared bed file 
Mv nq.txt nq.bed
cat nq.bed |awk -F "\t" '{if($9<$10) {print $2,$9-1,$10}}'>1.bed
cat nq.bed |awk -F "\t" '{if($9>$10) {print $2,$10-1,$9}}'>2.bed
cat 1.bed >>2.bed
### check if the file is tab delimited 
 cat -t 2.bed
### make it as a tab delimited file 
perl -p -i -e 's/ /\t/g' 2.bed
## Step4:  generate a file including annotated gene name 
bedtools intersect -a 2.bed -b ath.geneid.bed -wa -wb >a.bed 

#option2 blast smFISH probes with cDNA fasta file
## Step1 : prepare file for cDNA fasta 
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
gzip -d Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
mv Arabidopsis_thaliana.TAIR10.cdna.all.fa  Ath.cDNA.fa
### index file
Makeblastdb -in Ath.cDNA.fa -dbtype nucl -out Ath -parse_seqids
## Step2: prepare a fasta file including smFISH probe sequences
nano q.fa
## Step3 blast probe sequences with indexed cDNA databast 
Blastn -task blastn -db Ath -query q.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore" -out q.txt
less q.txt
wc q.txt


 





