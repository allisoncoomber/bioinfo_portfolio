#!/bin/bash

fastQCPATH="/usr/local/usrapps/jbr/FastQC"
export PATH="$PATH:$fastQCPATH"

for read1 in $PATH
do	
	read2=$(echo $read1| sed 's/sra_1.fastq/sra_2.fastq/');
	/usr/local/usrapps/jbr/TrimGalore-0.6.5/trim_galore -q 20 --fastqc --gzip --paired --trim-n --length 20 -o /rsstu/users/j/jbr/Phytophthora-database/modern --path_to_cutadapt /usr/local/usrapps/jbr/cutadapt $read1 $read2
done
