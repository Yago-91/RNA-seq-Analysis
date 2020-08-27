#!/bin/bash

#This script was written by Ivó Hernández Hernández (ivohh91@gmail.com)

#HISAT2 indexes were downloaded from https://cloud.biohpc.swmed.edu/index.php/s/grcm38/download
#STAR indexes were generated with the following command: STAR --runThreadN 12 --runMode genomeGenerate --genomeFastaFiles /path/to/Mus_musculus.GRCm38.dna.primary_assembly.fa --genomeDir STAR_ms_2.7 --sjdbGTFfile /path/to/Mus_musculus.GRCm38.100.gtf
#BarraCUDA indexes were generated with the following command:  barracuda index -a bwtsw -p GRCm38_barracuda /path/to/Mus_musculus.GRCm38.dna.primary_assembly.fa
#Whippet indexes were generated with the following command: julia /path/to/whippet-index.jl --fasta /path/to/Mus_musculus.GRCm38.dna.primary_assembly.fa --bam merged.sorted.rmdup.bam --gtf /path/to/Mus_musculus.GRCm38.100.gtf --bam-both-novel -x Whippet_index
#The file merged.sorted.rmdup.bam was generated using samtools merge with all the files in the experiment SRP116945.

#------------------------HISAT2------------------------------------#
echo "Hisat loop started" $(date +"%r") > alignment_trials_log.txt

for file in *.fastq
	do
		hisat2 -x /path/to/HISAT/genome -p 12 --qc-filter -U $file --remove-chrname 2>>alignment_trials_log.txt | samtools view -bS -> ${file}.bam
	done

echo "Hisat loop finished" $(date +"%r") $'\n'>> alignment_trials_log.txt
#------------------------------------------------------------------#

#--------------------------STAR------------------------------------#
echo "STAR loop started" $(date +"%r") >> alignment_trials_log.txt

for file in *.fastq
	do
		STAR --runThreadN 12 --genomeDir /path/to/Mouse/STAR_ms_index_2.7 --readFilesIn $file --outSAMtype BAM Unsorted \
                --sjdbGTFfile /path/to/Mus_musculus.GRCm38.100.gtf --genomeLoad LoadAndKeep --outFileNamePrefix ${file}_
	done

echo "STAR loop finished" $(date +"%r") $'\n'>> alignment_trials_log.txt
#------------------------------------------------------------------#

#--------------------------BarraCUDA-------------------------------#
barracuda index -a bwtsw -p GRCm38_barracuda /path/to/Mus_musculus.GRCm38.dna.primary_assembly.fa

echo "BarraCUDA loop started" $(date +"%r") >> alignment_trials_log.txt

for file in *.fastq
	do
		barracuda aln B_CUDA_index/GRCm38_barracuda $file > ${file}.sai
		barracuda samse B_CUDA_index/GRCm38_barracuda ${file}.sai $file | samtools view -bS -> ${file}.bam
	done

echo "BarraCUDA loop finished" $(date +"%r")$'\n' >> alignment_trials_log.txt
#------------------------------------------------------------------#

#--------------------------Vast-tools------------------------------#
echo "Vast-tools started" $(date +"%r") >> alignment_trials_log.txt

for file in *.fastq
	do
		~/bin/vast-tools align $file --sp mm10 -c 12
	done

echo "Vast-tools finished" $(date +"%r")$'\n' >> alignment_trials_log.txt
#------------------------------------------------------------------#

#--------------------------Whippet------------------------------#
echo "Whippet started" $(date +"%r") >> alignment_trials_log.txt

for file in *.fastq
	do
		julia ~/path/to/whippet-quant.jl $file -o $file -x Whippet_index.jls --sam > ${file}.sam
		samtools view -bS ${file}.sam -> ${file}.bam
	done

echo "Whippet finished" $(date +"%r")$'\n' >> alignment_trials_log.txt
#------------------------------------------------------------------#
