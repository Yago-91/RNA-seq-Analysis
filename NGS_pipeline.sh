#!/bin/bash

#This pipeline was developed by Ivó Hernández Hernández. ivohh91@gmail.com

#----------------------------------------------------------Manage errors----------------------------------------------------------------------# 
#exit when any command fails
set -e

#keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
#echo an error message before exiting
trap 'echo "\"${last_command}\" Script stopped with exit code $?."' EXIT
#---------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------Format stderror---------------------------------------------------------------------#
color()(set -o pipefail;"$@" 2>&1>&3|sed $'s,.*,\e[31m&\e[m,'>&2)3>&1
#Usage: color command-------------------------------------------------------------------------------------------------------------------------#


#-------------------------------------------------------Define help function------------------------------------------------------------------#
function help() {
    cat<<EOF
    Required arguments:

    -g|--GTF            Path to GTF file. This file must be in
                        the same release version as the genome
                        used in aligning step.

    --gff               Path to GFF file. This file is used for
                        the MAJIQ build step.

    -o|--outdir         The directory to where save output.

    -l|--libType        Either "SINGLE" or "PAIRED" library
                        formats are accedpted.

    -p|--Nthreads       Number of parallel threads to use.
                        Default is 6.
    
    -P|--Proyect        Proyect ID: SRP...

    -I|--Index          Path to STAR genomic indexes.
  
    -X|--tx2gene        File with the correspondance of transcripts
                        and genes for tximport step.

    -G|--gene           Name of the gene for  exon-wise analysis.

    -c|--config         Configuration file.

    -R|--read-length    Read length (bp).


    
    Optional arguments:

    --remove (opt)      If used, this option drives removing of 
                        opt-type files, being opt one of these:
                        [sra|fq|sra:fq|bam|all]. Default: none.
      
    -T|--transcripts    Fasta file with transcripts of the 
                        selected species for Salmon step.

    -F|--fasta          Genome fasta file for creation of the
                        transcripts file for Salmon step.

    
EOF
}
#---------------------------------------------------------------------------------------------------------------------------------------------#


#-----------Setting defaul parameters----------#
REMOVE="none"
#----------------------------------------------#

#---------------------------------------------------------------Parsing options---------------------------------------------------------------#
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -g|--GTF)
    GTF="$2"; shift; shift
    ;;
    --gff)
    GFF="$2"; shift; shift
    ;;
    -o|--outdir)
    OUT_dir="$2";shift; shift
    ;;
    -l|--libType)
    LIBTYPE="$2";shift; shift
    ;;
    --remove)
    REMOVE="$2";shift; shift
    ;;
    -p|--Nthreads)
    CORES="$2";shift; shift
    ;;
    -P|--Proyect)
    SRP="$2";shift; shift
    ;;
    -G|--gene)
    GENE="$2";shift; shift
    ;;
    -I|--Index)
    GENOME_INDEX="$2";shift; shift
    ;;
    -T|--transcripts)
    TRANSCRIPTS="$2";shift; shift
    ;;
    -F|--fasta)
    FASTA="$2";shift; shift
    ;;
    -X|--tx2gene)
    TX2GENE="$2";shift; shift
    ;;
    -c|--config)
    CONFIG="$2";shift; shift
    ;;
    -R|--read-length)
    READ_L="$2";shift; shift
    ;;
    -h|--help)
    help
    exit 0
    ;;
    *)    # unknown option
    echo "Unrecognized option"; help; exit 1
    ;;

esac
done
#---------------------------------------------------------------------------------------------------------------------------------------------#

dir_name=$(dirname $0)

#--------------------------------------------------------------Create log file----------------------------------------------------------------#
touch $OUT_dir/run_log.txt
echo "The script has been run with the following command: $0 "$@$'\n' > $OUT_dir/run_log.txt
#---------------------------------------------------------------------------------------------------#

#----------------------------------------------------------Parse configuration file-----------------------------------------------------------#
Accessions=$(cat $CONFIG | cut -d ' ' -f1)
cat $CONFIG | cut -d ' ' -f1 > $OUT_dir/SRA_file.tmp
Groups=$(cat $CONFIG | cut -d ' ' -f3 | uniq)
Group_number=()
Names=()

for i in $Groups; do Group_number+=($i); Group_number+=($(grep -cw $i $CONFIG)); done
#---------------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------Set up global variables------------------------------------------------------------#
#Extended date: $(date +"%x %r %Z")
#Just the time $(date +"%r")

declare -A genomes=(["hsapiens"]="hg38" ["mmusculus"]="mm10" ["rnorvegicus"]="rn6")

n_samples=$(cat $OUT_dir/SRA_file.tmp | wc -l)
to_species=$(head -1 $OUT_dir/SRA_file.tmp)
Species=$(esearch -db sra -query $to_species | efetch -format xml | xtract -pattern EXPERIMENT_PACKAGE -block SAMPLE -element SCIENTIFIC_NAME | cut -d ' ' -f1)
cat $GTF | grep -w $GENE | head -1 > $OUT_dir/Info_line.txt #Line with the gene information.
Gene_start=$(cat Info_line.txt | cut -d $'\t' -f4)
Gene_stop=$(cat Info_line.txt | cut -d $'\t' -f5)
Chrom=$(cat Info_line.txt | cut -d $'\t' -f1)
Genome="${genomes[$Species]}"
Quant_dirs=() #To save output directories names of salmon quant step.

echo "Job start: $(date +"%x %r %Z")"$'\n'
echo "Job start: $(date +"%x %r %Z")"$'\n' >> $OUT_dir/run_log.txt
echo "Info_line generated $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
cat $OUT_dir/Info_line.txt >> $OUT_dir/run_log.txt
#---------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------Create directories and files----------------------------------------------------------#
mkdir -p $OUT_dir/SRA
mkdir -p $OUT_dir/Fastq
mkdir -p $OUT_dir/BAM
mkdir -p $OUT_dir/Cropped_bams
mkdir -p $OUT_dir/Exon_quant
mkdir -p $OUT_dir/Salmon
mkdir -p $OUT_dir/MAJIQ

echo "Creating Exons table...$(date +"%r")"$'\n'
echo $'\n'"Creating Exons table... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt

awk -F $'\t' '$3 ~ /exon/' $GTF | grep -w $GENE | awk -F $'\t' '{print $4,$5,$9}'   | awk -F ';' '{print $1,$3,$5,$12}' \
| awk -F ' ' '{print $4,$6,$10,$7":"$8,$1"-"$2}' | sed 's/\"//g' | grep -v CCDS | sort -u -t ' ' -k5 > $OUT_dir/Exon_quant/TAF1-exons.table
#---------------------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------Prepare SRA acc=url---------------------------------------------------------------#
if [[ -f "$OUT_dir/${SRP}_metadata.tab" ]]; then
    echo "Found metadata for project $SRP"$'\n'
    echo "Found metadata for project $SRP"$'\n' >> $OUT_dir/run_log.txt

else
    echo "Retrieving metadata for project $SRP"$'\n'
    echo "Retrieving metadata for project $SRP"$'\n' >> $OUT_dir/run_log.txt
    pysradb metadata $SRP --detailed --saveto $OUT_dir/${SRP}_metadata.tab

fi

#Here we extract the accession number and url for each sample
#We use the awk script "get_colums"
awk -F $'\t' -f $dir_name/get_columns -v cols=run_accession,sra_url $OUT_dir/${SRP}_metadata.tab > $OUT_dir/sra_url.tmp
#Extract the SRR accessions of interest with the given file SRA_file
grep -Fwf $OUT_dir/SRA_file.tmp $OUT_dir/sra_url.tmp > $OUT_dir/acc_url.txt
cat $OUT_dir/acc_url.txt | cut -d ' ' -f 2 > $OUT_dir/url_list.tmp
rm -fr $OUT_dir/sra_url.tmp
#---------------------------------------------------------------------------------------------------------------------------------------------#


#1) In this step all the provided SRR accession are downloaded and dumped one by one to avoid sotrage problems (once dumped .sra files are erased)
# Inside the while loop samples are also aligned with STAR (needed for splicing analysis) using selected index (total or partial genome --> seed up)
# Then .fastq samples are compressed in parallel processing with pigz as the next step is able to use .gz files
# At the last step .fastq.gz files are provided to Salmon (in this case all the transcriptome is used to generate an index)

   #1.1----------------------------------------------------Download complete project---------------------------------------------------------#
   echo $'\n'"Downloading sra files... $(date +"%r")"$'\n'
   echo $'\n'"Downloading sra files...$(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
   while read url
        do
            wget -N --quiet -P $OUT_dir/SRA/ $url &
        done < url_list.tmp
   wait

   echo "Download completed  $(date +"%r")"$'\n'
   echo "Download completed  $(date +"%r")" >> $OUT_dir/run_log.txt
   rm -fr url_list.tmp
   rename 's/\..+//' $OUT_dir/SRA/SRR*
   #----------------------------------------------------------------------------------------------------------------------------------------#

#From this point we use perform the steps one file at a time to avoid storage problems as the files to be generated can be very big
for Acc in ${Accessions//,/ }
  do

    Name+=($(cat $CONFIG | grep -w $Acc | cut -d ' ' -f2))
    #1.2-------------------------------------Dumping SRR files with fasterq-dump and quality control------------------------------------------#
    echo "Processing sample ${Acc}...$(date +"%r")"
    if [[ $LIBTYPE = "SINGLE" ]]; then
        echo $'\n'"Dumping file $Acc... $(date +"%r")" >> $OUT_dir/run_log.txt
        fasterq-dump $OUT_dir/SRA/$Acc -e $CORES -t $OUT_dir/Fastq/tmp/Fasterq_dump -O $OUT_dir/Fastq/ --split-3 2>>$OUT_dir/run_log.txt
        
    else
        echo "Dumping file $Acc...$(date +"%r")" >> $OUT_dir/run_log.txt
        fasterq-dump $OUT_dir/SRA/$Acc --split-3 -e $CORES -t /tmp/Fasterq_dump -O $OUT_dir/Fastq/ 2>>$OUT_dir/run_log.txt                                                                                                                                                                            

    fi

    fastqc -o $OUT_dir/Fastq/ $OUT_dir/Fastq/${Acc}.fastq 
    #-----------------------------------------------------------------------------------------------------------------------------------------#

    #1.3---------------------------------------------------Removing .sra file-----------------------------------------------------------------#
    if [[ $REMOVE = "all" ]] || [[ $REMOVE = "sra" ]] || [[ $REMOVE = "sra:fq" ]]; then
        echo $'\n'"Removing ${Acc} sra file... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
        rm -fr $OUT_dir/SRA/*
    fi
    #-----------------------------------------------------------------------------------------------------------------------------------------#


    #1.4-------------------------------------------STAR 2.7 align to selected reference-------------------------------------------------------#
    #We are going to take advantage of the cache: in we wanted to configure this we just execute ./vdb-config -i
    echo $'\n'"Started alignment of sample $Acc with STAR... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt

    if [[ $LIBTYPE = "SINGLE" ]]; then
        STAR --runThreadN $CORES --genomeDir $GENOME_INDEX --readFilesIn $OUT_dir/Fastq/${Acc}.fastq --outSAMtype BAM Unsorted \
        --sjdbGTFfile $GTF --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --outFileNamePrefix $OUT_dir/BAM/${Acc}_ \
        > $OUT_dir/run_log.txt

    else
        STAR --runThreadN $CORES --genomeDir $GENOME_INDEX --readFilesIn $OUT_dir/Fastq/${Acc}_1.fastq $OUT_dir/Fastq/${Acc}_2.fastq --outSAMtype BAM Unsorted \
        --sjdbGTFfile $GTF --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend --outFileNamePrefix $OUT_dir/BAM/${Acc}_ \
        >$OUT_dir/run_log.txt
    fi
    #-----------------------------------------------------------------------------------------------------------------------------------------#


    #1.5--------------------------------------------Remove fastq files or compress them-------------------------------------------------------#
    if [[ $REMOVE = "all" ]] || [[ $REMOVE = "fq" ]] || [[ $REMOVE = "sra:fq" ]]; then
        rm -fr $OUT_dir/Fastq/*.fastq

        else 
            echo "Compressing fastq file/s of sample $Acc... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
            pigz -7 -p $CORES $OUT_dir/Fastq/*.fastq

    fi  
    #-----------------------------------------------------------------------------------------------------------------------------------------#

    #1.6------------------------------------------------Sort, index and crop bam files--------------------------------------------------------#
    echo "Started sorting of sample $Acc... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt

    color samtools sort -@ $CORES -o $OUT_dir/BAM/${Acc}.bam $OUT_dir/BAM/${Acc}_Aligned.out.bam 2>>$OUT_dir/run_log.txt
    color samtools index -@ $CORES $OUT_dir/BAM/${Acc}.bam $OUT_dir/BAM/${Acc}.bam.bai

    echo $'\n'"Preparing cropped bam file for sample ${Acc}...$(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
    samtools view -@ $CORES $OUT_dir/BAM/${Acc}.bam "$Chrom:$Gene_start-$Gene_stop" > $OUT_dir/Cropped_bams/${Acc}.cropped.bam
    
    #Remove unsorted bam
    rm -fr $OUT_dir/BAM/${Acc}_Aligned.out.bam
    #-------------------------------------------------------------------------------------------------------------------------------------------#

    #1.6----------------------------------------------------Exon quantification-----------------------------------------------------------------#
    echo "Started exon quantification of sample $Acc... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt

    while read exon
	  do
		  coord=$(echo $exon | cut -d ' ' -f 5)

		  samtools view -@ $CORES $OUT_dir/BAM/$Acc.bam "$Chrom:$coord" | wc -l >> $OUT_dir/Exon_quant/$Acc.counts.tmp

	  done < $OUT_dir/Exon_quant/TAF1-exons.table
    
    sed -i "1s/^/${Acc}\n/" $OUT_dir/Exon_quant/${Acc}.counts.tmp

    if [[ $REMOVE = "all" ]] || [[ $REMOVE = "bam" ]]; then
        rm -fr $OUT_dir/BAM/*
    fi
    #-------------------------------------------------------------------------------------------------------------------------------------------#


    #1.9------------------------------------Run Salmon with aligning mode and erase fastq.gz files----------------------------------------------#
    echo "Started Salmon quantification of sample $Acc... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
    if [[ ! $TRANSCRIPTS ]]; then 
        echo "Building transcripts.fa... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt
        gffread -w $OUT_dir/Salmon/transcripts.fa -g $FASTA $GTF 2>/dev/null
        TRANSCRIPTS="$OUT_dir/Salmon/transcripts.fa"
    fi

    salmon quant -t $TRANSCRIPTS -l A -a $OUT_dir/BAM/${Acc}_Aligned.toTranscriptome.out.bam -o $OUT_dir/Salmon/Salmon_quant_${Acc} 2>>$OUT_dir/run_log.txt
    
    Quant_dirs+=(Salmon/Salmon_quant_${Acc})    
    #-------------------------------------------------------------------------------------------------------------------------------------------#

    echo "Finished processing of sample ${Acc}...$(date +"%r")"$'\n'
    echo $'\n'"Finished processing of sample ${Acc}...$(date +"%r")"$'\n' >> $OUT_dir/run_log.txt

done

#1.9.1---------------------------------------Mergin output of salmon quant and pass to tximport-Deseq2---------------------------------------------# 
#There is not a good correspondance between GTF and cdna.fa files regarding tx_ID-Gene_ID so -g option in salmon quant is of no use
echo $'\n'"Performing differential expression analysis... $(date +"%r")"$'\n'
echo "Performing differential expression analysis... $(date +"%r")"$'\n' >> $OUT_dir/run_log.txt

salmon quantmerge --quants ${Quant_dirs[@]} --names $Names -o $OUT_dir/Salmon/Salmon_tx_merge.txt 2>>$OUT_dir/run_log.txt

Dirs=$(printf "%s," "${Quant_dirs[@]}" | cut -d "," -f 1-${#Quant_dirs[@]})
Titles=$(printf "%s," "${Names[@]}" | cut -d "," -f 1-${#Names[@]})
Groups_R=$(printf "%s," "${Group_number[@]}" | cut -d "," -f 1-${#Group_number[@]})

Rscript --vanilla $dir_name/Salmon_analysis.R $OUT_dir $Dirs $Species $Titles $TX2GENE $Groups_R 2>>$OUT_dir/run.log.txt
#-----------------------------------------------------------------------------------------------------------------------------------------------#



#2) Quantification by exon: reads per exon normalized by the total exonic count in each sample.
    
    echo "Performing exon-wise analysis...$(date +"%r")"$'\n'
    #2.1--------------------Set header of exons file----------------------------#
    sed -i "1s/^/Gene_ID Transcript_ID Exon_ID Exon_number Coordinates\n/" $OUT_dir/Exon_quant/TAF1-exons.table

    #2.2--------------------Join .counts.tmp files and exons file---------------#
    paste -d ' ' $OUT_dir/Exon_quant/TAF1-exons.table $OUT_dir/Exon_quant/*counts.tmp > $OUT_dir/Exon_quant/merge_exon_quant.table
    rm -fr $OUT_dir/Exon_quant/*counts.tmp

    #2.3-------------------------Run Rscript------------------------------------#
    Dir_ex="$OUT_dir/Exon_quant"

    Rscript --vanilla $dir_name/Exon_analysis.R $Dir_ex $n_samples $Names  2>>$OUT_dir/run.log.txt


#3) Splicing analysis step: in this step we call the software MAJIQ (Biochiphers) for which we have to pass a configuration file
#that is made on the fly

    echo "Started splicing analysis...$(date +"%r")"$'\n'
    #3.1 Building configuration file
    touch config.txt
    echo "[Info]"$'\n' > config.txt
    echo "readlen=$READ_L"$'\n' >> config.txt
    echo "genome=$Genome"$'\n' >> config.txt
    echo "bamdirs=$OUT_dir/BAM/"$'\n' >> config.txt
    echo "strandness=none"$'\n' >> config.txt
    echo "[experiments]"$'\n' >> config.txt
    echo ${Group_number[0]}=$(cat SRA_file.tmp | head -${Group_number[1]} | tr $'\n' ',') >> config.txt
    echo ${Group_number[2]}=$(cat SRA_file.tmp | head -${Gorup_number[3]} | tr $'\n' ',') >> config.txt

    mv config.txt $OUT_dir/MAJIQ/

    conda activate MAJIQ

    #3.2 Running MAJIQ build step
    majiq build $GFF -c $OUT_dir/MAJIQ/config.txt -j $CORES --min-experiments 1 --min-denovo 1 --minpos 1 --minreads 1 -o $OUT_dir/MAJIQ/ 2>/dev/null
    
    #3.3 Running MAJIQ deltapsi step
    majiq deltapsi -grp1 $OUT_dir/MAJIQ/${Group_number[0]}_*.majiq -grp2 ${Group_number[2]}_*.majiq -j $CORES -o DPSI \
    --min-experiments 1 --minreads 1 --minpos 1 --output-type voila --names ${Group_number[0]} ${Group_number[2]} 2>/dev/null

    #3.4 Running voila tsv step
    voila tsv $OUT_dir/MAJIQ/splicegraph.sql $OUT_dir/MAJIQ/DPSI/*.deltapsi.voila --show-all -j $CORES -f ${Group_number[0]}_vs_${Group_number[2]}.tsv 2>/dev/null

    rm -fr SRA_file.tmp


#4) Cleaning step

if [[ $REMOVE = "bam" ]] || [[ $REMOVE = "all" ]]; then
    rm -fr $OUT_dir/BAM
fi

rm -fr $OUT_dir/Info_line.txt

echo "Analysis completed successfully $(date +"%r")"

#####End of the journey#####




