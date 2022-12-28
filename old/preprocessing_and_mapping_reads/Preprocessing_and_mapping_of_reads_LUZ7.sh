#!/bin/bash
#needed in batch file to get parameters
while getopts 'i:f:' flag; do
  case "${flag}" in
      i) inputDirectory="${OPTARG}";;
      f) inputDirectoryFastaFiles=${OPTARG} ;;
      ?) echo "I don't know what flag is this" ;;
           esac
       done

#initialization of folder names etc.
output_folder_pychopper=output_pychopper
output_alignment=BAM_files

mkdir "$output_folder_pychopper"
for file in "$inputDirectory"/*.fastq
do
    filename_extended=${file##*/}
    filename=${filename_extended%.*}
    mkdir "$output_folder_pychopper"/$filename
    cdna_classifier.py \
        -r "$output_folder_pychopper"/"$filename"/${filename}_report.pdf \
        -u "$output_folder_pychopper"/"$filename"/${filename}_unclassified.fq \
        -S "$output_folder_pychopper"/"$filename"/${filename}_cdna_classifier_report.tsv \
        "$file" \
        "$output_folder_pychopper"/"$filename"/${filename}_full_length_output.fq
done

mkdir temp
mkdir temp/cutadapt
for file in "$output_folder_pychopper"/*/*_full_length_output.fq
do
    filename_extended=${file##*/}
    searchstring="_full_length_output.fq"
    filename=${filename_extended%$searchstring*}
    ~/.local/bin/cutadapt -a A{10} -e 1 -j 0 -o temp/cutadapt/${filename}_full_length_output_cutadapt.fq ${file}
    ~/.local/bin/cutadapt -g TTTCTGTTGGTGCTGATATTGCTGGG -e 1 -j 0 -o temp/cutadapt/${filename}_full_length_output_cutadapt_2.fq temp/cutadapt/${filename}_full_length_output_cutadapt.fq
    rm temp/cutadapt/${filename}_full_length_output_cutadapt.fq
done

mkdir temp/minimap2
for file in temp/cutadapt/*.fq
do
    for fastaFile in $inputDirectoryFastaFiles/*.fasta
    do
        filename_extended=${file##*/}
        searchstring="_full_length_output_cutadapt_2.fq"
        filename=${filename_extended%$searchstring*}
        fastaname_extended=${fastaFile##*/}
        fastaname=${fastaname_extended%.*}
        mkdir temp/minimap2/${fastaname}
        minimap2 -ax map-ont -k14 -t 8 ${fastaFile} ${file} > temp/minimap2/${fastaname}/"${filename}-${fastaname}_alignment.sam"
    done
done


for fastaFile in $inputDirectoryFastaFiles/*.fasta
do
    samtools faidx ${fastaFile}
done


mkdir ${output_alignment}
mkdir temp/SAM_and_BAM
for file in temp/minimap2/*/*.sam
do
    filename_extended=${file##*/}
    searchstring="_alignment.sam"
    filename=${filename_extended%$searchstring*}
    tempString=$(echo "$file" |awk -F'temp/minimap2/|/*.sam' '{print $2}')
    fastaName=$(echo $tempString | cut -f 1 -d"/")
    
    samclip --max 10 --ref  $inputDirectoryFastaFiles/${fastaName}.fasta < ${file} > temp/SAM_and_BAM/${filename}_clipped.sam
    samtools view -q 30 -S -b temp/SAM_and_BAM/${filename}_clipped.sam > temp/SAM_and_BAM/${filename}_clipped.bam 
    samtools sort temp/SAM_and_BAM/${filename}_clipped.bam -o  ${output_alignment}/${filename}.sorted.bam
    samtools index ${output_alignment}/${filename}.sorted.bam
done




