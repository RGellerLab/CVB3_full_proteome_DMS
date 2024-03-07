#!/bin/bash
#SBATCH --job-name=TrimFqBam
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 4 
#SBATCH --mem 12G
#SBATCH --time=0-06:00:00 
#SBATCH --output=archivo_%j.out 
#SBATCH --error=archivo_%j.err 

########### #######################################################
## trims the files using fastp to 150bp and does a bit of cleaning
## run from a directory where all the raw files are

### #######################before running this ###################################################
### split execution into chunks 4 batches to speed up
# ls -d */ |split -l 4

## now run as
# sbatch ~/final_scripts/s1_trim_makebam.bash xaa ## and then xab etc

## load modules

module load samtools/1.16
module load fastp/0.23 


while IFS= read -r d; do

cd $d
mkdir -p ./cut_files
mkdir -p ./bam_files

### make filename directory
find *_1.fq.gz >fnames.txt

## defines length to trim files to (and to remove smaller files)
MLength=150

### loop on each entry in file 
while IFS= read -r line; do
# trims the name at the _.fq.gz so that each have a unique ID with the lane number
	condition="$(echo $line | sed 's/_[0-9]\.fq\.gz//')"
	echo $condition
	 #name of file2 
	fileR2="$(echo $line |sed -e s/1.fq.gz/2.fq.gz/)" # replace R1 with R2
    o1=./cut_files/${condition}_cut_trimR1.fq.gz
    o2=./cut_files/${condition}_cut_trimR2.fq.gz

    # has -x for polyX tract, -Q for no qual, -A for no adaptor
    fastp -i $line -I $fileR2 -o $o1 -O $o2 \
	  --max_len1 $MLength --max_len2 $MLength \
	  --length_required $MLength \
	  -x -Q -A -w 8 
	  #-h ${workdir}/${condition}.html \

done < fnames.txt
rm fnames.txt

## enter cut files directory -2
cd ./cut_files/

# get all the R1 to know mates
find *R1.fq.gz >cutnames.txt

### loop on each entry in file 
while IFS= read -r line; do
sample="$(echo $line | sed 's/_cut_trimR1\.fq\.gz//')" # cuts by delim and join id
fileR2="$(echo $line |sed -e s/1.fq.gz/2.fq.gz/)" # replace R1 with R2
out="../bam_files/"${sample}".bam"
echo $sample
echo $fileR2
echo $out
java -jar ~/apps/picard.jar FastqToSam  F1=$line F2=$fileR2 O=$out SM=${sample}
done < cutnames.txt

## exit back to condition directory
cd ..
### generate merged bam
mkdir -p ./merged_bam
### -2
cd bam_files
ls *.bam > list_bam.txt
samtools cat -b list_bam.txt -o ../merged_bam/merged.bam
rm *.bam ## 
### go back out to main directory 
cd ..
### go back to start 
cd ..


done < $1