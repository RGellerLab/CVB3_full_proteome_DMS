#!/bin/bash
#SBATCH --job-name=sscs
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 2 # “-c”
#SBATCH --mem 24G
#SBATCH --time=3-00:00:00 #### 0-00:00:00
#SBATCH --output=archivo_%j.out
#SBATCH --error=archivo_%j.err


### modules to load

module load bwa/0.7.17   
module load samtools/1.16
module load biotools ### for virvarseq
module load variantbam/0.0.0

## perl modules for virvarseq
export PERL5LIB=/home/ronge/apps/VirVarSeq/lib
export R_LIBS_USER=/home/ronge/apps/VirVarSeq/R/lib

# now run as
# sbatch ~/final_scripts/i2sysbio_pipeline/s2_sscs.bash xaa ## and then xab etc

### input correct reference and wd as current wd
### reference must be indexed with bwa and picard

#P1
#ref="/home/ronge/refs/p1_dms/hifi.fa"
#mutstart=76 #
#lenAA=851

#p2
#ref="/home/ronge/refs/p2_dms/p2dms.fa"
#mutstart=1
#lenAA=628


## p3
#ref="/home/ronge/refs/p2_dms/p3dms.fa" # positions are 98 to 2473
#mutstart=98
#lenAA=770

# where is working dir
#get with pwd
conddir=$(pwd)

## this is where the final codon files will go to...
outputfs="/codontables/"
statsfs="/stats/"
codondir="$(echo $conddir | sed 's/raw_data//')"
finaloutput="$codondir$outputfs"
finalstats="$codondir$statsfs"
# where codon files will go
mkdir -p $finaloutput
# where stat files will go for coverage
mkdir -p $finalstats

## global variables for duplex
### DEPENDING ON THE ADAPTERS, WE NEED TO DEFINE DIFFERENT NUMBER OF Ns
### THE DUAL INDEX ARE 8 ON EACH SIDE
############ make sure you select length of index correctly
indlength=8 # index length 8 or 12 depend on version of adaptors used
tagrm="$(($indlength+5))"  ## length of index and invariable region
cutoff=0.7 # min cutoff fraction to be consensus
ReadPerDup=2 # min for consesus
trim=7 #trimming on termini to cleanup
ncut=0.5 # faction of N allowed


### filname to work with for looping, this is the xaa etc. 
fname=$1


while IFS= read -r line; do

# enter the directory and make the bam directory and merged bam directory
cd $line
##remove "/" from sample
sample="$(echo $line |tr -d "/")" 

#### make and enter duplex directory
mkdir duplex -p
cd ./duplex # in duplex directory

### generate sample name
name=${sample}_ind${indlength}_fam${ReadPerDup}_min${cutoff}



## run duplex script #######################################
echo "running duplex"
## the slowest step
python3.6 ~/apps/Duplex-Seq-Pipeline/scripts/UnifiedConsensusMaker.py --input ../merged_bam/merged.bam  --prefix ${sample} --tagstats --cutoff $cutoff --minmem $ReadPerDup --write-sscs --Ncutoff $ncut

# copy stats over to folder
cp *tagstats.txt $finalstats
cp *cmStats.txt $finalstats

## clean up some extra files we don't need
rm *_dcs.fq.gz
rm *aln_seq1.fq.gz
rm *aln_seq2.fq.gz

################### align  and rm unaligned ##################################
echo "align and sort"
bwa mem $ref *_read1_sscs.fq.gz *_read2_sscs.fq.gz | samtools sort - -o sort.bam 


echo "remove unaligned"
samtools view sort.bam -b -F 4 -o filt.bam

## rm intermediate bams
rm sort.bam #cut.bam

##index
samtools index filt.bam #-@ $cpus


## to look at error by pos if want. looks like trim 7 is enough
#java -jar ~/apps/fgbio-2.1.0.jar ErrorRateByReadPosition -i filt.bam 

samtools sort -n -u filt.bam -o sorted.bam

## upgrades soft clip (those don't match or have low quality)
## to hard clip, which actually removes them from file (seq AND QUAL)
## then cuts overlaps between mates so 50% in each
java -jar ~/apps/fgbio-2.1.0.jar ClipBam \
-i sorted.bam \
-o hclip_overlap.bam \
-c Hard \
--upgrade-clipping true \
-r $ref \
--clip-overlapping-reads true \
--read-one-five-prime $trim \
--read-one-three-prime $trim \
--read-two-five-prime $trim \
--read-two-three-prime $trim \
-m clip_info2.txt

rm sorted.bam
rm filt.bam


## Max 0 ins, 4 del, or 4 snps, and filter size and mate
variant hclip_overlap.bam -r '{"":{"rules":[{"ins":[0,0],"del":[0,4],"nm":[0,4],"mate_mapped":true,"fr":true,"length":[30,150]}]}}' -o indel_snp.bam -b


#####   codons with virvarseq
# have to install R library first as indicated in readme
#set up dirs
mkdir -p sscs_trim
mkdir -p sscs_trim/codon_table
mkdir -p sscs_trim/map_vs_consensus


# make sam file for program to work on if not very big############################
# samtools view -h indel_snp.bam -@ $cpus -m $mem > ./sscs_trim/map_vs_consensus/$name.sam
# echo $name > samples_sscs.txt


## option for splitting bams if very large to make faster##########################
cd ./sscs_trim/map_vs_consensus/
samtools view -H ../../indel_snp.bam > header
samtools view ../../indel_snp.bam | split - ${name}_ -l 8000000 --filter='cat header - > $FILE.sam'
rm header
ls | sed -e 's/\.sam$//'> ../../samples_sscs.txt
cd ../../ ## back to name directory of duplex
##

~/apps/VirVarSeq/codon_table.pl --samplelist ./samples_sscs.txt --ref $ref --outdir ./sscs_trim --start=$mutstart --len=$lenAA --trimming=0 --qual=20 >> ./VirVarSeq.log 2>&1

echo "done"
date

## remove extra files
rm ./sscs_trim/map_vs_consensus/*.sam
cp  sscs_trim/codon_table/*.codon $finaloutput

## clean up... 
rm *.temp.sort.bam
rm hclip_overlap.bam
rm indel_snp.bam



# exit to condition directory
cd $conddir



done < $fname
