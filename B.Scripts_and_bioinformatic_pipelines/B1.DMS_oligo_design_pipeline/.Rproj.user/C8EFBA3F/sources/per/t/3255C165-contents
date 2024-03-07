
### finds all matches for each number of changes and gets table
library(Biostrings)
ref=readDNAStringSet("./CVB3_nancy.fa")
res=data.frame(matrix(data=NA,nrow=8,ncol=5))
colnames(res)=c("tile","zero","one","two","three")

for (i in 9:18){
  # i=1
  rg=readDNAStringSet(paste0("./Oligos/Tile_", i ,"_oligos.fa"))
  # temp.res=data.frame(matrix(data=NA,nrow=length(rg),ncol=4))
  three=sapply(rg, function(x) vmatchPattern(x,ref,min.mismatch = 3,max.mismatch = 3))      
  two=sapply(rg, function(x) vmatchPattern(x,ref,min.mismatch = 2,max.mismatch = 2))      
  one=sapply(rg, function(x) vmatchPattern(x,ref,min.mismatch = 1,max.mismatch = 1))      
  zero=sapply(rg, function(x) vmatchPattern(x,ref,max.mismatch = 0))      
  
  res[i-8,]=c(i,length(unlist(lapply(zero,function(x) startIndex(x)))),
            length(unlist(lapply(one,function(x) startIndex(x)))),
            length(unlist(lapply(two,function(x) startIndex(x)))),
            length(unlist(lapply(three,function(x) startIndex(x)))))
  }

write.csv(res,"./mutation_numbers.csv",row.names = F)


########
## trims overlaps for hifi and then align/translate
library(Biostrings)
library(tidyverse)
library(readxl)
library(DECIPHER) 
###Read ORF
seq<-readDNAStringSet("./CVB3_nancy.fa")

output.dir="./Alignments/"
dir.create(output.dir,showWarnings = F)

###Read tile and overlaps coordinates
coord<-data.frame(read_xlsx("./Tiles_overlaps_coord.xlsx"))
tiles.seq<-data.frame("index"=character(), "overlap_5"=character(), "ref_seq"=character(), "overlap_3"=character(), "start_aa"=character())
tiles.seq=bind_cols(X5_start=coord$X5_start,
                    X5_end=coord$X5_end,
                    X3_start=coord$X3_start,
                    X3_end=coord$X3_end)

for (i in 1:10){
  # i=1
  rg=readDNAStringSet(paste0("./Oligos/Tile_", i+8,"_oligos.fa"))
  x5=as.character(subseq(seq,start = tiles.seq$X5_start[i],tiles.seq$X5_end[i]))
  x3=as.character(subseq(seq,start = tiles.seq$X3_start[i],tiles.seq$X3_end[i]))
  rg=trimLRPatterns(Lpattern = x5,Rpattern = x3,rg)
  ref2=subseq(seq,tiles.seq$X5_end[i]+1,tiles.seq$X3_start[i]-1)
  rg=c(ref2,rg)
  aarg=translate(rg)
  rg=AlignTranslation(rg)
  
  writeXStringSet(rg,paste0(output.dir, "tile_",i+8,"_trim_align.fa"))
  writeXStringSet( aarg,paste0(output.dir, "tile_",i+8,"_AA_trim_align.fa"))
  # rm(x5,x3,,rg,i,aarg)
}


### Count oligos
library(Biostrings)
num.oligos=data.frame(matrix(data=NA,nrow=10,ncol=2))
colnames(num.oligos)=c("tile","number_oligos")

for (i in 9:18){
        # i=1
        rg=readDNAStringSet(paste0("./Oligos/Tile_", i ,"_oligos.fa"))
        num.oligos[i-8,]=c(i,length(rg))
}

num.oligos[nrow(num.oligos)+1,]<-c("total", sum(num.oligos$number_oligos))
write.csv(num.oligos,"./number_tiles.csv",row.names = F)
