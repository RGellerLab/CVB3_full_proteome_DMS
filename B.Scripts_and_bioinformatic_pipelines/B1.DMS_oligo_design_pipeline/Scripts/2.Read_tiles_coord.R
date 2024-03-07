###Takes coordinates of tiles and overlaps and generates a table with the 
###corresponding sequences

###Requires:
###CVB3_nancy.fasta
###Tiles_overlaps_coord.xlsx

###Output:Tiles_overlaps_seq.csv

library(Biostrings)
library(tidyverse)
library(readxl)

###Read ORF
seq<-readDNAStringSet("./CVB3_nancy.fa")


###Read tile and overlaps coordinates
coord<-data.frame(read_xlsx("./Tiles_overlaps_coord.xlsx"))
n.tiles=nrow(coord)

###Create empty data frame where sequences will be stored
tiles.seq<-data.frame("index"=character(), "overlap_5"=character(), "ref_seq"=character(), "overlap_3"=character(), "start_aa"=character())

X5_start=coord$X5_start
X5_end=coord$X5_end
ref_seq_start=coord$ref_seq_start
ref_seq_end=coord$ref_seq_end
X3_start=coord$X3_start
X3_end=coord$X3_end


for (n in 1:n.tiles){
        # n=1
        overlap_5<-subseq(seq, start=X5_start[n], end=X5_end[n])
        ref_seq<-subseq(seq, start=ref_seq_start[n], end=ref_seq_end[n])
        overlap_3<-subseq(seq, start=X3_start[n], end=X3_end[n])
        start_aa<-(ref_seq_start[n]+2)/3
        
        tiles.seq[nrow(tiles.seq)+1,]<-c(paste0("Tile_",n+8), overlap_5, ref_seq, overlap_3, start_aa)
}

write.csv(tiles.seq,"Tiles_overlaps_seq.csv",row.names = F)

