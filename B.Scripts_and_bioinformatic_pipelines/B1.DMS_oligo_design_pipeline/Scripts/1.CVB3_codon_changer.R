###Generates a table that associates to each codon a new codon with:
### 1) max number of different bp
### then
### 2) max frequency in CVB3 ORF
### and
### empty codon (deletion)

###Requires:
###CVB3_codon_usage.xlsx 
###codon_identity_matrix.csv

###Output:CVB3_codon_changer_table.csv

library(Biostrings)
library(tidyverse)
library(readxl)

# codons=names(GENETIC_CODE)
aas=names(AMINO_ACID_CODE)
aas=aas[!aas %in% c("U","O","B","J","Z","X")]

###Define codon changer empty data frame
codon_changer=data.frame("codon"=character(),"new.codon"=character(),"aa"=character(),"new.aa"=character())

###Data frame to translate
codon.to.aa=as.data.frame(GENETIC_CODE)
codon.to.aa$codon=rownames(codon.to.aa)
colnames(codon.to.aa)<-c( "aa","codon")



###Read matrix of codon identity
codon.identity<-data.frame(read.csv("./codon_identity_matrix.csv"), row.names = TRUE)
codon.identity$ref_codon=rownames(codon.identity)
###Read CVB3 codon usage
codon.usage<-data.frame(read_xlsx("./CVB3_codon_usage.xlsx"))

codon.usage=left_join(codon.to.aa, codon.usage,by=c("codon"="Codon"))
rm(codon.to.aa)
for(c in names(GENETIC_CODE)){
        # c="TCA"
        #j=1
        
        ###translate codon
        aa=codon.usage$aa[codon.usage$codon==c]
        ###For each aa
        for(k in aas){
                #print(k)
                #k="S"
                #c="TCA"
                
                ###Find all syn codons for each AA to change
                possible.cod=filter(codon.usage,aa==k) 
                
                ## merge with # of changes relative to original codon
                ## sort by number of changes and by freq and select top
                changes=select(codon.identity,c,ref_codon) %>% 
                        filter(ref_codon%in%possible.cod$codon)
                changes=left_join(possible.cod,changes,by=c("codon"="ref_codon")) %>% 
                        rename(select.codon=4) %>% 
                        arrange(desc(select.codon),desc(Fraction))
                
                new.cod=changes$codon[1]
                
                
                ###Add to table: ref codon, new codon, ref aa, new aa
                codon_changer[nrow(codon_changer)+1,]<-c(c,new.cod,aa,k)
                
        }
        ###Add empty codon to table 
        
        codon_changer[nrow(codon_changer)+1,]<-c(c,"",aa,"")
        
}
write.csv(codon_changer,"CVB3_codon_changer_table.csv",row.names = F)

