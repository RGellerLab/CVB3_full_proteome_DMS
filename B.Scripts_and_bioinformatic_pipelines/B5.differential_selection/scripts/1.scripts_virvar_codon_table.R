## all dms scripts

#utilites
# strand odds ratio calculation based on GATK info
# june 22 2023, RG
SOR=function(wtF, wtR, mutF, mutR){
  # Add pseudocount
  wtF=wtF+1;  wtR=wtR+1  ; mutF=mutF+1;  mutR=mutR+1
  r= (wtF*mutR)/(wtR*mutF)
  R= r+1/r ## symetrical ratio
  refRatio=(min(wtF,wtR)/ max(wtF,wtR))
  mutRatio=(min(mutF,mutR)/ max(mutF,mutR))
  
  res=log(R)+log(refRatio)- log(mutRatio)
  return(res)
}

# june 22 2023, RG
virvarseq_to_codon_table=function(file, name, f.length, 
                                  strand.test=c("fisher","SOR","none"),
                                  SOR.cut=4,
                                  out.directory, ref){
  # takes vir var seq output file, filters it, removes results that fail strand bias, and writes both
  # cleaned file and a codon table of the results)
  
  library(Biostrings)
  
  ### read wt sequence and get codons
  refseq=as.character(readDNAMultipleAlignment(ref))
  refseq=strsplit(refseq,"")[[1]]
  ref.codons=sapply(seq(from = 1,to = (length(refseq)), by = 3), function(x) paste0(refseq[x:(x+2)],collapse = ""))
  rm(ref,refseq)
  
  ## get genetic code
  gen_code=as.data.frame(cbind(codons=as.character(names(GENETIC_CODE)),AA=GENETIC_CODE),stringsAsFactors=F)
  gen_code=gen_code[order(gen_code$codons),]
  
  ##make results df
  res.df=data.frame(matrix(data=0, nrow=f.length,ncol=66))
  colnames(res.df)=c("site","wildtype",gen_code$codons)
  res.df$site=1:f.length
  rm(f.length)
  
  #### read vir var seq data
  vvs=read.delim(file,stringsAsFactors = F) 
  vvs$REF_CODON=toupper(vvs$REF_CODON)
  
  # remove amino acids that are X
  vvs=vvs[vvs$CODON%in%gen_code$codons,]
  
  ######## remove failure from strand bias, need to recalculate denom after removing X
  positions=unique(vvs$POSITION)
  
  print(paste0("working on file ",name))
  
  #### redo denom
  for (i in positions){ 
    #print(i)
    # i=831############
    ## fix up denominator
    vvs[vvs$POSITION==i,"FWD_DENOM"]=sum(vvs[vvs$POSITION==i,"FWD_CNT"])
    vvs[vvs$POSITION==i,"REV_DENOM"]=sum(vvs[vvs$POSITION==i,"REV_CNT"])
    vvs[vvs$POSITION==i,"DENOM"]=sum(vvs[vvs$POSITION==i,"CNT"])
    rm(i)
    
  }
  rm(positions)  
  
  
  #### do strand test for strand bias
  positions=unique(vvs$POSITION)
  vvs$std.bias=1
  vvs$id=1:nrow(vvs) #give an ID for each row 
  
  #print(paste0("fixing strand bias ",name))
  for (i in positions){
    #print(i)
    vvs2=vvs[vvs$POSITION==i,] ## vvs of individual AA
    
    if(nrow(vvs2)>1){  ### if there are mutants to test get WT and mutant rows
      wt=vvs2[vvs2$CNT==max(vvs2$CNT),] ## WT is defined as most common
      muts=vvs2[vvs2$CNT!=max(vvs2$CNT),]  
      
      if (strand.test=="fisher"){
        ### iterate overall mutants to test vs wt
        for (j in 1:nrow(muts)){
          # j=1
          vvs$std.bias[vvs$id==muts$id[j]]=fisher.test(rbind(c(wt$FWD_CNT, wt$REV_CNT),  #### pvalue calculation for bias
                                                             c(muts$FWD_CNT[j], muts$REV_CNT[j])))$p.value
          rm(j)
        }
        rm(i,vvs2,wt,muts)
      }
      
      if (strand.test=="SOR"){
        ### iterate overall mutants to test vs wt
        for (j in 1:nrow(muts)){
          # j=1
          vvs$std.bias[vvs$id==muts$id[j]]=SOR(wt$FWD_CNT, wt$REV_CNT,muts$FWD_CNT[j], muts$REV_CNT[j])
          rm(j)
        }
        rm(i,vvs2,wt,muts)
      }
      
      if (strand.test=="none"){print("not fixing strand")
      }
      
    }
  }
  
  rm(positions)
  
  ## correct frequency
  vvs$FREQ=vvs$CNT/vvs$DENOM
  vvs=vvs[,-ncol(vvs)] ## remove ID column
  
  ## strand.test 
  
  if (strand.test==T){
    p.val_strand_bias=0.05/nrow(vvs)
    rejected=vvs[vvs$std.bias<=p.val_strand_bias,] # those that would get kicked out normally
    vvs=vvs[vvs$std.bias>p.val_strand_bias,] #strand bias if chosen
    rm(p.val_strand_bias)
  } 
  
  ## fdr correction
  if (strand.test=="SOR"){
    rejected=vvs[vvs$std.bias>=SOR.cut,] # those that would get kicked out normally
    vvs=vvs[vvs$std.bias<SOR.cut,] #strand bias if chosen
  } 
  
  
  
  # iterate over rows of data and put in for each position/codon put in count data
  for (r in 1:nrow(vvs)){res.df[vvs$POSITION[r],vvs$CODON[r]]=vvs$CNT[r]; rm(r)}
  
  ## add wildtype to codon file
  res.df$wildtype=ref.codons
  
  
  #### write results     ##########
  ### create directories
  dir.create(  out.directory, recursive = T)
  dir.create(  paste0(out.directory,"/codon_tables/"), recursive = T)
  dir.create(  paste0(out.directory,"/rejected_codon_tables/"), recursive = T)
  dir.create(  paste0(out.directory,"/filtered_tables/"), recursive = T)
  ### write files
  write.csv(res.df,paste0(out.directory,"/codon_tables/",name, "_codon_table.csv"),row.names = F)
  write.csv(rejected,paste0(out.directory,"/rejected_codon_tables/",name, "_rej_codon_table.csv"),row.names = F)
  write.csv(vvs,paste0(out.directory,"/filtered_tables/",name, "_filtered.csv"),row.names = F)
  write(c(SOR.cut,strand.test),paste0(out.directory,"/codon_tables/", "strandtest.txt"))
  
  # test if wt is not the max codon, and if so, write in file
  fix.mut=sum(ref.codons!=apply(res.df[3:ncol(res.df)], 1 , function (x) colnames(res.df[3:ncol(res.df)])[which(x== max(x))]))
  if (fix.mut>1){
    write(name,file = paste0(out.directory, "fixed_mutations.txt"),append = F)}
  
  
}


# june 22 2023, RG
step1.collapse_codon_files=function(codon.dir="../codontables/",
                                    out.dir="../results/unfiltered/collapsed_codon_files/",
                                    out.dir.dels="../results/unfiltered/collapsed_deletion_files/",
                                    prefix="sscs.codon")
{
  
  ## this function collapses all the codon files that belong to the same sample
  ## and writes them out as a single file for each sample
  ## into an unfiltered or a deletion file folder
  
  ## codon.dir="../codontables/" is where codon files are from virvarseq
  ##out.dir= is where mutation codon file go 
  ## out.dir.dels, where dels go
  ## the d
  library(tidyverse)
  ## just run for each folder to collapse... can also loop
  ### data should be in codon tables
  #wd="../codontables/"
  ## output will be in unfiltered collapsed
  #out.dir="../results/unfiltered/collapsed_codon_files/"
  dir.create(out.dir,recursive = T)
  # out.dir.dels="../results/unfiltered/collapsed_deletion_files/"
  dir.create(out.dir.dels,recursive = T)
  fs=list.files(codon.dir,full.names = T,include.dirs = T)
  
  ## get conditions
  conds <- unique(sapply(basename(fs), function(x) {
    if (grepl("-", x)) {
      unlist(strsplit(x, split = "-"))[1]
    } else {
      unlist(strsplit(x, split = "\\."))[1]
    }
  }))
  
  
  ### make a df to make sure name swap is ok
  df.files=data.frame(matrix(data=NA,nrow=0,ncol=3))
  colnames(df.files)=c("cond","fname","fname.cond")
  for (c in conds){
    # c=conds[2]
    
    ## finds condition name in list of files
    fs.cond=fs[grep(paste0(c),fs)]
    
    # gets name from files
    nm.cond=c
    
    # reads the first of these files
    df=read_delim(fs.cond[1],delim="\t")
    
    if("FWD_MEAN_MIN_QUAL"%in%names(df)){
    df=df %>% filter(AA!="X") %>% 
      select(-c(FWD_MEAN_MIN_QUAL:REV_STDDEV_MIN_QUAL,FREQ, SAMPLE))
    }else{df=df %>% filter(AA!="X")}
    df$POSITION=as.numeric(df$POSITION)
    
    ## add names to df of files to make sure ok
    df.files=rbind(df.files,data.frame(cond=c,fname=fs.cond[1],fname.cond=nm.cond))
    
    ## if there are more files b/c sam was split, then combine
    if (length(fs.cond)>1){
      for (f in 2:length(fs.cond)){
        # f=2
        ## add names to df of files to make sure ok
        df.files=rbind(df.files,data.frame(cond=c,fname=fs.cond[f],fname.cond=nm.cond))
        
        # read additional files and process
        df2=read_delim(fs.cond[f],delim="\t")
        df2$POSITION=as.numeric(df2$POSITION)
        df2=df2 %>% filter(AA!="X") %>% 
          select(-c(FWD_MEAN_MIN_QUAL:REV_STDDEV_MIN_QUAL,FREQ,SAMPLE))
        df=bind_rows(df,df2)
        ## some checks works well
        # df5=aggregate(. ~ POSITION+REF_CODON+CODON+REF_AA+AA,data=df3,FUN=sum) %>%  arrange(POSITION,CODON)
        # dupls=df3 %>% group_by(POSITION,REF_CODON,CODON,REF_AA,AA) %>% filter(n()>1) %>% arrange(POSITION,CODON)
        # df4[df4$POSITION==335 &df4$CODON=="ACA",]
        rm(df2)
      } 
    }
    df=df %>% group_by(POSITION,REF_CODON,CODON,REF_AA,AA) %>% 
      summarise_all(sum)  %>% arrange(POSITION,CODON)
    #separate dels into separate DF and write into sep directory
    dels=filter(df,AA=="-")
    write_delim(dels,paste0(out.dir.dels,nm.cond,prefix),delim = "\t")
    # remove dels and write into codon dir
    df=filter(df,AA!="-")
    write_delim(df,paste0(out.dir,nm.cond,prefix),delim = "\t")
    rm(df,c,f,fs.cond,nm.cond)
  }
}

# june 22 2023, RG
step2.get_codon_tables=function(data.dir= "../results/unfiltered/collapsed_codon_files/",
                                aa.fragment.length=628,
                                out.directory.results="../results/unfiltered/",
                                strand.test.method="SOR",
                                SOR.cut.value=4,
                                fasta.reference="./1.collapse_vivarseq_get_codon_tables/p1dms.fa",
                                rename=T,
                                name.file="./p2dms_library_info.xlsx"## needs name and ref.name
)
  
{ 
  
  
  ## get codon tables from virvarseq output from step 1
  #data.dir= "../results/unfiltered/collapsed_codon_files/"# input where the codon files are found
  fs=list.files(path = data.dir,pattern = "*.codon$",include.dirs = F)
  
  ## run script on all files
  ## produces folder with codon_tables (after filtering by SOR if selected)
  ## and folders with the filtered tables and rejected mutations
  for (f in fs){ 
    #f=fs[1]
    out.name=unlist(strsplit (f,split = "d3"))[1]
    if (rename==T){
      library(readxl)
      names.df=read_xlsx(name.file,sheet = 1)
      nm=gsub("sscs.codon","",basename(f))
      out.name=names.df$ref.name[names.df$name==nm]
    }
    if (rename==F){
      out.name=gsub("sscs.codon","",basename(f))
    }
    virvarseq_to_codon_table(file = paste0(data.dir,f) ,
                             name = out.name, 
                             out.directory=out.directory.results,
                             strand.test=strand.test.method, 
                             SOR.cut=SOR.cut.value,
                             ref=fasta.reference,
                             f.length=aa.fragment.length);
    rm (f)}
  
  
  rm(fs)
}


# june 22 2023, RG
step2b.filter_single_codons=function(min.cov=0,
                                     codon.identity.matrix="./codon_identity_matrix.csv", 
                                     data.dir="../results/unfiltered/codon_tables/",
                                     out.dir="../results/filter_single_mutants/codon_files/"){
  
  ## filter codon table for more than 1 mutation, uses output 
  ## from virvarseq_to_codon_table function
  
  # PARAMETERS##
  # min.cov=0 # minimum coverage for accepting mutation, if below put 0
  # codon.identity.matrix="./codon_identity_matrix.csv"
  # data.dir="./" # where the data for the codon files is located
  # out.dir="./filter_single_mutants/") # where you want to write the files
  
  ## load data on which codons are a single mutation away
  cod.tbl=read.csv(codon.identity.matrix,stringsAsFactors = F,row.names = 1)
  
  ## create output director
  dir.create(out.dir,recursive = T,showWarnings = F)
  ## iterate on all codon files in data directory
  fs=list.files(data.dir,pattern = "*.csv",full.names = T)
  for (f in fs){
    df=read.csv(f,stringsAsFactors = F)
    ## iterate on row, get columns that are single mutants and change to 0
    for( x in 1:nrow(df)){
      df[x, colnames(cod.tbl)[which(cod.tbl[df$wildtype[x],]==1)]]=0
      rm(x)
    }
    ### put at 0 all below threshold
    df2=df[,3:ncol(df)]
    df2[df2<min.cov]=0
    df[,3:ncol(df)]=df2
    rm(df2)
    
    # write results, name obtained by taking first element split by _ from
    # codon file name
    # write.csv(x = df, file = paste0(out.dir,unlist(strsplit(basename(f),"_"))[1],
    #                                 "_2_3mut_codon_table.csv"),row.names = F)
    write.csv(x = df, file = paste0(out.dir,basename(f)),row.names = F)
    rm(f,df)
  }
}
# june 22 2023, RG
step2b.filter_expected=function(
    data.dir="../results/unfiltered/codon_tables/",
    out.dir="../results/filter_expected/codon_files/",
    out.dir2="../results/filter_not_expected/codon_files/",
    muts.introduced="./1.collapse_vivarseq_get_codon_tables/mutations_introduced_in_P3.csv"){
  
  #filter deletions
  muts=read_csv(muts.introduced)%>% filter(!is.na(new.aa))
  dir.create(out.dir,recursive = T,showWarnings = F)
  dir.create(out.dir2,recursive = T,showWarnings = F)
  
  ## list of files
  fs=list.files(data.dir,full.names = T,pattern = "*_codon_table.csv")
  
  #### loops on files and filters
  for (f in 1:length(fs)){
    # f=1
    df=read_csv(fs[f])
    nm=gsub("_codon_table.csv","",basename(fs[f]))
    real.nm=nm#names.df$ref.name[names.df$name==nm]
    expec=df
    n.expec=df
    # names.df$check_name[names.df$name==nm]=real.nm
    #names.df$check_file[names.df$name==nm]=fs[f]
    ## go AA by AA and filter out all matching the changes or not matching
    for (c in unique(df$wildtype)){
      # c="ATG"
      temp.muts=filter(muts,codon==c)
      expec[expec$wildtype==c,!colnames(expec)%in%c(temp.muts$new.codon,"wildtype","site",c)]=0
      ### n.expect also removes WT
      n.expec[n.expec$wildtype==c,colnames(n.expec)%in%temp.muts$new.codon]=0
      rm(c,temp.muts)
    }
    write_csv(expec,paste0(out.dir,nm,"_codon_table.csv"))
    write_csv(n.expec,paste0(out.dir2,nm,"_codon_table.csv"))
    rm(df,nm,expec,n.expec)
  }
}

# june 22 2023, RG
step2c.filter_syn=function(data.dir="../results/filter_expected/codon_files/",
                           out.dir="../results/filter_expected/",
                           codon.identity="./1.collapse_vivarseq_get_codon_tables/codon_identity_matrix.csv"){
  
  library(tidyverse)
  library(Biostrings)
  
  dir.create(paste0(out.dir,"/syn/"),showWarnings = F)
  ### read data
  fnames=list.files(data.dir,full.names = T,
                    include.dirs = F)
  gcode=data.frame(GENETIC_CODE) 
  gcode=gcode%>% mutate(codon=toupper(rownames(gcode)),aa=toupper(GENETIC_CODE)) %>% select(-1)
  # identity matrix
  id=read_csv(codon.identity)
  colnames(id)[1]="codon"
  for (f in fnames){
    # f=fnames[1]
    df=read_csv(f)
    ns=single=double=triple=df
    s=df
    for (r in 1:nrow(df)){
      # r=1
      ## split into syn and ns
      wt=df$wildtype[r]
      aa=gcode$aa[gcode$codon==wt]
      #ns codons
      ns.cods=gcode$codon[gcode$codon!=wt &gcode$aa!=aa]
      #s codons
      ss.cods=gcode$codon[gcode$codon!=wt &gcode$aa==aa]
      ## if there are syn then make them 0 in NS
      if (length(ss.cods)>0){ns[r,ss.cods]=0}
      # if there are non-syn, then make the 0 in syn
      if (length(ns.cods)>0){s[r,ns.cods]=0}
      rm(ns.cods,ss.cods)
      ### now split ns by numer of mutants
      ids=id[id$codon==wt,] %>% select(-1)
      ones=colnames(ids)[ids!=1]
      twos=colnames(ids)[ids!=2]
      threes=colnames(ids)[ids!=3]
      if (length(ones)>0){single[r,ones]=0}
      if (length(twos)>0){double[r,twos]=0}
      if (length(threes)>0){triple[r,threes]=0}
      rm(ones,twos,threes,ids)
      rm(r)
      rm(wt,aa)
    }
    nm=basename(f)  
    write_csv(s,paste0(out.dir,"/syn/",nm))
    rm(df)
  }
}

# june 22 2023, RG
step2d.filter_123muts=function(data.dir="../results/filter_expected/codon_files/",
                               out.dir="../results/filter_expected/",
                               codon.identity="./1.collapse_vivarseq_get_codon_tables/codon_identity_matrix.csv"){
  
  library(tidyverse)
  library(Biostrings)
  
  dir.create(paste0(out.dir,"/single/"),showWarnings = F)
  dir.create(paste0(out.dir,"/double/"),showWarnings = F)
  dir.create(paste0(out.dir,"/triple/"),showWarnings = F)
  ### read data
  fnames=list.files(data.dir,full.names = T,
                    include.dirs = F)
  gcode=data.frame(GENETIC_CODE) 
  gcode=gcode%>% mutate(codon=toupper(rownames(gcode)),aa=toupper(GENETIC_CODE)) %>% select(-1)
  # identity matrix
  id=read_csv(codon.identity)
  colnames(id)[1]="codon"
  for (f in fnames){
    # f=fnames[1]
    df=read_csv(f)
    ns=single=double=triple=df
    s=df
    for (r in 1:nrow(df)){
      # r=1
      ## split into syn and ns
      wt=df$wildtype[r]
      aa=gcode$aa[gcode$codon==wt]
      #ns codons
      ns.cods=gcode$codon[gcode$codon!=wt &gcode$aa!=aa]
      #s codons
      ss.cods=gcode$codon[gcode$codon!=wt &gcode$aa==aa]
      ## if there are syn then make them 0 in NS
      if (length(ss.cods)>0){ns[r,ss.cods]=0}
      # if there are non-syn, then make the 0 in syn
      if (length(ns.cods)>0){s[r,ns.cods]=0}
      rm(ns.cods,ss.cods)
      ### now split ns by numer of mutants
      ids=id[id$codon==wt,] %>% select(-1)
      ones=colnames(ids)[ids!=1]
      twos=colnames(ids)[ids!=2]
      threes=colnames(ids)[ids!=3]
      if (length(ones)>0){single[r,ones]=0}
      if (length(twos)>0){double[r,twos]=0}
      if (length(threes)>0){triple[r,threes]=0}
      rm(ones,twos,threes,ids)
      rm(r)
      rm(wt,aa)
    }
    nm=basename(f)  
    write_csv(single,paste0(out.dir,"/single/",nm))
    write_csv(double,paste0(out.dir,"/double/",nm))
    write_csv(triple,paste0(out.dir,"/triple/",nm))
    rm(df)
  }
}
# june 22 2023, RG
step3.filter_dels=function(min.cov=1,
                           max.sor=4,
                           wd="../results/unfiltered/collapsed_deletion_files/",
                           wt.wd="../results/unfiltered/codon_tables/",
                           out.dir="../results/unfiltered//deletions/",
                           rename=T,
                           name.file="../p2dms_library_info.xlsx"){
  library(tidyverse)
  
  ### filters dels out of codon files
  ### adds SOR 
  ### filter by max SOR and min cov
  min.cov=min.cov
  max.sor=max.sor
  
  #
  wd=wd#"../results/unfiltered/collapsed_deletion_files/"
  wt.wd=wt.wd#"../results/unfiltered/codon_tables/"
  out.dir=out.dir#"../results/unfiltered//deletions/"
  
  dir.create(out.dir,recursive = T,showWarnings = F)
  fs=list.files(wd,full.names = T,include.dirs = T)
  fs.wt=list.files(wt.wd,full.names = T,include.dirs = T)
  
  
  if (rename==T){
    library(readxl)
    names.df=read_xlsx(name.file,sheet = 1)
  }
  
  #### iterate on files
  for (f in 1:length(fs)){
    # f=1
    ## read the deletion data
    df=read_tsv(fs[f])
    
    # get the name
    nm=unlist(strsplit(basename(fs[f]),split = "\\."))[1]
    nm=gsub("sscs","",nm)
    if (rename==T){
      nm=names.df$ref.name[names.df$name==nm]}
    
    ## get the relevant WT and its coverage
    df.wt=read_csv(fs.wt[grep(nm,fs.wt)]) %>% 
      mutate(coverage=rowSums(select(.,-c(1:2)))) %>% select(1,2,coverage)
    # wt.cnt=rowSums(select(df.wt,-c(1:2)))#sapply(1:nrow(df.wt),function(x) df.wt[x,df.wt$wildtype[x]])
    df$sor=NA
    ## filter only dels and min of 1 or of
    dels=df %>% filter(AA=="-") %>% filter(CNT>= min.cov)
    for (i in 1:nrow(dels)){
      # i=1
      dels$sor[i]=
        SOR(wtF = dels$FWD_DENOM[i], wtR = dels$REV_DENOM[i],
            mutF = dels$FWD_CNT[i], mutR = dels$REV_CNT[i])
      
      # if (dels$REV_DENOM/dels$FWD_DENOM>100 |dels$FWD_DENOM/dels$REV_DENOM>100)
      
      rm(i)
    }
    
    dels=filter(dels, sor<max.sor) %>% mutate(REF_CODON=toupper(REF_CODON))
    res=left_join(df.wt,dels,by=c("site"="POSITION", "wildtype"="REF_CODON"))
    
    write_csv(res,paste0(out.dir,nm,".csv"))
    rm(df,nm,f)
  }
}