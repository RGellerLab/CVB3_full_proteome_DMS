# running stats

# june 22 2023, RG
analyze_codon_table=function(data.dir= "../results/filter_single_mutants/",
                             out.dir="../results/analyze_codon_tables/",
                             filter.positions=T,
                             start=18,
                             end=605,
                             correct.position=T,
                             aa.pos.correct=847,
                             min_reads=0){
  ## this function takes the codon table file and calculates some statistics, looks at coverage, 
  ## and produces some graphs to analyze the data
  
  # file =fnames[1]
  # name=samples[1]
  # start=18
  # end=605
  # aa.pos.correct=847
  # data.dir= "../results/filter_single_mutants/"
  # out.dir="../results/analyze/"
  # min_reads=0
  
  library(tidyverse)
  ## read files
  fnames=list.files(data.dir,full.names = T,#pattern = "sscs.codon",
                    include.dirs = F)
  ### full dir for data
  dir.create(out.dir,recursive = T)
  for (f in fnames){
    # f=fnames[1]
    
  ## read file
  df=read_csv(f)
  
  # get name
  name=gsub("_codon_table.csv","",basename(f))
  
  # filter any positions not in the dms part
  if (filter.positions==T){
  df=filter(df,site%in%c(start:end))
  }
  
  # correct positions
  if (correct.position==T){
  df$site=(df$site-start)+aa.pos.correct
  

  }
  
  library(Biostrings)
  
  ### get aa tables
  aa.df=df%>% select(-wildtype) %>% 
    pivot_longer(AAA:TTT)
  aa.df$aa=as.character(translate(DNAStringSet(aa.df$name),no.init.codon = T))
  
  res.aa=aa.df %>% select(-name) %>% 
    group_by(site,aa) %>% summarize(sum=sum(value)) %>% 
    pivot_wider(names_from = aa,values_from=sum )
  rm(aa.df)
  ## add wt aa
  res.aa$wildtype=as.character(translate(DNAStringSet(df$wildtype),no.init.codon = T))
  res.aa=res.aa %>% select(1,wildtype,"A":"Y","*")
  ### write
  dir.create(paste0(out.dir,"/aa_tables/"),recursive = T,showWarnings = F)
  
  write_csv(res.aa,paste0(out.dir,"/aa_tables/",name,"_aa.csv"))
  ####
  
  # coverage
  wt.codon=df$wildtype
  wt.codon.counts=unlist(sapply(1:nrow(df),function(x) df[x,wt.codon[x]]))
  row.counts=rowSums(df[,3:ncol(df)],na.rm = T)
  
  ###########################
  ## get the col infor to know if it is syn or not and the wt aa
  col.aa=as.character(translate(DNAStringSet(colnames(select(df,-site,-wildtype))),no.init.codon = T))
  wt.aa.seq=res.aa$wildtype
  ### syn data
  df.syn=select(df,-site,-wildtype)
  
  for (r in 1:nrow(df.syn)){
    # r=1
    # NA anything not syn
    df.syn[r, which(col.aa!=wt.aa.seq[r])]=NA
    df.syn[r, df$wildtype[r]]=NA
    rm(r)
  }
  df.syn=data.frame(df[,1:2], wt.aa=wt.aa.seq,df.syn)
  
  ## same for NS
  df.ns=select(df,-site,-wildtype)
  for (r in 1:nrow(df.ns)){
    # r=1
    df.ns[r, which(col.aa==wt.aa.seq[r])]=NA
    rm(r)
  }
  
  
  # remove stops from data
  df.ns[, colnames(df.ns)%in% c("TAA","TGA","TAG")]=NA
  df.ns=data.frame(df[,1:2], wt.aa=wt.aa.seq, df.ns)
  rm(wt.aa.seq,col.aa,wt.codon)
  ####
  ##### stops rate
  stops=rowSums(df[,colnames(df)%in% c("TAA","TGA","TAG")])
  
  ## ###########################
  
  #### get position, letter infor
  letters=c("A","C","G","T")
  df.let=data.frame(matrix(data=0,ncol=4,nrow=4))
  rownames(df.let)=letters
  colnames(df.let)=letters
  rm(letters)
  
  ### pos matrix for where in codon are the mutations to collect data
  df.pos=data.frame(matrix(data=0,ncol=3,nrow = 1))
  colnames(df.pos)=1:3
  
  ### res_codon_muts to collect data
  df.cod.mut=data.frame(matrix(data=0,ncol=3,nrow = 1))
  colnames(df.cod.mut)=c("single","double","triple")
  
  # iterate
  mut.rows=which(apply(df[,-c(1:2)],1, function(x) sum(x>1))>1) # places with >=1 mutation
  df2=df%>% select(-1,-2)
  
  for (i in mut.rows){ # iterate on rows where have muts
    # i=8
    r=df2[i,which(df2[i,]>0)]  #### keep all columns where you have a mutation  
    
    nms=colnames(r) # get all codon names
    counts=as.numeric(r[1,]) ## get counts
    
    wt=which(r==max(r)) ## figure out which is wt
    wt.cod=nms[wt] # get wt codon
    wt.count=counts[wt] # get wt count
    
    c.wt=unlist(strsplit(wt.cod,"")) ## bases of wt cod
    ncols=1:ncol(r) 
    ncols=ncols[!ncols%in%wt] # keep columns not wt for rest
    
    for (j in ncols){
      m.c=unlist(strsplit(nms[j],""))
      m.count=counts[j] # count of times see mutant codon
      
      diffs= sapply(1:3, function(x) c.wt[x]!=m.c[x]) ## finds which cod pos has different letter
      
      df.cod.mut[1,sum(diffs)]= 
        df.cod.mut[1,sum(diffs)]+m.count ## adds counts to column corresponding to number of letter diffs  
      
      ### add counts to letter matrix 
      if ((c.wt[1]!=m.c[1])){ # if differ at first base
        df.let[c.wt[1],m.c[1]]=df.let[c.wt[1],m.c[1]]+m.count # find letter x/y coordin at put count
        df.pos[1,1]=df.pos[1,1]+m.count ## add counts to first pos mutations
      }
      if ((c.wt[2]!=m.c[2])){# as above for second
        df.let[c.wt[2],m.c[2]]=df.let[c.wt[2],m.c[2]]+m.count
        df.pos[1,2]=df.pos[1,2]+m.count
      }
      if ((c.wt[3]!=m.c[3])){#same for thrid position
        df.let[c.wt[3],m.c[3]]=df.let[c.wt[3],m.c[3]]+m.count
        df.pos[1,3]=df.pos[1,3]+m.count
      }
      rm(m.c,m.count,diffs)
    }
    rm(r,i,counts,nms,wt.cod,wt,wt.count,c.wt,ncols,j)
  }
  
  rm(mut.rows,df2)
  # pivot letters
  df.let=data.frame(wt=row.names(df.let), df.let)
  df.let=pivot_longer(df.let,A:T) 
  colnames(df.let)=c("wt","mut","count")
  dir.create(paste0(out.dir,"/mutation_pos_bases/"),recursive = T,showWarnings = F)
  write_csv(df.let,paste0(out.dir,"/mutation_pos_bases/",name,"_bases.csv"))
  write_csv(df.pos,paste0(out.dir,"/mutation_pos_bases/",name,"_pos.csv"))
  rm(df.let,df.pos)
  
  
  ## results of row cov, wt cov, all mut cov, syn.cov, ns.cov, stp cov
  site.data=bind_cols(coverage=row.counts,wt.counts=wt.codon.counts,
                      mut.count=(row.counts-wt.codon.counts),
                      s.count=rowSums(df.syn[,-c(1:3)],na.rm = T),
                      ns.count=rowSums(df.ns[,-c(1:3)],na.rm = T),
                      stp.count=stops,
                      aa.sampled=apply(res.aa[,-c(1:2)],1,function(x) sum(x>0))-1,
                      codons.sampled=apply(df[,-c(1:2)],1,function(x) sum(x>0))-1)
  
  rm(df.ns,df.syn,stops,row.counts,wt.codon.counts)
  dir.create(paste0(out.dir,"/site_data/"),recursive = T,showWarnings = F)
  write_csv(site.data,paste0(out.dir,"/site_data/",name,"_counts.csv"))
  
  ## rate data per pos
  rate.data=bind_cols(mut.rate=site.data$mut.count/site.data$coverage,
                      s.rate=site.data$s.count/site.data$coverage,
                      ns.rate=site.data$ns.count/site.data$coverage,
                      stop.rate=site.data$stp.count/site.data$coverage)
  
  write_csv(rate.data,paste0(out.dir,"/site_data/",name,"_rates.csv"))
  rm(rate.data)
  ####### global rates
  
  sums.df=rbind(coverage=sum(site.data$coverage,na.rm = T),
                min.coverage=min(site.data$coverage,na.rm = T),
                frct.codons=sum(site.data$codons.sampled,na.rm = T)/(63*nrow(df)),
                frct.aa=sum(site.data$aa.sampled,na.rm = T)/(20*nrow(df)),
                mut.rate=sum(site.data$mut.count,na.rm = T)/sum(site.data$coverage,na.rm = T),
                rate.single=sum(df.cod.mut$single,na.rm = T)/sum(site.data$coverage,na.rm = T),
                rate.double=sum(df.cod.mut$double,na.rm = T)/sum(site.data$coverage,na.rm = T),
                rate.triple=sum(df.cod.mut$triple,na.rm = T)/sum(site.data$coverage,na.rm = T),
                s.rate=sum(site.data$s.count,na.rm = T)/sum(site.data$coverage,na.rm = T),
                ns.rate=sum(site.data$ns.count,na.rm = T)/sum(site.data$coverage,na.rm = T),
                stop.rate=sum(site.data$stp.count,na.rm = T)/sum(site.data$coverage,na.rm = T))
  
  rm(df.cod.mut)
  colnames(sums.df)=name
  dir.create(paste0(out.dir,"/stats/"),recursive = T,showWarnings = F)
  write.csv(sums.df,paste0(out.dir,"/stats/",name,"_stats.csv"),row.names = T)
  rm(sums.df)
  
  
}
}

# june 22 2023, RG
combine_stats=function(dir.files="../results/analyze_codon_tables/",
                       out.dir="../results/analyze_codon_tables/summary/"){
  library(tidyverse)
  
  fs=list.files(dir.files,full.names = T)
  df=read_csv(fs[1])
  for (i in 2:length(fs)){
    # i=2
    df2=read_csv(fs[i])
    df=left_join(df,df2)
    rm(df2,i)
  }
  colnames(df)[1]="parameter"
  dir.create(out.dir,recursive = T)
  write_csv(df,paste0(out.dir,"/summary.stats.csv"))
  
}



## plot statistics ### added 19-7-23
plot_stats=function(stat.file="../results/sum/summary.stats.csv",
                    out.dir="../results/sum/",
                    clean.name.pattern=".csv" ## removes this from condition name
){
  library(tidyverse)
  stats=read_csv(stat.file) %>% pivot_longer(cols = !parameter) #%>% pivot_wider(names_from = parameter,values_from = value)
  stats$name=gsub(  clean.name.pattern,"",stats$name)
  empty=tapply(stats$value,stats$parameter,sum)
  empty=names(empty)[empty==0]
  if (length(empty)>0){stats=stats[stats$parameter!= empty,]}
  
  for (p in unique(stats$parameter)){
    # p="coverage"
    tmp=stats[stats$parameter==p,]
    g=ggplot(tmp,aes(x=name,y=value))+geom_point()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ggtitle(p)
    ggsave(plot = g,filename = paste0(out.dir,p,".pdf"),device="pdf",width = 8, height=5)
    
  }
  
}
