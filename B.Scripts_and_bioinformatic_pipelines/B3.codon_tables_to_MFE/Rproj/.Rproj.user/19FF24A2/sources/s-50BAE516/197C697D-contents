

# October 2023, BAR
get_enrich_ratio=function(condition.file="./ab_comparisons.csv",
                          data.dir="../results/filter_expected/codon_files/",
                          out.dir="../results/analyze/mfe/",
                          pseudo=1,
                          min.cov.pre.or.post=NA,
                          min.cov.pre=NA,
                          min.cov.pre.post=NA,
                          rm.stops=T,
                          filter.position=T,
                          start=1,
                          end=851,
                          correct.start.pos=F,
                          aa.pos.correct=1,
                          mark.deads=F){
  
  ### reads mut files and get MFE or differential selection by adding corrected
  ### pseudotyped by coverage
  ### min.cov is for the minimum counts to call any mutation in the PRE selection only
  ### and rm.pre.zero removes mutants where there are 0 counts in the pre selection
  ### rm. stops avoid analysis of stops. 
  
  
  # condition.file="./ab_comparisons.csv"
  # data.dir="../results/filter_expected/codon_files/"
  # out.dir="../results/analyze/mfe/"
  # pseudo=1
  # min.cov.pre.or.post=1
  # min.cov.pre=2
  # min.cov.pre.post=3
  # rm.pre.zero=T
  # rm.stops=T
  # filter.position=T
  # start=18
  # end=605
  # correct.start.pos=T
  # aa.pos.correct=847
  
  ## also, think about when to remove pre 0
  library(tidyverse)
  library(Biostrings)
  
  # create out dir, add a / just in case it is missing when typed in
  dir.create(paste0(out.dir,"/"),recursive = T,showWarnings = F)
  
  # read files paths, add a / just in case it is missing when typed in
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  # read comparison file to iterate on
  comp=read_csv(condition.file)
  
  ## iterate on each line of the comparison file and do rel enrich
  for (r in 1:nrow(comp)){
    # r=1
    pre.file=fs[grep(comp$pre[r],basename(fs))]
    post.file=fs[grep(comp$post[r],basename(fs))]
    nm=comp$name[r]
    
    ## gets genetic code so can translate
    gcode=as.data.frame(GENETIC_CODE)
    
    # read pre, long format
    pre=read_csv(pre.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,pre=value) 
    
    if(rm.stops==T){
      pre=pre%>% 
        filter(!wildtype%in%c("TAA","TGA","TAG"))# remove stops
    }
    
    # get WT aa and mut AA 
    pre$wtAA=sapply(pre$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    pre$mutAA=sapply(pre$mut,function(x) gcode[x,"GENETIC_CODE"])
    
    ###### same for post
    # read data, pivot longer
    post=read_csv(post.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,post=value)
    
    if(rm.stops==T){
      post=post%>% 
        filter(!wildtype%in%c("TAA","TGA","TAG"))# remove stops
    }
    
    # get wt and mut AA
    post$wtAA=sapply(post$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    post$mutAA=sapply(post$mut,function(x) gcode[x,"GENETIC_CODE"])
    rm(gcode)
    
    
    ## combine All ###########
    # combine count data
    res=full_join(pre,post,by = c("site", "wildtype", "mut", "wtAA", "mutAA")) #%>% filter(site%in%24:25)
    rm(pre,post)
    
    ## now aggregate by mut AA
    res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% 
      filter(mutAA!="*") %>% summarize(pre=sum(pre,na.rm = T),
                                       post=sum(post,na.rm = T))
    
    
    
    ### remove sites with 0 in pre selection
    ### now included in filter step below
    # if (rm.pre.zero==T){
    #   ## NA absent in plasmid so don't get fake values for pseudo only
    #   rm.pos=which(res$pre==0)
    #   res$pre[rm.pos]=NA
    #   res$post[rm.pos]=NA
    # }
    # 
    
    ######################## filter sites by min coverages#####################
    ### remove all sites that don't reach min count in at least
    ### one of the conditions
    if (is.numeric(min.cov.pre.or.post)){
      res[which(res$pre < min.cov.pre.or.post & res$post < min.cov.pre.or.post),c("pre","post")]=NA}
    
    ## in both condition
    if (is.numeric(min.cov.pre.post)){
      res[which(res$pre < min.cov.pre.post | res$post < min.cov.pre.post),c("pre","post")]=NA}
    
    ##  in the PRE ONLY condition
    if (is.numeric(min.cov.pre)){
      res[which(res$pre < min.cov.pre),c("pre","post")]=NA}
    
    if (mark.deads){
      dead.vec=res$post
      dead.vec[dead.vec==0]=-20
      }
    ##########################################################################
    
    # get coverage by grouping according to WT and site
    cov=res %>% group_by(site,wildtype,wtAA) %>% 
      summarize(pre.cov=sum(pre,na.rm=T),post.cov=sum(post,na.rm=T))
    
    
    ## get pseudo based on coverage total
    cov$pseudo.pre=sapply(1:nrow(cov),function(x) pseudo * cov$pre.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    cov$pseudo.post=sapply(1:nrow(cov),function(x) pseudo * cov$post.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    
    # ## now aggregate by mut AA###
    ##### keep this option if you want to calculate things based on codon and not AA
    # res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% 
    #   filter(mutAA!="*") %>% summarize(pre=sum(pre,na.rm = T),
    #                                    post=sum(post,na.rm = T))
    
    
    ################## pseudo #####################
    ## add pseudo to counts in res (AA/mut AND wts)
    for (s in unique(res$site)){ 
      #s=1
      res$pre[res$site==s]= 
        res$pre[res$site==s]+cov$pseudo.pre[cov$site==s]
      res$post[res$site==s]= 
        res$post[res$site==s]+cov$pseudo.post[cov$site==s]
      rm(s)
    }
    
    rm(cov)
    
    
    #######
    ### get preferences################################
    
    res$enrich.pre=NA
    res$enrich.post=NA
    
    # standardize to wt
    for (s in unique(res$site)){
      #  s=24
      res$enrich.pre[res$site==s]=res$pre[res$site==s]/res$pre[res$site==s &res$wtAA==res$mutAA]
      res$enrich.post[res$site==s]=res$post[res$site==s]/res$post[res$site==s &res$wtAA==res$mutAA]
      
      rm(s)
    }
    res$prefs=res$enrich.post/res$enrich.pre
    res$mfe=log2(res$prefs)
    
    res=select(res,site,wildtype,wtAA,mutAA,log2mfe=mfe)
    
    if(mark.deads){
      res[which(dead.vec %in% -20),"log2mfe"]=-20
    rm(dead.vec)
      }
    
    
    if (filter.position==T){
      res=filter(res,site%in%start:end)
    }
    
    ### changed by SVA
    if (correct.start.pos==T){
      res$site=(res$site-start)+aa.pos.correct
    }
    
    names(res)[names(res)=="log2mfe"]=nm
    write_csv(res,paste0(out.dir,nm,"_rel.enrich.csv"))
    rm(nm,r)}
}

###########################################################################

get_enrich_stops=function(condition.file="comparisons_p1dms.csv",
                          data.dir="../results/filter_single_mutants/codon_files/",
                          out.dir="../results/analyze/mfe/stops/",
                          pseudo=1,
                          min.cov.pre.or.post=5,
                          min.cov.pre=NA,
                          min.cov.pre.post=NA,
                          filter.position=T,
                          start=1,
                          end=851,
                          correct.start.pos=F,
                          aa.pos.correct=1){
  
  ### reads mut files and get MFE or differential selection by adding corrected
  ### pseudotyped by coverage
  ### min.cov is for the minimum counts to call any mutation in the PRE selection only
  ### and rm.pre.zero removes mutants where there are 0 counts in the pre selection
  ### rm. stops avoid analysis of stops. 
  
  
  # condition.file="comparisons_p1dms.csv"
  # data.dir="../results/filter_single_mutants/codon_files/"
  # out.dir="../results/analyze/mfe/stops/"
  # pseudo=1
  # min.cov.pre.or.post=5
  # min.cov.pre=NA
  # min.cov.pre.post=NA
  # filter.position=T
  # start=1
  # end=851
  # correct.start.pos=F
  # aa.pos.correct=1
  
  ## also, think about when to remove pre 0
  library(tidyverse)
  library(Biostrings)
  
  # create out dir, add a / just in case it is missing when typed in
  dir.create(paste0(out.dir,"/"),recursive = T,showWarnings = F)
  
  # read files paths, add a / just in case it is missing when typed in
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  # read comparison file to iterate on
  comp=read_csv(condition.file)
  
  ## iterate on each line of the comparison file and do rel enrich
  for (r in 1:nrow(comp)){
    # r=1
    pre.file=fs[grep(comp$pre[r],basename(fs))]
    post.file=fs[grep(comp$post[r],basename(fs))]
    nm=comp$name[r]
    
    ## gets genetic code so can translate
    gcode=as.data.frame(GENETIC_CODE)
    
    # read pre, long format
    pre=read_csv(pre.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,pre=value)
    
    # get WT aa and mut AA 
    pre$wtAA=sapply(pre$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    pre$mutAA=sapply(pre$mut,function(x) gcode[x,"GENETIC_CODE"])
    
    ###### same for post
    # read data, pivot longer
    post=read_csv(post.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,post=value)
    
    
    # get wt and mut AA
    post$wtAA=sapply(post$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    post$mutAA=sapply(post$mut,function(x) gcode[x,"GENETIC_CODE"])
    rm(gcode)
    
    
    ## combine All ###########
    # combine count data
    res=full_join(pre,post,by = c("site", "wildtype", "mut", "wtAA", "mutAA")) #%>% filter(site%in%24:25)
    rm(pre,post)
    
    ## now aggregate by mut AA
    res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% 
    summarize(pre=sum(pre,na.rm = T),
    post=sum(post,na.rm = T))
    
    
    
    ### remove sites with 0 in pre selection
    ### now included in filter step below
    # if (rm.pre.zero==T){
    #   ## NA absent in plasmid so don't get fake values for pseudo only
    #   rm.pos=which(res$pre==0)
    #   res$pre[rm.pos]=NA
    #   res$post[rm.pos]=NA
    # }
    # 
    
    ######################## filter sites by min coverages#####################
    ### remove all sites that don't reach min count in at least
    ### one of the conditions
    if (is.numeric(min.cov.pre.or.post)){
      res[which(res$pre < min.cov.pre.or.post & res$post < min.cov.pre.or.post),c("pre","post")]=NA}
    
    ## in both condition
    if (is.numeric(min.cov.pre.post)){
      res[which(res$pre < min.cov.pre.post | res$post < min.cov.pre.post),c("pre","post")]=NA}
    
    ##  in the PRE ONLY condition
    if (is.numeric(min.cov.pre)){
      res[which(res$pre < min.cov.pre),c("pre","post")]=NA}
    ##########################################################################
    
    # get coverage by grouping according to WT and site
    cov=res %>% group_by(site,wildtype,wtAA) %>% 
      summarize(pre.cov=sum(pre,na.rm=T),post.cov=sum(post,na.rm=T))
    
    
    ## get pseudo based on coverage total
    cov$pseudo.pre=sapply(1:nrow(cov),function(x) pseudo * cov$pre.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    cov$pseudo.post=sapply(1:nrow(cov),function(x) pseudo * cov$post.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    
    # ## now aggregate by mut AA###
    ##### keep this option if you want to calculate things based on codon and not AA
    # res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% 
    #   filter(mutAA!="*") %>% summarize(pre=sum(pre,na.rm = T),
    #                                    post=sum(post,na.rm = T))
    
    res=res%>% 
      filter(mutAA=="*"|wtAA==mutAA)# keep wt and stops
    
    
    ################## pseudo #####################
    ## add pseudo to counts in res (AA/mut AND wts)
    for (s in unique(res$site)){ 
      # s=5
      res$pre[res$site==s]= 
        res$pre[res$site==s]+cov$pseudo.pre[cov$site==s]
      res$post[res$site==s]= 
        res$post[res$site==s]+cov$pseudo.post[cov$site==s]
      rm(s)
    }
    
    rm(cov)
    #######
    ### get preferences################################
    
    res$enrich.pre=NA
    res$enrich.post=NA
    
    # standardize to wt
    for (s in unique(res$site)){
      #  s=24
      res$enrich.pre[res$site==s]=res$pre[res$site==s]/res$pre[res$site==s &res$wtAA==res$mutAA]
      res$enrich.post[res$site==s]=res$post[res$site==s]/res$post[res$site==s &res$wtAA==res$mutAA]
      
      rm(s)
    }
    
    res=res%>% 
      filter(!wtAA==mutAA)# keep wt and stops
    res$prefs=res$enrich.post/res$enrich.pre
    res$mfe=log2(res$prefs)
    
    res=select(res,site,wildtype,wtAA,mutAA,log2mfe=mfe)
    
    if (filter.position==T){
      res=filter(res,site%in%start:end)
    }
    
    ### changed by SVA
    if (correct.start.pos==T){
      res$site=(res$site-start)+aa.pos.correct
    }
    
    names(res)[names(res)=="log2mfe"]=nm
    write_csv(res,paste0(out.dir,nm,"_rel.enrich.csv"))
    rm(nm,r)}
}


# june 22 2023, RG
get_enrich_syn=function(condition.file="../p2dms_comparisons.csv",
                        data.dir="../results/filter_expected/syn/",
                        out.dir="../results/analyze/mfe/syn/",
                        pseudo=1,
                        min.cov.pre.or.post=NA,
                        min.cov.pre=1,
                        min.cov.pre.post=NA,
                        rm.stops=T,
                        rm.pre.zero=T,
                        filter.position=T,
                        start=18,
                        end=605,
                        correct.start.pos=T,
                        aa.pos.correct=847){
  
  
  library(tidyverse)
  library(Biostrings)
  
  # condition.file="../p2dms_comparisons.csv"
  # data.dir="../results/unfiltered/deletions/"
  # out.dir="../results/analyze/mfe/dels/"
  # pseudo=1
  # min.cov=1
  # rm.pre.zero=T
  # rm.stops=T
  # filter.position=T
  # start=18
  # end=605
  # correct.start.pos=T
  # aa.pos.correct=847
  
  
  # read comparison file to iterate on
  comp=read_csv(condition.file)
  
  # create out dir, add a / just in case it is missing when typed in
  dir.create(paste0(out.dir,"/"),recursive = T,showWarnings = F)
  
  # read files paths, add a / just in case it is missing when typed in
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  ## iterate on each line of the comparison file and do rel enrich
  for (r in 1:nrow(comp)){
    # r=4
    pre.file=fs[grep(comp$pre[r],basename(fs))]
    post.file=fs[grep(comp$post[r],basename(fs))]
    nm=comp$name[r]
    
    
    # # gets genetic code so can translate
    gcode=as.data.frame(GENETIC_CODE)
    
    # read pre, long format
    pre=read_csv(pre.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,pre=value) 
    
    if(rm.stops==T){
      pre=pre%>% 
        filter(!wildtype%in%c("TAA","TGA","TAG"))# remove stops
    }
    
    # get WT aa and mut AA 
    pre$wtAA=sapply(pre$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    pre$mutAA=sapply(pre$mut,function(x) gcode[x,"GENETIC_CODE"])
    
    ###### same for post
    # read data, pivot longer
    post=read_csv(post.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,post=value)
    
    if(rm.stops==T){
      post=post%>% 
        filter(!wildtype%in%c("TAA","TGA","TAG"))# remove stops
    }
    
    
    # get wt and mut AA
    post$wtAA=sapply(post$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    post$mutAA=sapply(post$mut,function(x) gcode[x,"GENETIC_CODE"])
    rm(gcode)
    
    ## combine All ###########
    # combine count data
    res=full_join(pre,post,by = c("site", "wildtype", "mut", "wtAA", "mutAA")) #%>% filter(site%in%24:25)
    rm(pre,post)
    
    ####### keep wt codon data as it may change when group with different AA
    ## combine wt data #####
    # get wt count by filtering on wt==mut
    wts=filter(res,wildtype==mut)
    
    ### remove sites with 0 in pre selection
    if (rm.pre.zero==T){
      ## NA absent in plasmid so don't get fake values for pseudo only
      res$pre[res$pre==0]=NA
      res$post[res$pre==0]=NA
    }
    
    ######################## filter sites by min coverages#####################
    ### remove all sites that don't reach min count in at least
    ### one of the conditions
    if (is.numeric(min.cov.pre.or.post)){
      res[which(res$pre < min.cov.pre.or.post & res$post < min.cov.pre.or.post),c("pre","post")]=NA}
    
    ## in both condition
    if (is.numeric(min.cov.pre.post)){
      res[which(res$pre < min.cov.pre.post | res$post < min.cov.pre.post),c("pre","post")]=NA}
    
    ##  in the PRE ONLY condition
    if (is.numeric(min.cov.pre)){
      res[which(res$pre < min.cov.pre),c("pre","post")]=NA}
    ##########################################################################
    
    
    ## get coverage after removing things from filter
    # get coverage by grouping according to WT and site
    cov=res %>% group_by(site,wildtype,wtAA) %>% 
      summarize(pre.cov=sum(pre,na.rm=T),post.cov=sum(post,na.rm=T))
    
    ## get pseudo based on coverage total
    cov$pseudo.pre=sapply(1:nrow(cov),function(x) pseudo * cov$pre.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    cov$pseudo.post=sapply(1:nrow(cov),function(x) pseudo * cov$post.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    
    ## now aggregate by mut AA
    res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% filter(mutAA!="*") %>% summarize(pre=sum(pre,na.rm = T),
                                                                                        post=sum(post,na.rm = T))
    ################## pseudo
    
    ## add pseudo to counts in res (AA/mut AND wts)
    for (s in unique(res$site)){ 
      # s=5
      res$pre[res$site==s]= 
        res$pre[res$site==s]+cov$pseudo.pre[cov$site==s]
      res$post[res$site==s]= 
        res$post[res$site==s]+cov$pseudo.post[cov$site==s]
      
      wts$pre[wts$site==s]= 
        wts$pre[wts$site==s]+cov$pseudo.pre[cov$site==s]
      
      wts$post[wts$site==s]= 
        wts$post[wts$site==s]+cov$pseudo.post[cov$site==s]
      
      rm(s)
    }
    
    rm(cov)
    
    #######
    ### get preferences
    
    res$enrich.pre=NA
    res$enrich.post=NA
    
    # standardize to wt
    for (s in unique(res$site)){
      #  s=24
      res$enrich.pre[res$site==s]=res$pre[res$site==s]/wts$pre[wts$site==s]
      res$enrich.post[res$site==s]=res$post[res$site==s]/wts$post[wts$site==s]
      
      rm(s)
    }
    rm(wts)
    res$prefs=res$enrich.post/res$enrich.pre
    res$mfe=log2(res$prefs)
    
    res=filter(res,wtAA==mutAA)
    res=select(res,site,wildtype,wtAA,mutAA,log2mfe=mfe)# prefs,
    
    if (filter.position==T){
      res=filter(res,site%in%start:end)
    }
    if (correct.start.pos==T){
      res$site=(res$site-start)+aa.pos.correct
    }
    
    names(res)[names(res)=="log2mfe"]=nm
    write_csv(res,paste0(out.dir,nm,"_rel.enrich.csv"))
    rm(nm,r)}
}

# june 22 2023, RG
get_enrich_del=function(condition.file="../p2dms_comparisons.csv",
                        data.dir="../results/unfiltered/deletions/",
                        out.dir="../results/analyze/mfe/dels/",
                        pseudo=1,
                        min.cov.pre.or.post=NA,
                        min.cov.pre=1,
                        min.cov.pre.post=NA,
                        filter.position=T,
                        start=18,
                        end=605,
                        correct.start.pos=T,
                        aa.pos.correct=847,
                        mark.deads=F){
  
  # out.dir="../results/analyze/mfe/dels/"
  # condition.file=comparison.file
  # pseudo=pseudo.count
  # min.cov.pre.or.post = min.cov.pre.or.post
  # min.cov.pre=min.cov.pre
  # filter.position=filter.positions
  # start=AA.start.mut
  # end=AA.end.mut
  # correct.start.pos=correct.position
  # aa.pos.correct=startAApos
  # mark.deads=T
  library(tidyverse)
  # read comparison file to iterate on
  comp=read_csv(condition.file)
  
  # create out dir, add a / just in case it is missing when typed in
  dir.create(paste0(out.dir,"/"),recursive = T,showWarnings = F)
  
  # read files paths, add a / just in case it is missing when typed in
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  ## iterate on each line of the comparison file and do rel enrich
  for (r in 1:nrow(comp)){
    # r=1
    pre.file=fs[grep(comp$pre[r],basename(fs))]
    post.file=fs[grep(comp$post[r],basename(fs))]
    nm=comp$name[r]
    
    pre=read_csv(pre.file)  %>% select(site,
                                       wildtype,
                                       n.ref.pre=coverage,n.del.pre=CNT)
    
    
    post=read_csv(post.file) %>% select(site,
                                        wildtype,
                                        n.ref.post=coverage,n.del.post=CNT)
    
    
    
    res=full_join(pre,post) %>% arrange(site)
    rm(pre,post)
    
    ## turn NA to 0 and then filter as desired
    res$n.del.pre[is.na(res$n.del.pre)]=0
    res$n.del.post[is.na(res$n.del.post)]=0
    
    ######################## filter sites by min coverage#####################
    ### remove all sites that don't reach min count in at least
    ### one of the conditions
    if (is.numeric(min.cov.pre.or.post)){
      res[which(res$n.del.pre < min.cov.pre.or.post & res$n.del.post < min.cov.pre.or.post),c("n.del.pre","n.del.post")]=NA}
    
    ## in both condition
    if (is.numeric(min.cov.pre.post)){
      res[which(res$n.del.pre < min.cov.pre.post | res$n.del.post < min.cov.pre.post),c("n.del.pre","n.del.post")]=NA}
    
    ##  in the PRE ONLY condition
    if (is.numeric(min.cov.pre)){
      res[which(res$n.del.pre < min.cov.pre),c("n.del.pre","n.del.post")]=NA}
    
    if (mark.deads){
      dead.vec=res$n.del.post
      dead.vec[dead.vec==0]=-20
    }
    ########################################################################
    
    ## get adj.pseudo
    adj.pseudo=pseudo* select(res,n.ref.pre,n.ref.post)/apply(select(res,n.ref.pre,n.ref.post),1,min)
    
    ## add pseudo
    # according to adj pseudo for each
    res$n.ref.pre= res$n.ref.pre+adj.pseudo$n.ref.pre
    res$n.del.pre= res$n.del.pre+adj.pseudo$n.ref.pre
    
    res$n.ref.post= res$n.ref.post+adj.pseudo$n.ref.post
    res$n.del.post= res$n.del.post+adj.pseudo$n.ref.post
    
    # get log 2
    
    res=res %>% mutate(log2mfe= log2 ((n.del.post/n.ref.post)/(n.del.pre/n.ref.pre)))
    
    if(mark.deads){
      res[which(dead.vec %in% -20),"log2mfe"]=-20
      rm(dead.vec)
    }
    
    if (filter.position==T){
      res=filter(res,site%in%start:end)
    }
    if (correct.start.pos==T){
      res$site=(res$site-start)+aa.pos.correct
    }
    
    res=select(res,site ,wildtype,log2mfe)
    names(res)[names(res)=="log2mfe"]=nm
    write_csv(res,paste0(out.dir,nm,"_rel.enrich.csv"))
    rm(nm,r)}
}

##############################################################

get_ratio_tables=function(condition.file="./comparisons_p1dms.csv",
                          data.dir="../results/filter_single_mutants/codon_files/",
                          out.dir="../results/analyze/ratios/",
                          pseudo=1,
                          min.cov.pre.or.post=NA,
                          min.cov.pre=NA,
                          min.cov.pre.post=NA,
                          rm.stops=T,
                          filter.position=T,
                          start=1,
                          end=851,
                          correct.start.pos=F,
                          aa.pos.correct=1){
  
  ### reads mut files and get ratio by adding corrected
  ### pseudotyped by coverage
  ### min.cov is for the minimum counts to call any mutation in the PRE selection only
  ### and rm.pre.zero removes mutants where there are 0 counts in the pre selection
  ### rm. stops avoid analysis of stops. 
  
  
  # condition.file="./ab_comparisons.csv"
  # data.dir="../results/filter_expected/codon_files/"
  # out.dir="../results/analyze/mfe/"
  # pseudo=1
  # min.cov.pre.or.post=1
  # min.cov.pre=2
  # min.cov.pre.post=3
  # rm.pre.zero=T
  # rm.stops=T
  # filter.position=T
  # start=18
  # end=605
  # correct.start.pos=T
  # aa.pos.correct=847
  
  ## also, think about when to remove pre 0
  library(tidyverse)
  library(Biostrings)
  
  # create out dir, add a / just in case it is missing when typed in
  dir.create(paste0(out.dir,"/"),recursive = T,showWarnings = F)
  
  # read files paths, add a / just in case it is missing when typed in
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  # read comparison file to iterate on
  comp=read_csv(condition.file)
  
  ## iterate on each line of the comparison file and do rel enrich
  for (r in 1:nrow(comp)){
    # r=1
    pre.file=fs[grep(comp$pre[r],basename(fs))]
    post.file=fs[grep(comp$post[r],basename(fs))]
    nm.post=comp$name[r]
    nm.pre=comp$pre[r]
    ## gets genetic code so can translate
    gcode=as.data.frame(GENETIC_CODE)
    
    # read pre, long format
    pre=read_csv(pre.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,pre=value) 
    
    if(rm.stops==T){
      pre=pre%>% 
        filter(!wildtype%in%c("TAA","TGA","TAG"))# remove stops
    }
    
    # get WT aa and mut AA 
    pre$wtAA=sapply(pre$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    pre$mutAA=sapply(pre$mut,function(x) gcode[x,"GENETIC_CODE"])
    
    ###### same for post
    # read data, pivot longer
    post=read_csv(post.file)  %>% pivot_longer(AAA:TTT) %>% mutate(wtAA=NA,mutAA=NA,mut=name) %>% 
      select(1,2,mut,wtAA,mutAA,post=value)
    
    if(rm.stops==T){
      post=post%>% 
        filter(!wildtype%in%c("TAA","TGA","TAG"))# remove stops
    }
    
    # get wt and mut AA
    post$wtAA=sapply(post$wildtype,function(x) gcode[x,"GENETIC_CODE"])
    post$mutAA=sapply(post$mut,function(x) gcode[x,"GENETIC_CODE"])
    rm(gcode)
    
    
    ## combine All ###########
    # combine count data
    res=full_join(pre,post,by = c("site", "wildtype", "mut", "wtAA", "mutAA")) #%>% filter(site%in%24:25)
    rm(pre,post)
    
    ## now aggregate by mut AA
    res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% 
      filter(mutAA!="*") %>% summarize(pre=sum(pre,na.rm = T),
                                       post=sum(post,na.rm = T))
    
    
    
    ### remove sites with 0 in pre selection
    ### now included in filter step below
    # if (rm.pre.zero==T){
    #   ## NA absent in plasmid so don't get fake values for pseudo only
    #   rm.pos=which(res$pre==0)
    #   res$pre[rm.pos]=NA
    #   res$post[rm.pos]=NA
    # }
    # 
    
    ######################## filter sites by min coverages#####################
    ### remove all sites that don't reach min count in at least
    ### one of the conditions
    if (is.numeric(min.cov.pre.or.post)){
      res[which(res$pre < min.cov.pre.or.post & res$post < min.cov.pre.or.post),c("pre","post")]=NA}
    
    ## in both condition
    if (is.numeric(min.cov.pre.post)){
      res[which(res$pre < min.cov.pre.post | res$post < min.cov.pre.post),c("pre","post")]=NA}
    
    ##  in the PRE ONLY condition
    if (is.numeric(min.cov.pre)){
      res[which(res$pre < min.cov.pre),c("pre","post")]=NA}
    ##########################################################################
    
    # get coverage by grouping according to WT and site
    cov=res %>% group_by(site,wildtype,wtAA) %>% 
      summarize(pre.cov=sum(pre,na.rm=T),post.cov=sum(post,na.rm=T))
    
    
    ## get pseudo based on coverage total
    cov$pseudo.pre=sapply(1:nrow(cov),function(x) pseudo * cov$pre.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    cov$pseudo.post=sapply(1:nrow(cov),function(x) pseudo * cov$post.cov[x]/min(cov[x,c("pre.cov","post.cov")]))
    
    # ## now aggregate by mut AA###
    ##### keep this option if you want to calculate things based on codon and not AA
    # res=res %>% group_by(site,wildtype,wtAA,mutAA) %>% 
    #   filter(mutAA!="*") %>% summarize(pre=sum(pre,na.rm = T),
    #                                    post=sum(post,na.rm = T))
    
    
    ################## pseudo #####################
    ## add pseudo to counts in res (AA/mut AND wts)
    for (s in unique(res$site)){ 
      # s=5
      res$pre[res$site==s]= 
        res$pre[res$site==s]+cov$pseudo.pre[cov$site==s]
      res$post[res$site==s]= 
        res$post[res$site==s]+cov$pseudo.post[cov$site==s]
      rm(s)
    }
    
    rm(cov)
    #######
    ### get preferences################################
    
    res$enrich.pre=NA
    res$enrich.post=NA
    
    # standardize to wt
    for (s in unique(res$site)){
      #  s=24
      res$enrich.pre[res$site==s]=res$pre[res$site==s]/res$pre[res$site==s &res$wtAA==res$mutAA]
      res$enrich.post[res$site==s]=res$post[res$site==s]/res$post[res$site==s &res$wtAA==res$mutAA]
      
      rm(s)
    }
    
    
    res=select(res,site,wildtype,wtAA,mutAA,enrich.pre,enrich.post)
    
    if (filter.position==T){
      res=filter(res,site%in%start:end)
    }
    
    ### changed by SVA
    if (correct.start.pos==T){
      res$site=(res$site-start)+aa.pos.correct
    }
    
    names(res)[names(res)=="enrich.pre"]=nm.pre
    names(res)[names(res)=="enrich.post"]=nm.post
    write_csv(res,paste0(out.dir,nm.post,"_ratios.csv"))
    rm(nm.pre, nm.post,r)}
}

###############################################################################
combine_ratio_tables=function(data.dir="../results/analyze/ratios",
                     out.dir="../results/analyze/ratios/"){
  
  
  library(tidyverse)
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  fs=fs[grep("_ratios.csv",basename(fs))]
  
  # iterate and join all data
  df=read_csv(fs[1])
  for (r in 2:length(fs)){
    # r=2
    df2=read_csv(fs[r])
    df=df %>% left_join(df2, by=c("site","wildtype","wtAA","mutAA"))
    rm(df2,r)
  }
  df=df %>% select(!ends_with(".y"))
  write.csv(df, file= paste0(out.dir,"ratios_all.csv"), row.names = F)
  
}

################################################################################
