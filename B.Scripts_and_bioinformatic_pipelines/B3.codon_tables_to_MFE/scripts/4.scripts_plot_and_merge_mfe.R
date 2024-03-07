########## combine and plot MFEs per mut and per site###########################

#BAR FUNCTIONS
combine_mutmfe=function(data.dir="../results/analyze/mfe",
                        condition.file="./comparisons_p1dms.csv",
                        out.dir="../results/combined_mfe/mutmfe/",
                        out.name="all",
                        group="per_passage"){
  
  # data.dir=paste0("../results_",region,"/analyze/mfe/deads_not_marked/")
  # condition.file=comparison.file
  # out.dir=paste0("../results_",region,"/combined_mfe/mutmfe/")
  # group="per_region"

  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  
  dir.create(paste0("../results_",region,"/combined_mfe/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
 
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=2
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=select(df, site:mutAA)
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:mutAA,cond$name[cond$grp==g]) 
      # find length of 0 to not average single values
      # zeros=apply(select(tmp,-(site:mutAA)),1,function(x) sum(x<1e-5 &x > 0))
      # get mean
      if(ncol(tmp)>5){
        

      ####### PLOT CORRPLOT #######
        
      toplot=select(tmp,!site:mutAA)
      col_names <- colnames(toplot)
      combinations <- combn(col_names, 2)
      num_combinations <- ncol(combinations)
      plots_list <- list()
      
      # Loops through each combination and plot the correlation
      for (i in 1:num_combinations) {
        #i=1
        col1 <- combinations[1, i]
        col2 <- combinations[2, i]
        # Create the plot using ggscatter
        plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                          add = "reg.line", conf.int = TRUE,
                          cor.coef = TRUE, cor.method = "pearson")
        
        # Adds all the plots to the empty list 
        plots_list[[i]] <- plot
        }
      combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
      plot_width = 3 +  2*3
      plot_height = 3+2 * (num_combinations/3)
      
      ggsave(plot = combined_plots, paste0("correl_plots_mutmfe_",g, ".pdf"),device = "pdf", 
             path = paste0(out.dir,group,"/"), width=plot_width, height = plot_height)
      
      ###### PLOT CORRMATRIX ######

      col=rev(mako(10))
      corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
      corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
        scale_fill_gradientn(limit = c(0:1),colors=col)+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(),
              axis.ticks = element_line()
        )
      
      ggsave(plot = corrmatrix, paste0("correl_matrix_mutmfe_",g, ".pdf"),device = "pdf",  
             path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)}
      
      ######## GET AND SAVE AVERAGE ########

     #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(tmp$site))
      mutdf = tmp %>% cbind(tempIndex)
      
      # Creates a list with the sum of NA in each row 
      mut.nas=apply(select(tmp,5:(ncol(tmp))),1,function(x) sum(is.na(x)))
      
      
      # Gets the mean of each row 
      mutdf$avg= apply(select(mutdf,5:(ncol(tmp))),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(mutdf[,3:(ncol(mutdf))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      mutdf$avg[mutdf$tempIndex%in% which(mut.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      mutdf = select(mutdf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(mutdf,site:mutAA,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(tmp,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  ## now analyze by analyses column to get correlations
  
  ### Modified by SVA swince the function is not workinf properly when sourced. Check 
  
  write.csv(df, file= paste0(out.dir,group,"/mutmfe_",group,".csv"), row.names = F)
  
}


###############################################################################


combine_mutmfe_deads_marked=function(data.dir="../results/analyze/mfe",
                        condition.file="./comparisons_p1dms.csv",
                        out.dir="../results/combined_mfe/mutmfe/",
                        out.name="all",
                        group="per_passage"){
  # 
  # data.dir="../results/analyze/mfe/deads_marked/"
  # condition.file=comparison.file
  # out.dir="../results/combined_mfe_deads_marked/"
  # group="per_passage"
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  
  dir.create(paste0("../results_",region,"/combined_mfe_deads_marked/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=1
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=select(df, site:mutAA)
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=2
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:mutAA,cond$name[cond$grp==g]) 
      # find length of 0 to not average single values
      # zeros=apply(select(tmp,-(site:mutAA)),1,function(x) sum(x<1e-5 &x > 0))
      # get mean
      if(ncol(tmp)>5){
        
        ####### PLOT CORRPLOT #######
        
        toplot=select(tmp,!site:mutAA)
        toplot <- mutate_all(toplot, ~ifelse(. == -20, NA, .))
        col_names <- colnames(toplot)
        combinations <- combn(col_names, 2)
        num_combinations <- ncol(combinations)
        plots_list <- list()
        
        # Loops through each combination and plot the correlation
        for (i in 1:num_combinations) {
          col1 <- combinations[1, i]
          col2 <- combinations[2, i]
          # Create the plot using ggscatter
          plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson")
          
          # Adds all the plots to the empty list 
          plots_list[[i]] <- plot
        }
        combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
        plot_width = 3 +  2*3
        plot_height = 3+2 * (num_combinations/3)
        ggsave(plot = combined_plots, paste0("correl_plots_mutmfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"),  width=plot_width, height = plot_height)
        
        
        ####### PLOT CORRMATRIX #######
        
        col=rev(mako(10))
        corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
        corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
          scale_fill_gradientn(limit = c(0:1),colors=col)+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                axis.ticks = element_line()
          )
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_mutmfe_",g, ".pdf"),device = "pdf",   
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)
      }
      
      ####### GET AND SAVE AVERAGE #######
      
      #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(tmp$site))
      mutdf = tmp %>% cbind(tempIndex)
      
      # Creates a list with the sum of NA in each row 
      mut.nas=apply(select(tmp,5:(ncol(tmp))),1,function(x) sum(is.na(x)))
      
      # Creates a list with the sum of -20s in each row 
      mut.deads=apply(select(tmp,5:(ncol(tmp))),1,function(x) sum(x==-20))
      
      # Remove -20s and do average
      mutdf[mutdf==-20]<-NA
      
      # Gets the mean of each row excluding -20s
      mutdf$avg= apply(select(mutdf,5:(ncol(tmp))),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(mutdf[,5:(ncol(mutdf))])/2)+1
      
      # For each row in which the sum of NA is greater than cutoff, change the avg value to NA
      mutdf$avg[mutdf$tempIndex%in% which(mut.nas>=na.cutoff)]=NA
      
      # #deads cutoff depending on number of conditions to merge. When we have 4,
      # #we want the cutoff to be =>3, when we have 3 =>2, when we have 8 =>5 .
      deads.cutoff=floor(ncol(tmp[,5:(ncol(tmp))])/2)+1
      mutdf$avg[mutdf$tempIndex%in% which(mut.deads>=deads.cutoff)]=-20
      # 
      # # If a row contains all -20s, put a -20.
      # mutdf$avg[mutdf$tempIndex%in% which(mut.deads==ncol(tmp[,5:(ncol(tmp))]))]=-20
      
      #Eliminates temIndex column
      mutdf = select(mutdf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(mutdf,site:mutAA,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(tmp,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  ## now analyze by analyses column to get correlations
  
  ### Modified by SVA swince the function is not workinf properly when sourced. Check 
  
  write.csv(df, file= paste0(out.dir,group,"/mutmfe_",group,".csv"), row.names = F)
  
}

#########################


combine_sitemfe=function(data.dir="../results/analyze/mfe",
                         condition.file="./comparisons_p1dms.csv",
                         out.dir="../results/combined_mfe/sitemfe/",
                         group="per_passage"){
  
  
  # data.dir="../results/analyze/mfe"
  # condition.file="./comparisons_p1dms.csv"
  # out.dir="../results/combined_mfe/sitemfe/"
  # group="per_passage"
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  
  dir.create(paste0("../results_",region,"/combined_mfe/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=1
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=df %>% group_by(site) %>% summarise()
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:mutAA,cond$name[cond$grp==g]) 
      sitedf = filter(tmp, !wtAA == mutAA)
      sitedf=sitedf %>% group_by(site, wtAA) %>%  summarize(across(3:(length(tmp)-2),~ mean(.x,na.rm=T)))
      
      if(ncol(sitedf)>3){
        
        ####### PLOT CORRPLOT #######
        toplot=sitedf %>% ungroup() %>% select(!site:wtAA)
        col_names <- colnames(toplot)
        combinations <- combn(col_names, 2)
        num_combinations <- ncol(combinations)
        plots_list <- list()
        
        # Loops through each combination and plot the correlation
        for (i in 1:num_combinations) {
          col1 <- combinations[1, i]
          col2 <- combinations[2, i]
          # Create the plot using ggscatter
          plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson")
          
          # Adds all the plots to the empty list 
          plots_list[[i]] <- plot
        }
        combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
        plot_width = 3 +  2*3
        plot_height = 3+2 * (num_combinations/3)
        ggsave(plot = combined_plots, paste0("correl_plots_sitemfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_height)
        
        
        ####### PLOT CORRMATRIX #######
        col=rev(mako(10))
        corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
        corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
          scale_fill_gradientn(limit = c(0:1),colors=col)+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                axis.ticks = element_line()
          )
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_sitemfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)
      }
      
      ####### GET AND SAVE AVERAGE #######
      
      #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(sitedf$site))
      sitedf$tempIndex=tempIndex
      
      # Creates a list with the sum of NA in each row 
      sitedf=sitedf %>% ungroup()
      site.nas=apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) sum(is.na(x)))
      
      
      # Gets the mean of each row 
      sitedf$avg= apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(sitedf[,3:(ncol(sitedf))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      sitedf$avg[sitedf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      sitedf = select(sitedf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(sitedf,site:wtAA,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(sitedf,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  lineplot=df %>% pivot_longer(cols=c(3:length(df)), names_to="sample", values_to="sitemfe")
  sitemfeplot<- ggplot(data = lineplot, aes(x = site, y=sitemfe)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = 0), colour="darkgray",linetype="dashed", lwd=0.3)+
    facet_wrap2(.~sample,  axes = "all",
                #remove_labels = "x", 
                ncol = 1)+
    ylab("site mfe") +
    theme(text = element_text(size = 12),
          panel.background = element_blank(),axis.line = element_line(colour = "black", size=0.4),
          axis.title.y = element_text(margin = margin(r = 10)), axis.title.x = element_text(margin = margin(t = 5)),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          strip.background =element_rect(fill="white"))+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(df)/2)
  ggsave(plot = sitemfeplot, paste0("lineplot_sitemfe_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/"), width=plot_width,
         height=plot_height)
  
  write.csv(df, file= paste0(out.dir,group,"/sitemfe_",group,".csv"), row.names = F)
  
}

###################################################

combine_delsmfe=function(data.dir="../results/analyze/mfe/dels/",
                         condition.file="./comparisons_p1dms.csv",
                         out.dir="../results/combined_mfe/delsmfe/",
                         group="per_rep"){
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  
  dir.create(paste0("../results_",region,"/combined_mfe/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=2
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=df %>% group_by(site) %>% summarise()
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:wildtype,cond$name[cond$grp==g]) 
      sitedf=tmp
      
      if(ncol(sitedf)>3){
        
        ####### PLOT CORRPLOT####
        toplot=sitedf %>% ungroup() %>% select(!site:wildtype)
        col_names <- colnames(toplot)
        combinations <- combn(col_names, 2)
        num_combinations <- ncol(combinations)
        plots_list <- list()
        
        # Loops through each combination and plot the correlation
        for (i in 1:num_combinations) {
          col1 <- combinations[1, i]
          col2 <- combinations[2, i]
          # Create the plot using ggscatter
          plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson")
          
          # Adds all the plots to the empty list 
          plots_list[[i]] <- plot
        }
        combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
        plot_width = 3 +  2*3
        plot_height = 3+2 * (num_combinations/3)
        ggsave(plot = combined_plots, paste0("correl_plots_delsmfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_height)
        
        
        ####### PLOT CORRMATRIX####
        col=rev(mako(10))
        corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
        corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
          scale_fill_gradientn(limit = c(0:1),colors=col)+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                axis.ticks = element_line()
          )
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_delsmfe_",g, ".pdf"),device = "pdf",  
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)}
      
      #### GET AND SAVE AVERAGE ####
      
      #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(sitedf$site))
      sitedf$tempIndex=tempIndex
      
      # Creates a list with the sum of NA in each row 
      sitedf=sitedf %>% ungroup()
      site.nas=apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) sum(is.na(x)))
      
      
      # Gets the mean of each row 
      sitedf$avg= apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(sitedf[,3:(ncol(sitedf))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      sitedf$avg[sitedf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      sitedf = select(sitedf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(sitedf,site:wildtype,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(sitedf,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  lineplot=df %>% pivot_longer(cols=c(3:length(df)), names_to="sample", values_to="delsmfe")
  delsmfeplot<- ggplot(data = lineplot, aes(x = site, y=delsmfe)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = 0), colour="darkgray",linetype="dashed", lwd=0.3)+
    facet_wrap2(.~sample,  axes = "all",
                #remove_labels = "x", 
                ncol = 2)+
    ylab("site mfe") +
    theme(text = element_text(size = 12),
          panel.background = element_blank(),axis.line = element_line(colour = "black", size=0.4),
          axis.title.y = element_text(margin = margin(r = 10)), axis.title.x = element_text(margin = margin(t = 5)),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          strip.background =element_rect(fill="white"))+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(df)/2)
  
  ggsave(plot = delsmfeplot, paste0("lineplot_delsmfe_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/"), width=plot_width,
         height=plot_height)
  
  write.csv(df, file= paste0(out.dir,group,"/delsmfe_",group,".csv"), row.names = F)
  
}
#################################################################################
combine_delsmfe_deads_marked=function(data.dir="../results/analyze/mfe/dels_deads_marked/",
                         condition.file="./comparisons_p1dms.csv",
                         out.dir="../results/combined_mfe/delsmfe_deads_marked/",
                         group="per_rep"){
  # data.dir="../results/analyze/mfe/dels_deads_marked/"
  # condition.file=comparison.file
  # out.dir="../results/combined_mfe/delsmfe/"
  # group="per_passage"
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  
  dir.create(paste0("../results_",region,"/combined_mfe/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=2
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=df %>% group_by(site) %>% summarise()
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:wildtype,cond$name[cond$grp==g]) 
      sitedf=tmp
      
      if(ncol(sitedf)>3){
        
        ####### PLOT CORRPLOT####
        toplot=sitedf %>% ungroup() %>% select(!site:wildtype)
        toplot <- mutate_all(toplot, ~ifelse(. == -20, NA, .))
        col_names <- colnames(toplot)
        combinations <- combn(col_names, 2)
        num_combinations <- ncol(combinations)
        plots_list <- list()
        
        # Loops through each combination and plot the correlation
        for (i in 1:num_combinations) {
          col1 <- combinations[1, i]
          col2 <- combinations[2, i]
          # Create the plot using ggscatter
          plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson")
          
          # Adds all the plots to the empty list 
          plots_list[[i]] <- plot
        }
        combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
        plot_width = 3 +  2*3
        plot_height = 3+2 * (num_combinations/3)
        ggsave(plot = combined_plots, paste0("correl_plots_delsmfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_height)
        
        
        ####### PLOT CORRMATRIX####
        col=rev(mako(10))
        corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
        corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
          scale_fill_gradientn(limit = c(0:1),colors=col)+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                axis.ticks = element_line()
          )
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_delsmfe_",g, ".pdf"),device = "pdf",  
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)}
      
      #### GET AND SAVE AVERAGE ####
      
      #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(sitedf$site))
      sitedf$tempIndex=tempIndex
      
      # Creates a list with the sum of NA in each row 
      sitedf=sitedf %>% ungroup()
      site.nas=apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) sum(is.na(x)))
      
      # Creates a list with the sum of -20s in each row 
      site.deads=apply(select(tmp,3:(ncol(tmp))),1,function(x) sum(x==-20))
      
      # Remove deads and get the mean of each row 
      sitedf[sitedf==-20]<-NA
      sitedf$avg= apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,3:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      sitedf$avg[sitedf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      # #deads cutoff depending on number of conditions to merge. When we have 4,
      # #we want the cutoff to be =>3, when we have 3 =>2, when we have 8 =>5 .
      deads.cutoff=floor(ncol(tmp[,3:(ncol(tmp))])/2)+1
      sitedf$avg[sitedf$tempIndex%in% which(site.deads>=deads.cutoff)]=-20
      
      #Eliminates temIndex column
      sitedf = select(sitedf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(sitedf,site:wildtype,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(sitedf,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  lineplot=df %>% pivot_longer(cols=c(3:length(df)), names_to="sample", values_to="delsmfe")
  delsmfeplot<- ggplot(data = lineplot, aes(x = site, y=delsmfe)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = 0), colour="darkgray",linetype="dashed", lwd=0.3)+
    facet_wrap2(.~sample,  axes = "all",
                #remove_labels = "x", 
                ncol = 2)+
    ylab("site mfe") +
    theme(text = element_text(size = 12),
          panel.background = element_blank(),axis.line = element_line(colour = "black", size=0.4),
          axis.title.y = element_text(margin = margin(r = 10)), axis.title.x = element_text(margin = margin(t = 5)),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          strip.background =element_rect(fill="white"))+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(df)/2)
  
  ggsave(plot = delsmfeplot, paste0("lineplot_delsmfe_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/"), width=plot_width,
         height=plot_height)
  
  write.csv(df, file= paste0(out.dir,group,"/delsmfe_",group,".csv"), row.names = F)
  
}
#################################################################################


combine_synmfe=function(data.dir="../results/analyze/mfe/syn/",
                        condition.file="./comparisons_p1dms.csv",
                        out.dir="../results/combined_mfe/synmfe/",
                        group="per_rep"){
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  
  dir.create(paste0("../results_",region,"/combined_mfe/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=2
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=df %>% group_by(site) %>% summarise()
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:wildtype,cond$name[cond$grp==g]) 
      sitedf=tmp %>% group_by(site, wildtype) %>%  summarize(across(1:(length(tmp)-2),~ mean(.x,na.rm=T)))
      
      if(ncol(sitedf)>3){
        
        ####### PLOT CORRPLOT####
        toplot=sitedf %>% ungroup() %>% select(!site:wildtype)
        col_names <- colnames(toplot)
        combinations <- combn(col_names, 2)
        num_combinations <- ncol(combinations)
        plots_list <- list()
        
        # Loops through each combination and plot the correlation
        for (i in 1:num_combinations) {
          col1 <- combinations[1, i]
          col2 <- combinations[2, i]
          # Create the plot using ggscatter
          plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson")
          
          # Adds all the plots to the empty list 
          plots_list[[i]] <- plot
        }
        combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
        plot_width = 3 +  2*3
        plot_height = 3+2 * (num_combinations/3)
         ggsave(plot = combined_plots, paste0("correl_plots_synmfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_height)
        
        
        ####### PLOT CORRMATRIX####
        col=rev(mako(10))
        corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
        corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
          scale_fill_gradientn(limit = c(0:1),colors=col)+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                axis.ticks = element_line()
          )
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_synmfe_",g, ".pdf"),device = "pdf",  
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)}
      
      #### GET AND SAVE AVERAGE ####
      
      #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(sitedf$site))
      sitedf$tempIndex=tempIndex
      
      # Creates a list with the sum of NA in each row 
      sitedf=sitedf %>% ungroup()
      site.nas=apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) sum(is.na(x)))
      
      
      # Gets the mean of each row 
      sitedf$avg= apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(sitedf[,3:(ncol(sitedf))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      sitedf$avg[sitedf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      sitedf = select(sitedf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(sitedf,site:wildtype,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(sitedf,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  lineplot=df %>% pivot_longer(cols=c(3:length(df)), names_to="sample", values_to="synmfe")
  synmfeplot<- ggplot(data = lineplot, aes(x = site, y=synmfe)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = 0), colour="darkgray",linetype="dashed", lwd=0.3)+
    facet_wrap2(.~sample,  axes = "all",
                #remove_labels = "x", 
                ncol = 2)+
    ylab("site mfe") +
    theme(text = element_text(size = 12),
          panel.background = element_blank(),axis.line = element_line(colour = "black", size=0.4),
          axis.title.y = element_text(margin = margin(r = 10)), axis.title.x = element_text(margin = margin(t = 5)),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          strip.background =element_rect(fill="white"))+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(df)/2)
  ggsave(plot = synmfeplot, paste0("lineplot_synmfe_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/"), width=plot_width,
         height=plot_height)
  
  write.csv(df, file= paste0(out.dir,group,"/synmfe_",group,".csv"), row.names = F)
  
}



#####################################################################


combine_stopsmfe=function(data.dir="../results/analyze/mfe/stops/",
                        condition.file="./comparisons_p1dms.csv",
                        out.dir="../results/combined_mfe/stopsmfe/",
                        group="per_rep"){
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  
  dir.create(paste0("../results_",region,"/combined_mfe/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  
  # iterate and join all data
  df=read_csv(fs[grep(paste0(cond$name[1],"_rel.enrich.csv"),basename(fs))])
  for (r in 2:nrow(cond)){
    # r=2
    df2=read_csv(fs[grep(paste0(cond$name[r],"_rel.enrich.csv"),basename(fs))])
    df=full_join(df,df2)
    rm(df2,r)
  }
  
  
  ### if there is a grouping defined
  if (!is.na(group)){
    # read condition file
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=df %>% group_by(site) %>% summarise()
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:wildtype,cond$name[cond$grp==g]) 
      sitedf=tmp %>% group_by(site, wildtype) %>%  summarize(across(1:(length(tmp)-2),~ mean(.x,na.rm=T)))
      
      if(ncol(sitedf)>3){
        
        ####### PLOT CORRPLOT####
        toplot=sitedf %>% ungroup() %>% select(!site:wildtype)
        col_names <- colnames(toplot)
        combinations <- combn(col_names, 2)
        num_combinations <- ncol(combinations)
        plots_list <- list()
        
        # Loops through each combination and plot the correlation
        for (i in 1:num_combinations) {
          col1 <- combinations[1, i]
          col2 <- combinations[2, i]
          # Create the plot using ggscatter
          plot <- ggscatter(na.omit(toplot), x = col1, y = col2, 
                            add = "reg.line", conf.int = TRUE,
                            cor.coef = TRUE, cor.method = "pearson")
          
          # Adds all the plots to the empty list 
          plots_list[[i]] <- plot
        }
        combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
        plot_width = 3 +  3*2
        plot_height = 3+2 * (num_combinations/3)
        ggsave(plot = combined_plots, paste0("correl_plots_stopsmfe_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_height)
        
        
        ####### PLOT CORRMATRIX #######
        col=rev(mako(10))
        corr <- round(cor(toplot, use = "complete.obs",method = "pearson"), 2)
        corrmatrix <- ggcorrplot(corr,lab = TRUE,type ="upper",lab_col = "white", lab_size =4)+
          scale_fill_gradientn(limit = c(0:1),colors=col)+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(),
                axis.ticks = element_line()
          )
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_stopsmfe_",g, ".pdf"),device = "pdf",   
               path = paste0(out.dir,group,"/"), width=plot_width, height = plot_width)}
      
      ####### GET AND SAVE AVERAGE #######
      
      #Temporal index to assign the row number to each site and adds it to the dataframe.
      
      tempIndex = as.integer(1:length(sitedf$site))
      sitedf$tempIndex=tempIndex
      
      # Creates a list with the sum of NA in each row 
      sitedf=sitedf %>% ungroup()
      site.nas=apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) sum(is.na(x)))
      
      
      # Gets the mean of each row 
      sitedf$avg= apply(select(sitedf,3:(ncol(sitedf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(sitedf[,3:(ncol(sitedf))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      sitedf$avg[sitedf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      sitedf = select(sitedf, -"tempIndex")
      
      rm(tempIndex)
      
      df2=left_join(df2,select(sitedf,site:wildtype,!!g:= avg))
      # rm(zeros,tmp,g)
      rm(sitedf,g,n)
      
    }
    ## now remake df with the collapsed data
    df=df2
  }
  
  lineplot=df %>% pivot_longer(cols=c(3:length(df)), names_to="sample", values_to="stopsmfe")
  stopsmfeplot<- ggplot(data = lineplot, aes(x = site, y=stopsmfe)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = 0), colour="darkgray",linetype="dashed", lwd=0.3)+
    facet_wrap2(.~sample,  axes = "all",
                #remove_labels = "x", 
                ncol = 2)+
    ylab("site mfe") +
    theme(text = element_text(size = 12),
          panel.background = element_blank(),axis.line = element_line(colour = "black", size=0.4),
          axis.title.y = element_text(margin = margin(r = 10)), axis.title.x = element_text(margin = margin(t = 5)),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          strip.background =element_rect(fill="white"))+
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(df)/2)
  ggsave(plot = stopsmfeplot, paste0("lineplot_stopsmfe_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/"),width=plot_width,
         height=plot_height)
  
  write.csv(df, file= paste0(out.dir,group,"/stopsmfe_",group,".csv"), row.names = F)
  
}
