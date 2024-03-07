########## combine and plot difsel per mut and per site#########################

#BAR FUNCTIONS
get_mutdifsel=function(data.dir="../results/analyze/mfe/deads_not_marked/",
                       condition.file="./comparisons_p1dms.csv",
                       out.dir="../results/combined_difsel/mutdifsel/",
                       group="per_passage"){
# 
# data.dir="../results/analyze/mfe/deads_not_marked/"
# condition.file="./comparisons_p1dms.csv"
# out.dir="../results/combined_difsel/mutdifsel/"
# group="per_passage"

  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  
  dir.create(paste0("../results_",region,"/combined_difsel/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  dir.create(paste0(out.dir,group,"/negdifsel/"))
  dir.create(paste0(out.dir,group,"/posdifsel/"))
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
    posdf2=df2
    negdf2=df2
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:mutAA,cond$name[cond$grp==g])
      
      
      #for each rep get the pos and neg sitedifsel
      posdf <- tmp %>%
        group_by(site, wtAA, mutAA) %>%
        summarise(
          across(2:(length(.)-3), 
                 ~ ifelse(sum(.[. > 0], na.rm = TRUE) == 0, NA, sum(.[. > 0], na.rm = TRUE)), 
                 .names = "posdifsel_{col}")
        )
      
      negdf <- tmp %>%
        group_by(site, wtAA, mutAA) %>%
        summarise(
          across(2:(length(.)-3), 
                 ~ ifelse(sum(.[. < 0], na.rm = TRUE) == 0, NA, sum(.[. < 0], na.rm = TRUE)), 
                 .names = "negdifsel_{col}")
        )
      
      #########################################################################
      #plot posdifsel corplots and cormatrix 
      if(ncol(posdf)>4){
        
        ####### PLOT CORRPLOT #######
        
        toplot=posdf %>% ungroup() %>% select(!site:mutAA)
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
        ggsave(plot = combined_plots, paste0("correl_plots_posdf_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/posdifsel/"), width=plot_width, height = plot_height)
        
        
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
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_posdf_",g, ".pdf"),device = "pdf",   
               path = paste0(out.dir,group,"/posdifsel/"), width=plot_width, height = plot_width)}
      
      #########################################################################
      #plot negdifsel corplots and cormatrix
      if(ncol(negdf)>4){
        
        ####### PLOT CORRPLOT #######
        
        toplot=negdf %>% ungroup() %>% select(!site:mutAA)
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
        ggsave(plot = combined_plots, paste0("correl_plots_negdf_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/negdifsel/"), width=plot_width, height = plot_height)
        
        
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
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_negdf_",g, ".pdf"),device = "pdf",   
               path = paste0(out.dir,group,"/negdifsel/"), width=plot_width, height = plot_width)}
      
      #########################################################################
      #### GET AND SAVE AVERAGE ####
      ###For posdifsel
      
      # #Temporal index to assign the row number to each site and adds it to the dataframe.
      tempIndex = as.integer(1:length(posdf$site))
      posdf$tempIndex=tempIndex
      # 
      # # Creates a list with the sum of NA in each row 
      posdf=posdf %>% ungroup()
      site.nas=apply(select(posdf,3:(ncol(posdf)-1)),1,function(x) sum(is.na(x)))
      # 
      # 
      # Gets the mean of each row
      posdf$posdifsel= apply(select(posdf,4:(ncol(posdf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,5:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      posdf$posdifsel[posdf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      posdf = select(posdf, -"tempIndex")
      
      rm(tempIndex)
      
      
      posdf2=left_join(posdf2,select(posdf,site:mutAA,!!g:= posdifsel))
      # rm(zeros,tmp,g)
      rm(posdf)
      
      ###For negdifsel
      
      # #Temporal index to assign the row number to each site and adds it to the dataframe.
      tempIndex = as.integer(1:length(negdf$site))
      negdf$tempIndex=tempIndex
      # 
      # # Creates a list with the sum of NA in each row 
      negdf=negdf %>% ungroup()
      site.nas=apply(select(negdf,3:(ncol(negdf)-1)),1,function(x) sum(is.na(x)))
      # 
      # 
      # Gets the mean of each row
      negdf$negdifsel= apply(select(negdf,4:(ncol(negdf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,5:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      negdf$negdifsel[negdf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      negdf = select(negdf, -"tempIndex")
      
      rm(tempIndex)
      
      
      negdf2=left_join(negdf2,select(negdf,site:mutAA,!!g:= negdifsel))
      # rm(zeros,tmp,g)
      rm(negdf)
      
    }
    ## now remake df with the collapsed data
    posdf=posdf2
    negdf=negdf2
  }
  
  
  ##save positive difsel
  write.csv(posdf, file= paste0(out.dir,group,"/posdifsel/posdifsel_",group,".csv"), row.names = F)
  
  
  #save negative difsel
  write.csv(negdf, file= paste0(out.dir,group,"/negdifsel/negdifsel_",group,".csv"), row.names = F)
  
}





#########################


get_sitedifsel=function(data.dir="../results/analyze/mfe",
                        condition.file="./comparisons_p1dms.csv",
                        out.dir="../results/combined_difsel/sitedifsel/",
                        group="per_passage",
                        mean.or.median="mean",
                        sd.times=2){

  # data.dir="../results/analyze/mfe/deads_not_marked/"
  # condition.file="./comparisons_p1dms.csv"
  # out.dir="../results/combined_difsel/sitedifsel/"
  # group="per_rep"
  # mean.or.median="mean"
  # sd.times=2
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  library(RColorBrewer)
  
  
  dir.create(paste0("../results_",region,"/combined_difsel/"))
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  dir.create(paste0(out.dir,group,"/negdifsel/"))
  dir.create(paste0(out.dir,group,"/posdifsel/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  #read protein info
  prot.info=read.csv(ref.file)%>% filter(genome_region==region)
  
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
    posdf2=df2
    negdf2=df2
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="L1_virus_passage1"
      # select columns by group
      tmp= select(df,site:mutAA,cond$name[cond$grp==g]) 
      sitedf = filter(tmp, !wtAA == mutAA)
      
      
      #for each rep get the pos and neg sitedifsel
      posdf=sitedf %>%
        group_by(site, wtAA) %>%
        summarise(
          across(3:(length(.)-2), 
                 ~ ifelse(sum(.[. > 0], na.rm = TRUE) == 0, NA, sum(.[. > 0], na.rm = TRUE)),
                 .names = "posdifsel_{col}"))
      
      negdf=sitedf %>%
        group_by(site, wtAA) %>%
        summarise(
          across(3:(length(.)-2), 
                 ~ ifelse(sum(.[. < 0], na.rm = TRUE) == 0, NA, sum(.[. < 0], na.rm = TRUE)),
                 .names = "negdifsel_{col}"))
      
      
      
      
      #plot posdifsel corplots and cormatrix 
      if(ncol(posdf)>3){
        
        ####### PLOT CORRPLOT #######
        toplot=posdf %>% ungroup() %>% select(!site:wtAA)
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
        ggsave(plot = combined_plots, paste0("correl_plots_posdf_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/posdifsel/"), width=plot_width, height = plot_height)
        
        
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
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_posdf_",g, ".pdf"),device = "pdf",   
               path = paste0(out.dir,group,"/posdifsel/"), width=plot_width, height = plot_width)}
      
      #plot negdifsel corplots and cormatrix
      if(ncol(negdf)>3){
        
        ####### PLOT CORRPLOT #######
        toplot=negdf %>% ungroup() %>% select(!site:wtAA)
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
        ggsave(plot = combined_plots, paste0("correl_plots_negdf_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/negdifsel/"), width=plot_width, height = plot_height)
        
        
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
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_negdf_",g, ".pdf"),device = "pdf",  
               path = paste0(out.dir,group,"/negdifsel/"), width=plot_width, height = plot_width)}
      
      
      ####### GET AND SAVE AVERAGE #######
      ###For posdifsel
      
      # #Temporal index to assign the row number to each site and adds it to the dataframe.
      tempIndex = as.integer(1:length(posdf$site))
      posdf$tempIndex=tempIndex
      # 
      # # Creates a list with the sum of NA in each row 
      posdf=posdf %>% ungroup()
      site.nas=apply(select(posdf,3:(ncol(posdf)-1)),1,function(x) sum(is.na(x)))
      # 
      # 
      # Gets the mean of each row
      posdf$posdifsel= apply(select(posdf,3:(ncol(posdf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,5:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      posdf$posdifsel[posdf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      posdf = select(posdf, -"tempIndex")
      
      rm(tempIndex)
      
      posdf2=left_join(posdf2,select(posdf,site:wtAA,!!g:= posdifsel))
      # rm(zeros,tmp,g)
      rm(posdf)
      
      ###For negdifsel
      
      # #Temporal index to assign the row number to each site and adds it to the dataframe.
      tempIndex = as.integer(1:length(negdf$site))
      negdf$tempIndex=tempIndex
      # 
      # # Creates a list with the sum of NA in each row 
      negdf=negdf %>% ungroup()
      site.nas=apply(select(negdf,3:(ncol(negdf)-1)),1,function(x) sum(is.na(x)))
      # 
      # 
      # Gets the mean of each row
      negdf$negdifsel= apply(select(negdf,3:(ncol(negdf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,5:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      negdf$negdifsel[negdf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      negdf = select(negdf, -"tempIndex")
      
      rm(tempIndex)
      
      negdf2=left_join(negdf2,select(negdf,site:wtAA,!!g:= negdifsel))
      # rm(zeros,tmp,g)
      rm(negdf)
      
      
    }
    ## now remake df with the collapsed data
    posdf=posdf2
    negdf=negdf2
  }
  
  proteins=unique(prot.info$protein)
  set2_palette <- brewer.pal(length(proteins), "Set2")
  
  ##positive difsel
  
  lineplot=posdf %>% pivot_longer(cols=c(3:length(posdf)), names_to="sample", values_to="posdifsel")
  
  #add cutoff to plot
  if(mean.or.median=="mean"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=mean(posdifsel, na.rm = TRUE) + sd.times * sd(posdifsel, na.rm = TRUE))
    
  }
  if(mean.or.median=="median"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=median(posdifsel, na.rm = TRUE) + sd.times * sd(posdifsel, na.rm = TRUE))
    
  }
  
  if(mean.or.median=="NA"){
    lineplot$cutoff=0
  }
  
  lineplot$posdifsel <- ifelse(is.na(lineplot$posdifsel), 0, lineplot$posdifsel) 
  posdifselplot<- ggplot(data = lineplot, aes(x = site, y=posdifsel)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = lineplot$cutoff), colour="darkgray",linetype="dashed", lwd=0.3)+
    annotate("text",
             x = max(lineplot$site-100),
             y = mean(lineplot$cutoff),
             label = paste0(mean.or.median,"+",sd.times,"*SD"),
             color = "black",
             size = 3)+
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
    scale_y_continuous(expand=c(0,0))+
    lapply(seq_along(proteins), function(i) {
    annotate("text",
             x = mean(c(min(prot.info$site[prot.info$protein == proteins[i]]), max(prot.info$site[prot.info$protein == proteins[i]]))),
             y = mean(c(min(lineplot$posdifsel, na.rm = TRUE), max(lineplot$posdifsel, na.rm = TRUE))),
             label = proteins[i],
             color = "black",
             size = 3)})+
    lapply(seq_along(proteins), function(i) {
      annotate("rect", 
               xmin = min(prot.info$site[prot.info$protein == proteins[i]]),
               xmax = max(prot.info$site[prot.info$protein == proteins[i]]),
               ymin = min(lineplot$posdifsel, na.rm = TRUE),
               ymax = max(lineplot$posdifsel, na.rm = TRUE),
               alpha = 0.3,
               fill = set2_palette[i])})
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(posdf)/2)
  ggsave(plot = posdifselplot, paste0("lineplot_posdifsel_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/posdifsel/"), width=plot_width,
         height=plot_height)
  
  write.csv(posdf, file= paste0(out.dir,group,"/posdifsel/posdifsel_",group,".csv"), row.names = F)
  
  
  # negative difsel
  
  lineplot=negdf %>% pivot_longer(cols=c(3:length(negdf)), names_to="sample", values_to="negdifsel")
  
  #add cutoff to plot
  if(mean.or.median=="mean"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=mean(negdifsel, na.rm = TRUE) - sd.times * sd(negdifsel, na.rm = TRUE))
    
  }
  if(mean.or.median=="median"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=median(negdifsel, na.rm = TRUE) - sd.times * sd(negdifsel, na.rm = TRUE))
    
  }
  
  if(mean.or.median=="NA"){
    posdifsel$cutoff=0
  }
  
  lineplot$negdifsel <- ifelse(is.na(lineplot$negdifsel), 0, lineplot$negdifsel) 
  negdifselplot<- ggplot(data = lineplot, aes(x = site, y=negdifsel)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = lineplot$cutoff), colour="darkgray",linetype="dashed", lwd=0.3)+
    annotate("text",
             x = max(lineplot$site-100),
             y = mean(lineplot$cutoff),
             label = paste0(mean.or.median,"+",sd.times,"*SD"),
             color = "black",
             size = 3)+
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
    scale_y_continuous(expand=c(0,0))+
    lapply(seq_along(proteins), function(i) {
      annotate("text",
               x = mean(c(min(prot.info$site[prot.info$protein == proteins[i]]), max(prot.info$site[prot.info$protein == proteins[i]]))),
               y = mean(c(min(lineplot$negdifsel, na.rm = TRUE), max(lineplot$negdifsel, na.rm = TRUE))),
               label = proteins[i],
               color = "black",
               size = 3)})+
    lapply(seq_along(proteins), function(i) {
      annotate("rect", 
               xmin = min(prot.info$site[prot.info$protein == proteins[i]]),
               xmax = max(prot.info$site[prot.info$protein == proteins[i]]),
               ymin = min(lineplot$negdifsel, na.rm = TRUE),
               ymax = max(lineplot$negdifsel, na.rm = TRUE),
               alpha = 0.3,
               fill = set2_palette[i])})
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(negdf)/2)
  ggsave(plot = negdifselplot, paste0("lineplot_negdifsel_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/negdifsel/"), width=plot_width,
         height=plot_height)
  
  write.csv(negdf, file= paste0(out.dir,group,"/negdifsel/negdifsel_",group,".csv"), row.names = F)
  
}

###################################################
#### BAR 18 feb 2024 ### get and combine del difsel

get_delsdifsel=function(data.dir="../results/analyze/mfe/dels",
                        condition.file="./comparisons_p1dms.csv",
                        out.dir="../results/combined_difsel/delsdifsel/",
                        group="per_passage",
                        mean.or.median="mean",
                        sd.times=2){
  
  # data.dir="../results/analyze/mfe/dels/"
  # condition.file=comparison.file
  # out.dir="../results/combined_difsel/delsdifsel/"
  # group="per_passage"
  # mean.or.median=mean.or.median.cutoff
  # sd.times=sd.times.cutoff
  
  library(tidyverse)
  library(gridExtra)
  library(ggpubr)
  library(ggcorrplot)
  library(viridis)
  library(ggh4x)
  library(RColorBrewer)
  
  
  dir.create("../results/combined_difsel/")
  dir.create(out.dir)
  dir.create(paste0(out.dir,group,"/"))
  dir.create(paste0(out.dir,group,"/negdifsel/"))
  dir.create(paste0(out.dir,group,"/posdifsel/"))
  # read condition file
  cond=read.csv(condition.file,na.strings="")
  # read files
  fs=list.files(paste0(data.dir,"/"),full.names = T)
  
  
  #read protein info
  prot.info=read.csv(ref.file)%>% filter(genome_region==region)
  
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
    # group="per_rep"
    colnm=paste0(group)
    cond$grp=cond[,which(colnames(cond)==group)]
    groups=unique(cond$grp[!is.na(cond$grp)])
    df2=df %>% group_by(site) %>% summarise()
    posdf2=df2
    negdf2=df2
    cond=cond %>% filter(!is.na(grp))
    for (n in 1:length(groups)){
      #n=1
      g=groups[n]
      # g="per_rep"
      # select columns by group
      tmp= select(df,site:wildtype,cond$name[cond$grp==g]) 
      sitedf = tmp
      
      
      #for each rep get the pos and neg sitedifsel
      posdf=sitedf %>%
        group_by(site, wildtype) %>%
        summarise(
          if (ncol(.) > 3) {  # Check if there are enough columns to subset
            across(
              3:(ncol(.) - 2),  # Adjusted range to dynamically calculate number of columns
              ~ ifelse(sum(.[. > 0], na.rm = TRUE) == 0, NA, sum(.[. > 0], na.rm = TRUE)),
              .names = "posdifsel_{col}"
            )
          } else {
            across(
              everything(),
              ~ ifelse(sum(.[. > 0], na.rm = TRUE) == 0, NA, sum(.[. > 0], na.rm = TRUE)),
              .names = "posdifsel_{col}"
            )
          }
        )
      
      negdf=sitedf %>%
        group_by(site, wildtype) %>%
        summarise(
          if (ncol(.) < 3) {  # Check if there are enough columns to subset
            across(
              3:(ncol(.) - 2),  # Adjusted range to dynamically calculate number of columns
              ~ ifelse(sum(.[. < 0], na.rm = TRUE) == 0, NA, sum(.[. < 0], na.rm = TRUE)),
              .names = "negdifsel_{col}"
            )
          } else {
            across(
              everything(),
              ~ ifelse(sum(.[. < 0], na.rm = TRUE) == 0, NA, sum(.[. < 0], na.rm = TRUE)),
              .names = "negdifsel_{col}"
            )
          }
        )
      
      
      #plot posdifsel corplots and cormatrix 
      if(ncol(posdf)>3){
        
        ####### PLOT CORRPLOT #######
        toplot=posdf %>% ungroup() %>% select(!site:wildtype)
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
        ggsave(plot = combined_plots, paste0("correl_plots_posdf_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/posdifsel/"), width=plot_width, height = plot_height)
        
        
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
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_posdf_",g, ".pdf"),device = "pdf",   
               path = paste0(out.dir,group,"/posdifsel/"), width=plot_width, height = plot_width)}
      
      #plot negdifsel corplots and cormatrix
      if(ncol(negdf)>3){
        
        ####### PLOT CORRPLOT #######
        toplot=negdf %>% ungroup() %>% select(!site:wildtype)
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
        ggsave(plot = combined_plots, paste0("correl_plots_negdf_",g, ".pdf"),device = "pdf", 
               path = paste0(out.dir,group,"/negdifsel/"), width=plot_width, height = plot_height)
        
        
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
        
        ggsave(plot = corrmatrix, paste0("correl_matrix_negdf_",g, ".pdf"),device = "pdf",  
               path = paste0(out.dir,group,"/negdifsel/"), width=plot_width, height = plot_width)}
      
      
      ####### GET AND SAVE AVERAGE #######
      ###For posdifsel
      
      # #Temporal index to assign the row number to each site and adds it to the dataframe.
      tempIndex = as.integer(1:length(posdf$site))
      posdf$tempIndex=tempIndex
      # 
      # # Creates a list with the sum of NA in each row 
      posdf=posdf %>% ungroup()
      site.nas=apply(select(posdf,3:(ncol(posdf)-1)),1,function(x) sum(is.na(x)))
      # 
      # 
      # Gets the mean of each row
      posdf$posdifsel= apply(select(posdf,3:(ncol(posdf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,3:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      posdf$posdifsel[posdf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      posdf = select(posdf, -"tempIndex")
      
      rm(tempIndex)
      
      posdf2=left_join(posdf2,select(posdf,site:wildtype,!!g:= posdifsel))
      # rm(zeros,tmp,g)
      rm(posdf)
      
      ###For negdifsel
      
      # #Temporal index to assign the row number to each site and adds it to the dataframe.
      tempIndex = as.integer(1:length(negdf$site))
      negdf$tempIndex=tempIndex
      # 
      # # Creates a list with the sum of NA in each row 
      negdf=negdf %>% ungroup()
      site.nas=apply(select(negdf,3:(ncol(negdf)-1)),1,function(x) sum(is.na(x)))
      # 
      # 
      # Gets the mean of each row
      negdf$negdifsel= apply(select(negdf,3:(ncol(negdf)-1)),1,function(x) mean(x, na.rm=T))
      
      #NA cutoff depending on number of conditions to merge. When we have 4,
      #we want the cutoff to be >=2, when we have 3 >=2, when we have 8 >=5.
      na.cutoff=floor(ncol(tmp[,3:(ncol(tmp))])/2)+1
      
      # For each row in which the sum of NA is greater than 1, change the avg value to 0
      negdf$negdifsel[negdf$tempIndex%in% which(site.nas>=na.cutoff)]=NA
      
      #Eliminates temIndex column
      negdf = select(negdf, -"tempIndex")
      
      rm(tempIndex)
      
      negdf2=left_join(negdf2,select(negdf,site:wildtype,!!g:= negdifsel))
      # rm(zeros,tmp,g)
      rm(negdf)
      
      
    }
    ## now remake df with the collapsed data
    posdf=posdf2
    negdf=negdf2
  }
  
  proteins=unique(prot.info$protein)
  set2_palette <- brewer.pal(length(proteins), "Set2")
  
  ##positive difsel
  
  lineplot=posdf %>% pivot_longer(cols=c(3:length(posdf)), names_to="sample", values_to="posdifsel")
  
  #add cutoff to plot
  if(mean.or.median=="mean"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=mean(posdifsel, na.rm = TRUE) + sd.times * sd(posdifsel, na.rm = TRUE))
    
  }
  if(mean.or.median=="median"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=median(posdifsel, na.rm = TRUE) + sd.times * sd(posdifsel, na.rm = TRUE))
    
  }
  
  if(mean.or.median=="NA"){
    lineplot$cutoff=0
  }
  
  lineplot$posdifsel <- ifelse(is.na(lineplot$posdifsel), 0, lineplot$posdifsel) 
  posdifselplot<- ggplot(data = lineplot, aes(x = site, y=posdifsel)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = lineplot$cutoff), colour="darkgray",linetype="dashed", lwd=0.3)+
    annotate("text",
             x = max(lineplot$site-100),
             y = mean(lineplot$cutoff),
             label = paste0(mean.or.median,"+",sd.times,"*SD"),
             color = "black",
             size = 3)+
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
    scale_y_continuous(expand=c(0,0))+
    lapply(seq_along(proteins), function(i) {
      annotate("text",
               x = mean(c(min(prot.info$site[prot.info$protein == proteins[i]]), max(prot.info$site[prot.info$protein == proteins[i]]))),
               y = mean(c(min(lineplot$posdifsel, na.rm = TRUE), max(lineplot$posdifsel, na.rm = TRUE))),
               label = proteins[i],
               color = "black",
               size = 3)})+
    lapply(seq_along(proteins), function(i) {
      annotate("rect", 
               xmin = min(prot.info$site[prot.info$protein == proteins[i]]),
               xmax = max(prot.info$site[prot.info$protein == proteins[i]]),
               ymin = min(lineplot$posdifsel, na.rm = TRUE),
               ymax = max(lineplot$posdifsel, na.rm = TRUE),
               alpha = 0.3,
               fill = set2_palette[i])})
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(posdf)/2)
  ggsave(plot = posdifselplot, paste0("lineplot_posdifsel_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/posdifsel/"), width=plot_width,
         height=plot_height)
  
  write.csv(posdf, file= paste0(out.dir,group,"/posdifsel/posdifsel_",group,".csv"), row.names = F)
  
  
  # negative difsel
  
  lineplot=negdf %>% pivot_longer(cols=c(3:length(negdf)), names_to="sample", values_to="negdifsel")
  
  #add cutoff to plot
  if(mean.or.median=="mean"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=mean(negdifsel, na.rm = TRUE) - sd.times * sd(negdifsel, na.rm = TRUE))
    
  }
  if(mean.or.median=="median"){
    lineplot=lineplot %>% group_by(sample) %>% 
      mutate("cutoff"=median(negdifsel, na.rm = TRUE) - sd.times * sd(negdifsel, na.rm = TRUE))
    
  }
  
  if(mean.or.median=="NA"){
    posdifsel$cutoff=0
  }
  
  lineplot$negdifsel <- ifelse(is.na(lineplot$negdifsel), 0, lineplot$negdifsel) 
  negdifselplot<- ggplot(data = lineplot, aes(x = site, y=negdifsel)) +
    geom_line(linewidth=0.4)+
    geom_hline(aes(yintercept = lineplot$cutoff), colour="darkgray",linetype="dashed", lwd=0.3)+
    annotate("text",
             x = max(lineplot$site-100),
             y = mean(lineplot$cutoff),
             label = paste0(mean.or.median,"+",sd.times,"*SD"),
             color = "black",
             size = 3)+
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
    scale_y_continuous(expand=c(0,0))+
    lapply(seq_along(proteins), function(i) {
      annotate("text",
               x = mean(c(min(prot.info$site[prot.info$protein == proteins[i]]), max(prot.info$site[prot.info$protein == proteins[i]]))),
               y = mean(c(min(lineplot$negdifsel, na.rm = TRUE), max(lineplot$negdifsel, na.rm = TRUE))),
               label = proteins[i],
               color = "black",
               size = 3)})+
    lapply(seq_along(proteins), function(i) {
      annotate("rect", 
               xmin = min(prot.info$site[prot.info$protein == proteins[i]]),
               xmax = max(prot.info$site[prot.info$protein == proteins[i]]),
               ymin = min(lineplot$negdifsel, na.rm = TRUE),
               ymax = max(lineplot$negdifsel, na.rm = TRUE),
               alpha = 0.3,
               fill = set2_palette[i])})
  
  plot_width = 3 +  2*2
  plot_height = 3+2 * (length(negdf)/2)
  ggsave(plot = negdifselplot, paste0("lineplot_negdifsel_",group, ".pdf"),
         device = "pdf", path = paste0(out.dir,group,"/negdifsel/"), width=plot_width,
         height=plot_height)
  
  write.csv(negdf, file= paste0(out.dir,group,"/negdifsel/negdifsel_",group,".csv"), row.names = F)
  
}






