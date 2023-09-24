#Derive functional group corrrealtions (group level rather than species-species)
#Jordan DiNardo

rm(list=ls())

#Develop generalized parrotfish cpue timeseries by region
library(data.table);library(tidyverse);library(scales);library(cowplot)

load(file="./model_output/regional_species_results_final.RData")
source("./functions/reef_functions.R") #data prep

regions <- c("Bonaire", "Key Largo", "Cozumel", "Little Cayman", "Roatan","Tortolla") #define region
pf_commonnames <- c("Redband","Redtail","Yellowtail (Redfin)","Stoplight","Rainbow","Midnight","Queen","Princess","Blue","Striped") #define sp common name
pf_speciesid <- c(152,154,153,155,151,148,150,149,144,156) #define species id
N_R <- length(pf_speciesid) #define number of species

#derive timeseries estimates
for (i in 1:length(regions)){
  pf_ts[[i]] <- ts_reef(pf_occs[[i]])
}

#create list to populate with model outputs (population and observation estimates)
x_trends <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){
  
  x_trends[[i]] <- estimate_trends(sp=pf_commonnames,region=regions[i],grp=pf_commonnames,params1=pf_param[[i]],TT=max(pf_occs[[i]]$year)-min(pf_occs[[i]]$year)+1,ts=pf_ts[[i]],n.groups=N_R)
}

x_trends_new <- x_trends

for (i in 1:length(regions)){ #for each region
  for(j in 1:length(pf_commonnames)){ #for each species
    for(k in 1:(max(pf_occs[[i]]$year)-min(pf_occs[[i]]$year)+1)){ #for each year
      
      x_trends_new[[i]][[j]][[k]]$species <-sort(pf_commonnames)[j] #create species var
      x_trends_new[[i]][[j]][[k]]$year <-(min(pf_occs[[i]]$year):max(pf_occs[[i]]$year))[k] #create year var
      x_trends_new[[i]][[j]][[k]]$region <- regions[i] #create region var
    }
    x_trends_new[[i]][[j]] <- rbindlist(x_trends_new[[i]][[j]]) #combine lists
  }
  
  x_trends_new[[i]] <- rbindlist(x_trends_new[[i]]) #combine lists
}

scrapers <- c("Princess","Striped","Queen","Blue") #define species to scrapers functional group
excavators <- c("Stoplight","Rainbow","Midnight") #define species to excavators functional group
browsers <- c("Redband", "Redtail","Redfin") #define species to browsers functional group

x_trends_new2 <- rbindlist(x_trends_new) #combine lists

#create lsit 
x_trends_new3 <- vector(mode="list",length=length(regions))

for (i in 1:length(regions)){
  
  x_trends_new3[[i]] <- x_trends_new[[i]]%>%
    #filter(year>1993)%>%
    mutate(
      species=ifelse(species=="Yellowtail (Redfin)","Redfin",species), #adjust name of species
      functional_group=NA, #create functional group var
      functional_group=ifelse(species %in% scrapers, "Scraper",functional_group), #match species to functional group
      functional_group=ifelse(species %in% excavators, "Excavator",functional_group), #match species to functional group
      functional_group=ifelse(species %in% browsers, "Browser",functional_group) #match species to functional group
    )
  
  #merge nsurv from pf_ts df to results
  x_trends_new3[[i]] <- merge(x_trends_new3[[i]],pf_ts[[i]][,c("year","n.survs")],by="year",all.x=T)
  
  
}

#derive regional aggregated relative abundance of pf by functional group
sums_regions <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){
  
  sums_regions[[i]] <- x_trends_new3[[i]]%>%
    group_by(region,functional_group,year,n.survs,iter)%>%
    summarise(
      sum_x=sum(lambda.x,na.rm=T),
      sum_y=sum(lambda.y,na.rm=T))%>%
    ungroup()%>%
    group_by(region,functional_group,year,n.survs)%>%
    summarize(
      med_x= median(sum_x,na.rm=T),
      med_y=median(sum_y,na.rm=T),
      cpue_lower=quantile(sum_x,0.05),
      cpue_upper=quantile(sum_x,0.95))%>%
    ungroup()
}


#Create compare dfs to run correlation tests using cor.test 
pair_1 <- c("Scraper","Excavator","Scraper")
pair_2 <-c("Excavator","Browser","Browser")
pair1_df <- vector(mode="list",length=length(regions))
pair2_df <- vector(mode="list",length=length(regions))
compare_df <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){ #for each region
  
  for(j in 1:length(pair_1)){ #for each comparison
    
    pair1_df[[i]][[j]] <- sums_regions[[i]]%>%
      filter(functional_group==pair_1[j])%>% 
      select(functional_group,year,iter,sum_x)
    
    pair2_df[[i]][[j]] <- sums_regions[[i]]%>%
      filter(functional_group==pair_2[j])%>%
      select(functional_group,year,iter,sum_x)
    
    compare_df[[i]][[j]] <- merge(pair1_df[[i]][[j]],pair2_df[[i]][[j]],by=c("year","iter")) #merge comparison dfs by year and iteration
    compare_df[[i]][[j]] <- compare_df[[i]][[j]]%>%
      rename(
        group_1=functional_group.x,
        group_2=functional_group.y,
        group_1_trend=sum_x.x,
        group_2_trend=sum_x.y,
      )
  }   
}

#create lists to populate with correlations
cor_df <- vector(mode="list",length=length(regions))
cor_summary <- vector(mode="list",length=length(regions))

for (i in 1:length(regions)){ #for each region
  for(j in 1:length(pair_1)){ #for each comparison
    
      cor_df[[i]][[j]] <- compare_df[[i]][[j]]%>%
        group_by(iter)%>% #for each iteration
        summarize(
          cor=as.numeric(cor.test(group_1_trend,group_2_trend)$estimate))%>% #derive correlation between each functional group abundance trend
        ungroup()
      
      
      cor_summary[[i]][[j]] <- cor_df[[i]][[j]]%>%
        summarize(
          med_cor=median(cor)) #derive median correlation
      
      
      cor_df[[i]][[j]]<- cor_df[[i]][[j]]%>%
        mutate(
          comparison=paste(unique(compare_df[[i]][[j]]$group_1),unique(compare_df[[i]][[j]]$group_2),sep="-"), #create comparion var
          sign=ifelse(cor_summary[[i]][[j]]$med_cor>=0,"positive","negative")) #derive sign of corr
      
  }
  
  cor_df[[i]] <- rbindlist(cor_df[[i]]) #combine lists 
}


#create lists to populate with proportions of pos and negative corr
rollup_df <- vector(mode="list",length=length(regions))
rollup_fig <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){
  #calculate number (proportion) of pairwise comparisons in each comparison type that are positive vs negative
  rollup_df[[i]] <- cor_df[[i]]%>%
    dplyr::group_by(sign)%>% 
    dplyr::summarize(
      count=n())%>%
    ungroup()%>%
    mutate(
      sum_count=sum(count),
      frac=count/sum_count)
  
  #visualize a rollup barplot that depicts proportion of pos and neg pairwise comparison in each comparison type
  rollup_fig[[i]] <- ggplot(rollup_df[[i]])+
    geom_bar(mapping = aes(x = 1, y = frac, fill = sign), 
             stat = "identity")+
    theme_void()+
    scale_fill_manual(values = c("negative"="#6699CC","positive"="#DDCC77"))+
    theme(panel.spacing.x = unit(1, "mm"),
          legend.position = "none",
          strip.text.x = element_blank())+
    coord_flip()
}

#create list to populate with figures
cor_fig <- vector(mode="list",length=length(regions))
title <- "" #define title 
text_format=element_text(face="bold",size=25,color="black") #define general text format

#create legend for corr figure
corr_legend <- get_legend(ggplot(data=cor_df[[2]])+
  geom_violin(data=cor_df[[2]],aes(x=comparison,y=cor,fill=factor(sign,levels=c("positive","negative"))))+
  geom_boxplot(aes(x=comparison,y=cor),width=0.1, color="black", alpha=0.2) +
  geom_hline(aes(yintercept = 0),color="black",size=1)+
  theme_classic()+
  labs(x="Functional Group Pair",y="Correlation")+
  ylim(c(-1,1.2))+
  scale_fill_manual(values = c("positive"="#DDCC77","negative"="#6699CC"),labels=c("Positive Correlation","Negative Correlation"))+
  theme(strip.text = text_format,
        axis.title.x = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(face="bold",size=25,color="black",vjust=0.2),
        axis.text.y = text_format,
        plot.margin = unit(c(1,1,1,1),'cm'),
        panel.border = element_rect(colour = "black",fill=NA,size=2),
        strip.background = element_rect(fill = NA, colour = "black", size = 2),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text =  element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 0, 0, 0), "mm"))))

for(i in 1:length(regions)){ #for each region
  #visualize correlations for each region separately
  cor_fig[[i]] <- ggplot(data=cor_df[[i]])+
    geom_violin(data=cor_df[[i]],aes(x=comparison,y=cor,fill=sign))+
    geom_boxplot(aes(x=comparison,y=cor),width=0.1, color="black", alpha=0.2) +
    geom_hline(aes(yintercept = 0),color="black",size=1)+
    theme_classic()+
    labs(x="Functional Group Pair",y="Correlation")+
    ylim(c(-1,1.2))+
    scale_fill_manual(values = c("negative"="#6699CC","positive"="#DDCC77"))+
    theme(strip.text = text_format,
          axis.title.x = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(5, 0, 0, 0), "mm")),
          axis.title.y = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
          axis.text.x = element_text(face="bold",size=25,color="black",vjust=0.2),
          axis.text.y = text_format,
          legend.position = "none",
          plot.margin = unit(c(1,1,1,1),'cm'),
          panel.border = element_rect(colour = "black",fill=NA,size=2),
          strip.background = element_rect(fill = NA, colour = "black", size = 2))
}

fxnlgrp_cols <- c("Scraper"="#882255","Excavator"="#332288","Browser"="#999933") #define functional group colors 

#create legend for trends figure
trend_legend <- get_legend(ggplot(data=sums_regions2)+
                             geom_point(aes(x=year,y=med_x,color=functional_group))+
                             geom_line(aes(x=year,y=med_x,color=functional_group))+
                             geom_line(aes(x=year,y=med_y,color=functional_group),linetype="dashed")+
                             facet_wrap(~region,ncol=1,scales="free")+
                             geom_ribbon(aes(x=year,ymin=cpue_lower,ymax=cpue_upper,fill=functional_group),alpha=0.5)+
                             theme_classic()+
                             labs(y="Median parrotfish per survey")+
                             scale_x_continuous(breaks=seq(1993, 2020, 2))+
                             scale_color_manual(values=fxnlgrp_cols)+
                             scale_y_continuous(labels=comma)+
                             scale_fill_manual(values=fxnlgrp_cols)+
                             theme(axis.title.x = element_blank(),
                                   axis.title.y = element_text(face="bold",size=20,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
                                   axis.text.x = element_text(face="bold",size=15,color="black",angle=90,hjust=0.95,vjust=0.2),
                                   axis.text.y = text_format,
                                   strip.text = element_text(face="bold",size=15,color="black"),
                                   legend.title = element_blank(),
                                   legend.text =  element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
                                   legend.position = "bottom",
                                   legend.background = element_blank(),
                                   plot.margin = unit(c(1,1,1,1),'cm')))

trends_figs <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){ #for each regions
  #visualize abundance trends of pf aggregated by functional group
  trends_figs[[i]] <- ggplot(data=sums_regions[[i]])+
    geom_point(aes(x=year,y=med_x,color=functional_group))+
    geom_line(aes(x=year,y=med_x,color=functional_group))+
    geom_line(aes(x=year,y=med_y,color=functional_group),linetype="dashed")+
    geom_ribbon(aes(x=year,ymin=cpue_lower,ymax=cpue_upper,fill=functional_group),alpha=0.5)+
    theme_classic()+
    labs(y="Median parrotfish per survey")+
    scale_x_continuous(breaks=seq(1993, 2020, 2))+
    scale_color_manual(values=fxnlgrp_cols)+
    scale_y_continuous(labels=comma)+
    scale_fill_manual(values=fxnlgrp_cols)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
          axis.text.x = element_text(face="bold",size=25,color="black",angle=90,hjust=0.95,vjust=0.2),
          axis.text.y = text_format,
          strip.text = element_text(face="bold",size=15,color="black"),
          legend.position = "none",
          plot.margin = unit(c(0.7,0.7,0.7,0.7),'cm'))
}


#organize corr figures by row
cor_row_1 <- plot_grid(title,rollup_fig[[1]],cor_fig[[1]]+theme(axis.text.x = element_blank(),axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15))
cor_row_2 <- plot_grid(title,rollup_fig[[2]],cor_fig[[2]]+theme(axis.text.x = element_blank(),
                                                                axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15))
cor_row_3 <- plot_grid(title,rollup_fig[[3]],cor_fig[[3]]+theme(axis.text.x = element_blank(),
                                                                axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15))
cor_row_4 <- plot_grid(title,rollup_fig[[4]],cor_fig[[4]]+theme(axis.text.x = element_blank(),
                                                                axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15))
cor_row_5 <- plot_grid(title,rollup_fig[[5]],cor_fig[[5]]+theme(axis.text.x = element_blank(),
                                                                axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15))
cor_row_6 <- plot_grid(title,rollup_fig[[6]],cor_fig[[6]], axis = "lrtb",align = "v",ncol = 1, rel_heights = c(0.4,0.5, 15))

#organize last trend figure with legend (for aesthetic purposes)
trend_6 <- plot_grid(trends_figs[[6]],trend_legend,ncol=1,rel_heights = c(1,0.2))

row_1 <- plot_grid(trends_figs[[1]]+theme(axis.text.x = element_blank(),axis.title.x = element_blank()), cor_row_1,nrow=1,ncol=2,labels=c(paste(paste(LETTERS[1],".",sep=""),regions[1],sep=" "),paste(LETTERS[2],".",sep="")),label_size = 25,hjust = 0, label_x = 0.01,label_y=1,vjust=1,align = "h", axis = "b",rel_heights=c(1.5,1))
row_2 <- plot_grid(trends_figs[[2]]+theme(axis.text.x = element_blank(),axis.title.x = element_blank()),cor_row_2,nrow=1,ncol=2,labels=c(paste(paste(LETTERS[3],".",sep=""),regions[2],sep=" "),paste(LETTERS[4],".",sep="")),label_size = 25,hjust = 0, label_x = 0.01,label_y=1,vjust=1,align = "h", axis = "b",rel_heights=c(1.5,1))
row_3 <- plot_grid(trends_figs[[3]]+theme(axis.text.x = element_blank(),axis.title.x = element_blank()),cor_row_3,nrow=1,ncol=2,labels=c(paste(paste(LETTERS[5],".",sep=""),regions[3],sep=" "),paste(LETTERS[6],".",sep="")),label_size = 25,hjust = 0, label_x = 0.01,label_y=1,vjust=1,align = "h", axis = "b",rel_heights=c(1.5,1))
row_4 <- plot_grid(trends_figs[[4]]+theme(axis.text.x = element_blank(),axis.title.x = element_blank()),cor_row_4,nrow=1,ncol=2,labels=c(paste(paste(LETTERS[7],".",sep=""),regions[4],sep=" "),paste(LETTERS[8],".",sep="")),label_size = 25,hjust = 0, label_x = 0.01,label_y=1,vjust=1,align = "h", axis = "b",rel_heights=c(1.5,1))
row_5 <- plot_grid(trends_figs[[5]]+theme(axis.text.x = element_blank(),axis.title.x = element_blank()),cor_row_5,nrow=1,ncol=2,labels=c(paste(paste(LETTERS[9],".",sep=""),regions[5],sep=" "),paste(LETTERS[10],".",sep="")),label_size = 25,hjust = 0, label_x = 0.01,label_y=1,vjust=1,align = "h", axis = "b",rel_heights=c(1.5,1))
row_6 <- plot_grid(trends_figs[[5]],cor_row_6,nrow=1,ncol=2,labels=c(paste(paste(LETTERS[11],".",sep=""),regions[6],sep=" "),paste(LETTERS[12],".",sep="")),label_size = 25,hjust = 0, label_x = 0.01,label_y=1,vjust=1,align = "h", axis = "b",rel_heights=c(0.1,1))
row_7 <- plot_grid(trend_legend,corr_legend,ncol=2)

pdf('./figures/regional_functionalgrp_trends_corr.pdf',height=50,width=30)
plot_grid(row_1,row_2,row_3,row_4,row_5,row_6,row_7,ncol=1,rel_heights = c(rep(1,times=6),0.3))          
dev.off()

    
    