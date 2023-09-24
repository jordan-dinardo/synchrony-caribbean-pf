#  Visualize regional interspecies correlations
# Jordan DiNardo

rm(list=ls())

#load necessary libraries
library(data.table);library(tidyverse);library(plyr);library(ggpubr);library(ggh4x);
library(matrixStats);library(cowplot);library(stringi)

load(file="./model_output/regional_species_results_final.RData") #load region-specific model results

regions <- c("Bonaire", "Key Largo", "Cozumel", "Little Cayman", "Roatan","Tortolla") #define regions
pf_commonnames <- c("Redband","Redtail","Yellowtail (Redfin)","Stoplight","Rainbow","Midnight","Queen","Princess","Blue","Striped") #define species common names
pf_speciesid <- c(152,154,153,155,151,148,150,149,144,156) #define species id
N_R <- length(pf_speciesid) #define number of species species

#create lists to populate
corrs <- vector(mode="list",length=length(regions))
cor_med <- vector(mode="list",length=length(regions))
cor_mat <- vector(mode="list",length=length(regions))
cor_mat_df <- vector(mode="list",length=length(regions))
cor_mat_table <- vector(mode="list",length=length(regions))
domains <- vector(mode="list",length=length(regions))
pf_names_mat <- vector(mode="list",length=length(regions))
pf_pair_names <- vector(mode="list",length=length(regions))
corrs_long <- vector(mode="list",length=length(regions))
corr_df <- vector(mode="list",length=length(regions))

##Extract correlations from model output and populate matrices
for(i in 1:length(regions)){ #for each region
  corrs[[i]] <- round(pf_param[[i]]$cor_1,4) #extract correlation from model estimates (round to 4 digits)
  cor_med[[i]] <- round(colMedians(pf_param[[i]]$cor_1),4) #calculate the median correlation for a given pairwise species (round to 4 digits)
  cor_mat[[i]] <- matrix(NA,N_R,N_R) #create a matrix N_RxN_R species
  colnames(cor_mat[[i]])<-sort(pf_commonnames) #define column names as (common) species names
  rownames(cor_mat[[i]])<-sort(pf_commonnames) #define row names as (common) species names
  
  pf_names_mat[[i]] <- matrix(NA,N_R,N_R) #create a matrix N_RxN_R species
  colnames(pf_names_mat[[i]])<-sort(pf_commonnames) #define column names as (common) species names
  rownames(pf_names_mat[[i]])<-sort(pf_commonnames) #define row names as (common) species names
  
  for (k in 1:N_R){ # for each species
    for (j in 1:(k - 1)) { 
      cor_mat[[i]][j, k] = cor_med[[i]][choose(k - 1, 2) + j] #populate matrix with median corr
    }
  }
  
  cor_mat[[i]][1,1] <- NA #add NA to 1,1 position
  cor_mat[[i]][is.na(cor_mat[[i]])] <- t(cor_mat[[i]])[is.na(cor_mat[[i]])] #transpose matrix
  
  
  cor_mat[[i]][lower.tri(cor_mat[[i]])] <- NA #populate lower triangular part of matrix with NAs
  
  cor_mat_df[[i]] <- data.frame(cor_mat[[i]]) #convert matrix to dataframe
  colnames(cor_mat_df[[i]])<- sort(pf_commonnames) #set col names to species names
  
}

scrapers <- c("Princess","Striped","Queen","Blue") #define species to scrapers 
excavators <- c("Stoplight","Rainbow","Midnight") #define species to excavators 
browsers <- c("Redband", "Redtail","Redfin") #define species to browsers 


for(i in 1:length(regions)){ 
  for (k in 1:N_R){ # rows
    for (j in 1:N_R) { #cols
      pf_names_mat[[i]][j,k] <- paste(colnames(pf_names_mat[[i]])[j],colnames(pf_names_mat[[i]])[k],sep="-") #define names of interspecies  comparison (species 1- species 2)
    }
  }
  
  pf_names_mat[[i]][lower.tri(pf_names_mat[[i]], diag = T)] <- NA #transpose matrix
  pf_pair_names[[i]] <-pf_names_mat[[i]][upper.tri(pf_names_mat[[i]],diag=F)] #extract upper triangular part
  corrs[[i]] <- data.frame(corrs[[i]]) #convert matrix to df
  colnames(corrs[[i]]) <- pf_pair_names[[i]] #define col names to interspecies comparison ids
  corrs_long[[i]] <- pivot_longer(corrs[[i]],cols=1:ncol(corrs[[i]]),names_to = "pair",values_to="corr") #convert (wide) df to long df
  corrs_long[[i]]$pair <- stringi::stri_replace_all_fixed(corrs_long[[i]]$pair,'Yellowtail (Redfin)','Redfin') #rename sp to shorten name (for aesthetics)
}


for(i in 1:length(regions)){ #for each region
  
  corr_df[[i]] <- data.frame(species_1=rownames(cor_mat[[i]])[row(cor_mat[[i]])], species_2=colnames(cor_mat[[i]])[col(cor_mat[[i]])], med_corr=c(cor_mat[[i]]))  #define vars in df
  corr_df[[i]] <- na.omit(corr_df[[i]]) #omit NAs
  corr_df[[i]]$region = regions[i] #define region
  corr_df[[i]]$species_1 <- ifelse(corr_df[[i]]$species_1=="Yellowtail (Redfin)","Redfin",corr_df[[i]]$species_1) #change species common name to Redfin
  corr_df[[i]]$species_2 <- ifelse(corr_df[[i]]$species_2=="Yellowtail (Redfin)","Redfin",corr_df[[i]]$species_2) #change species common name to Redfin
  corr_df[[i]]$pair1 <- paste(corr_df[[i]]$species_1,corr_df[[i]]$species_2,sep="-") #define pair
  corr_df[[i]]$pair2 <- paste(corr_df[[i]]$species_2,corr_df[[i]]$species_1,sep="-") #define pair (opp order of sp) in case original ordder is not recognized
  corr_df[[i]]$sign <- ifelse(corr_df[[i]]$med_corr>=0,"positive","negative") #define sign of corr
  
  corr_df[[i]]$fnxl_grp <- NA #create functional group var
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_1%in% scrapers & corr_df[[i]]$species_2%in% scrapers,"Scraper",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_1%in% browsers & corr_df[[i]]$species_2%in% browsers,"Browser",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_1%in% excavators & corr_df[[i]]$species_2%in% excavators,"Excavator",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group
  
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_2%in% scrapers & corr_df[[i]]$species_1%in% excavators,"Scraper-Excavator",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group comparison
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_2%in% scrapers & corr_df[[i]]$species_1%in% browsers,"Scraper-Browser",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group comparison
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_2%in% excavators & corr_df[[i]]$species_1%in% browsers,"Excavator-Browser",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group comparison
  #repeat with opposite order of species 
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_1%in% scrapers & corr_df[[i]]$species_2%in% excavators,"Scraper-Excavator",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group comparison
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_1%in% scrapers & corr_df[[i]]$species_2%in% browsers,"Scraper-Browser",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group comparison
  corr_df[[i]]$fnxl_grp <- ifelse(corr_df[[i]]$species_1%in% excavators & corr_df[[i]]$species_2%in% browsers,"Excavator-Browser",corr_df[[i]]$fnxl_grp) #assign pairwise comparison to corresponding functional group comparison
  
  corr_df[[i]]$comparison_type <- NA #create comparison type var
  corr_df[[i]]$comparison_type <- ifelse(corr_df[[i]]$species_1%in% scrapers & corr_df[[i]]$species_2%in% scrapers,"Within Functional Groups",corr_df[[i]]$comparison_type) #assign pairwise comparison to corresponding comparison type
  corr_df[[i]]$comparison_type <- ifelse(corr_df[[i]]$species_1%in% browsers & corr_df[[i]]$species_2%in% browsers,"Within Functional Groups",corr_df[[i]]$comparison_type) #assign pairwise comparison to corresponding comparison type
  corr_df[[i]]$comparison_type <- ifelse(corr_df[[i]]$species_1%in% excavators & corr_df[[i]]$species_2%in% excavators,"Within Functional Groups",corr_df[[i]]$comparison_type) #assign pairwise comparison to corresponding comparison type
  corr_df[[i]]$comparison_type <- ifelse(is.na(corr_df[[i]]$comparison_type),"Across Functional Groups",corr_df[[i]]$comparison_type) #assign pairwise comparison to corresponding comparison type
}

#create lists to populate
corrs_long_merge1 <- vector(mode="list",length=length(regions)) 
corrs_long_merge2 <- vector(mode="list",length=length(regions))
corrs_long_merge <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){ #for each region
  #merge corrs_long and corr_df (done twice in case order of species is not recognized)
  corrs_long_merge1[[i]] <- merge(corrs_long[[i]],corr_df[[i]],by.x="pair",by.y="pair1") 
  corrs_long_merge2[[i]] <- merge(corrs_long[[i]],corr_df[[i]],by.x="pair",by.y="pair2")
  corrs_long_merge1[[i]] <- subset(corrs_long_merge1[[i]],select=-c(pair2)) #remove duplicates
  corrs_long_merge2[[i]] <- subset(corrs_long_merge2[[i]],select=-c(pair1)) #remove duplicates
  corrs_long_merge[[i]] <- rbind(corrs_long_merge1[[i]],corrs_long_merge2[[i]]) #combine data
}

#calculate summary stats by functional group (within functional groups) 
  #and functional group comparisons (across functional groups)
corr_summary <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){ #for each region
  
  corr_summary[[i]] <- corrs_long_merge[[i]]%>%
    dplyr::group_by(fnxl_grp)%>% 
    dplyr::summarise(
      group_corr = median(corr) #calc median corr for a given functional group and across functional group
    )%>%
    mutate(
      corr_sign = ifelse(group_corr>=0,"positive","negative") #define sign of median corr
    )
  
  corrs_long_merge[[i]] <- merge(corrs_long_merge[[i]],corr_summary[[i]],by="fnxl_grp") #merge functional group summary stats with corrs_long_merge df
}


#calculate summary stats for comparison type (within and across functional groups)
rollup_df <- vector(mode="list",length=length(regions))
rollup_fig <- vector(mode="list",length=length(regions))

for(i in 1:length(regions)){ #for each region
  #calculate number (proportion) of interspecies comparisons in each comparison type that are positive vs negative
  rollup_df[[i]] <- corrs_long_merge[[i]]%>%
    dplyr::group_by(comparison_type,sign)%>%
    dplyr::summarize(
      count=n())%>%
    ungroup()%>%
    group_by(comparison_type)%>%
    mutate(
      sum_count=sum(count),
      frac=count/sum_count)
  
  #visualize a rollup barplot that depicts proportion of pos and neg pairwise comparison in each comparison type
  rollup_fig[[i]] <- ggplot(rollup_df[[i]])+
    geom_bar(mapping = aes(x = 1, y = frac, fill = sign), 
             stat = "identity")+
    theme_void()+
    scale_fill_manual(values = c("positive"="#DDCC77","negative"="#6699CC"),labels=c("Positive","Negative"))+
    theme(panel.spacing.x = unit(1, "mm"),
          legend.position = "none",
          strip.text.x = element_blank())+
    facet_grid(.~factor(comparison_type,levels=c("Within Functional Groups","Across Functional Groups")), scales = "free_x")+
    coord_flip()
}

#create a legend for figure
corr_legend <- get_legend(ggplot(rollup_df[[1]])+
  geom_bar(mapping = aes(x = 1, y = frac, fill = sign), 
           stat = "identity")+
  theme_void()+
  scale_fill_manual(values = c("positive"="#6699CC","negative"="#DDCC77"),labels=c("Positive","Negative"))+
  theme(panel.spacing.x = unit(1, "mm"),
        legend.position = "bottom",
        legend.title=element_blank(),
        legend.text=element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
        strip.text.x = element_blank())+
  facet_grid(.~factor(comparison_type,levels=c("Within Functional Groups","Across Functional Groups")), scales = "free_x")+
  coord_flip())

#create lists for figures
corr_fig <- vector(mode="list",length=length(regions))
fnxl_fig <- vector(mode="list",length=length(regions))
title <- ""

#define general text format
text_format=element_text(face="bold",size=25,color="black") #set general text format


#Visualize correlations for each region separately
for(i in 1:length(regions)){ #for each region
  
  corr_fig[[i]] <- ggplot(data=corrs_long_merge[[i]])+
    geom_violin(data=corrs_long_merge[[i]],aes(x=pair,y=corr,fill=sign))+
    geom_boxplot(aes(x=pair,y=corr),width=0.1, color="black", alpha=0.2) +
    geom_hline(aes(yintercept = 0),color="black",size=1)+
    geom_hline(aes(yintercept = group_corr,color=corr_sign),size=1)+
    facet_nested(~factor(comparison_type,levels=c("Within Functional Groups","Across Functional Groups"))+fnxl_grp,scales="free_x")+
    theme_classic()+
    labs(x="Species Pair",y="Species Correlation")+
    ylim(c(-1,1.2))+
    scale_fill_manual(values = c("positive"="#DDCC77","negative"="#6699CC"),labels=c("Positive","Negative"))+
    scale_color_manual(values = c("positive"="#DDCC77","negative"="#6699CC"))+
    theme(strip.text = text_format,
          axis.title.x = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(5, 0, 0, 0), "mm")),
          axis.title.y = element_text(face="bold",size=25,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
          axis.text.x = element_text(face="bold",size=20,color="black",angle=90,hjust=0.95,vjust=0.2),
          axis.text.y = text_format,
          legend.position = "none",
          plot.margin = unit(c(1,1,1,1),'cm'),
          panel.border = element_rect(colour = "black",fill=NA,size=2),
          strip.background = element_rect(fill = NA, colour = "black", size = 2))
  fnxl_fig[[i]] <- plot_grid(title,rollup_fig[[i]], corr_fig[[i]], align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15))
  pdf(paste("./figures/regional_interspecies_corr_",regions[i],".pdf",sep=""),height=10,width=25)
  print(fnxl_fig[[i]])
  dev.off()
  
}


#arrange plots for each region to visualize in a single figure
pdf("./figures/regional_interspecies_corr.pdf",height=70,width = 40)
plot_grid(
  plot_grid(title,rollup_fig[[1]],corr_fig[[1]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15)),
  plot_grid(title,rollup_fig[[2]],corr_fig[[2]]+theme(axis.text.x = element_blank(),
                                                               axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15)),
  plot_grid(title,rollup_fig[[3]],corr_fig[[3]]+theme(axis.text.x = element_blank(),
                                                                axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15)),
  plot_grid(title,rollup_fig[[4]],corr_fig[[4]]+theme(axis.text.x = element_blank(),
                                                              axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15)),
  plot_grid(title,rollup_fig[[5]],corr_fig[[5]]+theme(axis.text.x = element_blank(),
                                                              axis.title.x = element_blank()),align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15)),
  plot_grid(title,rollup_fig[[6]],corr_fig[[6]],align = "v", axis = "lrtb",ncol = 1, rel_heights = c(0.4,0.5, 15)),
  corr_legend,ncol=1,nrow = 7,labels = c(paste(paste(LETTERS[1:6],".",sep=""),regions[1:6],sep=" "),""),label_size = 25,hjust = 0,label_x = 0.01,label_y=1.01,rel_heights=c(rep(1,times=6),0.2))

dev.off()
