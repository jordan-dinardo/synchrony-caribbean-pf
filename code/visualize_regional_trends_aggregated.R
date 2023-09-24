#Visualize aggregated parrotfish relative abundance timeseries by region
#Jordan DiNardo

rm(list=ls())

#load necessary libraries
library(data.table);library(tidyverse);library(scales)

load(file="./model_output/regional_species_results_final.RData") #load region-specific model results
source("./functions/reef_functions.R") #load functions

regions <- c("Bonaire", "Key Largo", "Cozumel", "Little Cayman", "Roatan","Tortolla") #define regions
pf_commonnames <- c("Redband","Redtail","Yellowtail (Redfin)","Stoplight","Rainbow","Midnight","Queen","Princess","Blue","Striped") #define species common name
pf_speciesid <- c(152,154,153,155,151,148,150,149,144,156) #define species id
N_R <- length(pf_speciesid) #set no species

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

#derive regional aggregated relative abundance of pf 
x_trends_new2 <- rbindlist(x_trends_new)

x_trends_new2 <- x_trends_new2%>%
  filter(year>1993)%>%
  group_by(region,year,iter)%>% #for each region, year, and iteration
  summarise(
    sum_x=sum(lambda.x,na.rm=T), #sum over species abundance (process)
    sum_y=sum(lambda.y,na.rm=T))%>% #sum over species abundance (observation)
    ungroup()%>%
    group_by(region,year)%>%
    summarise(
      med_x= median(sum_x), #derive median abundance of process
      med_y=median(sum_y), #derive median abundance of observation
      cpue_lower=quantile(sum_x,0.05), #derive 0.05 quantile
      cpue_upper=quantile(sum_x,0.95))%>% #derive 0.95 quantile
    ungroup()

#define regional colors
region_cols <- c("Bonaire"="#E69F00","Key Largo"="#009E73","Cozumel"="#D55E00",
                 "Little Cayman"="#CC79A7","Roatan"="#0072B2","Tortolla"="#661100")
#define text format
text_format <- element_text(face="bold",size=20,color="black")

pdf('./figures/regional_trends.pdf',height=5,width=12)
ggplot(data=x_trends_new2)+
  geom_point(aes(x=year,y=med_x,color=region))+
  geom_line(aes(x=year,y=med_x,color=region))+
  geom_line(aes(x=year,y=med_y,color=region),linetype="dashed")+
  geom_ribbon(aes(x=year,ymin=cpue_lower,ymax=cpue_upper,fill=region),alpha=0.5)+
  theme_classic()+
  labs(y="Median parrotfish per survey")+
  scale_x_continuous(breaks=seq(1994, 2019, 1))+
  scale_color_manual(values=region_cols)+
  scale_y_continuous(labels=comma)+
  scale_fill_manual(values=region_cols)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=20,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(face="bold",size=15,color="black",angle=90,hjust=0.95,vjust=0.2),
        axis.text.y = text_format,
        legend.title = element_blank(),
        legend.text =  element_text(face="bold",size=20,color="black",vjust=0.5,margin = unit(c(0, 5, 0, 0), "mm")),
        legend.position = c(0.15,0.8),
        legend.background = element_blank(),
        plot.margin = unit(c(1,1,1,1),'cm'))
dev.off()


