# Run state-space model on 10 pf species by region
#Jordan DiNardo

rm(list=ls())

#load necessary libraries
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);
library(stringr);library(lubridate);library(rlist)library(here); library(sp); 
library(shinystan);library(data.table)
here()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load necessary functions
source("./functions_models/reef_functions.R") #data prep
source("./functions_models/regional_ss_model.R") #multi_sp ss model

#load prepared data
load("./data_prep/regional_species_surveys_final.RData") # read in pf occurrence data developed from data_cleaning_casestudy_pf_WCR.R file

#create a species var and ind species vars
for(i in 1:length(regions)){
  pf_occs[[i]]$species<- NA
  pf_occs[[i]]$species.rb<- NA
  pf_occs[[i]]$species.rt<- NA
  pf_occs[[i]]$species.y<- NA
  pf_occs[[i]]$species.sto<- NA
  pf_occs[[i]]$species.ra<- NA
  pf_occs[[i]]$species.m<- NA
  pf_occs[[i]]$species.q<- NA
  pf_occs[[i]]$species.p<-NA
  pf_occs[[i]]$species.b<- NA
  pf_occs[[i]]$species.str<- NA
  
  pf_occs[[i]]<- pf_occs[[i]]%>% 
    mutate(
      species=ifelse(speciesid==pf_speciesid[1],pf_commonnames[1],species),
      species=ifelse(speciesid==pf_speciesid[2],pf_commonnames[2],species),
      species=ifelse(speciesid==pf_speciesid[3],pf_commonnames[3],species),
      species=ifelse(speciesid==pf_speciesid[4],pf_commonnames[4],species),
      species=ifelse(speciesid==pf_speciesid[5],pf_commonnames[5],species),
      species=ifelse(speciesid==pf_speciesid[6],pf_commonnames[6],species),
      species=ifelse(speciesid==pf_speciesid[7],pf_commonnames[7],species),
      species=ifelse(speciesid==pf_speciesid[8],pf_commonnames[8],species),
      species=ifelse(speciesid==pf_speciesid[9],pf_commonnames[9],species),
      species=ifelse(speciesid==pf_speciesid[10],pf_commonnames[10],species),
      species.rb=ifelse(speciesid==pf_speciesid[1],1,0), 
      species.rt=ifelse(speciesid==pf_speciesid[2],1,0),
      species.y=ifelse(speciesid==pf_speciesid[3],1,0),
      species.sto=ifelse(speciesid==pf_speciesid[4],1,0),
      species.ra=ifelse(speciesid==pf_speciesid[5],1,0),
      species.m=ifelse(speciesid==pf_speciesid[6],1,0),
      species.q=ifelse(speciesid==pf_speciesid[7],1,0),
      species.p=ifelse(speciesid==pf_speciesid[8],1,0),
      species.b=ifelse(speciesid==pf_speciesid[9],1,0),
      species.str=ifelse(speciesid==pf_speciesid[10],1,0)
      
    )
  
  pf_occs[[i]]$R_year<- paste(pf_occs[[i]]$species,pf_occs[[i]]$year,sep='_') #create a species_year var
  pf_occs[[i]]<- subset(pf_occs[[i]],is.na(species)==F)}

N_R <- length(pf_speciesid) #define number of species

#create list to populate
pf_ts <- vector(mode="list",length=length(regions))

# use fxn to create timeseries estimates
for (i in 1:length(regions)){
  pf_ts[[i]] <- ts_reef_sp(pf_occs[[i]]) # use ts_reef_sp fxn to create timeseries df
}

X <- vector(mode="list",length=length(regions))

pf_SS_mr<- vector(mode="list",length=length(regions)) #create list for models for each region
pf_param<- vector(mode="list",length=length(regions)) #create list for model outputs for each region

#filter data to surveys conducted before covid (before 2020; keep data up until 2019)
for(i in 1:length(regions)){
  pf_occs[[i]]<- filter(pf_occs[[i]],year<2020)
}

#run multi-species ss model for each region 
for (i in 1:length(regions)){
  
  X[[i]]<- matrix(data=c(scale(as.numeric(pf_occs[[i]]$btime)),scale(as.numeric(pf_occs[[i]]$averagedepth)),scale(as.numeric(pf_occs[[i]]$visibility)),scale(as.numeric(pf_occs[[i]]$current)),pf_occs[[i]]$exp_binary),ncol=5,nrow=nrow(pf_occs[[i]]))
  
  pf_SS_mr[[i]]<- rstan::stan(model_code = SS_trend_ord_multi, data = list(y=pf_occs[[i]]$abundance2,
                                                                           N = nrow(pf_occs[[i]]),
                                                                           site=as.numeric(factor(pf_occs[[i]]$geogr)),
                                                                           N_site=length(unique(pf_occs[[i]]$geogr)),
                                                                           diver=as.numeric(factor(pf_occs[[i]]$memberid)),
                                                                           N_dv=length(unique(pf_occs[[i]]$memberid)),
                                                                           dmy=as.numeric(factor(pf_occs[[i]]$site_dmy)),
                                                                           N_dmy=length(unique(pf_occs[[i]]$site_dmy)),
                                                                           K=length(unique(pf_occs[[i]]$abundance)),
                                                                           X=X[[i]],
                                                                           Z=ncol(X[[i]]),
                                                                           TT=max(pf_occs[[i]]$year)-min(pf_occs[[i]]$year)+1,
                                                                           N_yr=length(unique(pf_occs[[i]]$year)),
                                                                           yr_index=sort(unique(as.numeric(factor(pf_occs[[i]]$year)))),
                                                                           year_id=as.numeric(factor(pf_occs[[i]]$year)),
                                                                           N_R=length(unique(pf_occs[[i]]$species)),
                                                                           R=as.numeric(factor(pf_occs[[i]]$species)),
                                                                           NC_1=45,#set number of correlations
                                                                           year_id_R=as.numeric(factor(pf_occs[[i]]$R_year))),
                              pars = c('c','sd_site','sd_dmy','sd_dv','sd_r','sd_q','x','a_yr','beta','cor_1','q_mat','z_r','z_q'),
                              control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 1000, chains = 4, iter = 2000, refresh=250,thin = 1,init=0)
}

#extract params for each region
for (i in 1:length(regions)){
  pf_param[[i]]<- rstan::extract(pf_SS_mr[[i]])
}

save(pf_occs,pf_param,pf_SS_mr,X,pf_ts,file="./model_output/regional_species_results_final.RData") #save outputs 

