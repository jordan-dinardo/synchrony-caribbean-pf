# Clean & prep data for regional species ss modeling
# Jordan DiNardo

rm(list=ls())

#Load libraries
library(rstan);library(loo);library(bayesplot);library(dplyr);library(tidyverse);
library(stringr);library(lubridate);library(rlist);library(here); library(sp)
here()

#Load functions
source("./functions_models/reef_functions.R") #data prep

## data cleaning ## 

#taxonomic data
pf_commonnames <- c("Redband","Redtail","Yellowtail (Redfin)","Stoplight","Rainbow","Midnight","Queen","Princess","Blue","Striped") #vector of pf common names
pf_speciesid <- c(152,154,153,155,151,148,150,149,144,156) #vector of pf species ids
regions <- c("Bonaire", "Key Largo", "Cozumel", "Little Cayman", "Roatan","Tortolla") #vector of regions
gr_codes<- c(8503,3403,5402,5202,5702,6501) #vector of region codes
#---------------------------------------
fish_reef<- read.csv(here('data_raw2','parrotfish_species.csv')) #read in fish species found in the Tropical Western Atlantic
parrotfish<- subset(fish_reef,speciesid %in% pf_speciesid) #subset pf species for study
pf_scientificnames <- unique(parrotfish$scientificname) #vector of pf sci names


#survey data
REEF<- read.csv(here("data_raw2", "parrotfish_021023.csv")) #read in REEF survey sighting data 
REEF_survey<- read.table(here("data_raw2","twa_surveys_021023.txt"),header=T) #read in survey metadata
m<- match(REEF$Form,REEF_survey$formid) #match sightings to surveys
REEF[,17:31]<- REEF_survey[m,2:16] #merge survey level data

#geog data
reef_geog<- read.csv(here("data_raw2","REEF_TWAgeog.csv")) #read in REEF region metadata
site_surveys<-REEF %>% group_by(geogr) %>% summarize(n=n_distinct(Form))
m<- match(reef_geog$geogid,site_surveys$geogr) #find matches between reef geog data and sites seen in reef survey sighting data
reef_geog$no.surveys<- site_surveys$n[m] #calc no of surveys
reef_geog$region.id<- substr(reef_geog$geogid, 1,4) #get the region id (first four digits)
reef_geog$lat_full<- reef_geog$lat #copy of the full latitude
reef_geog$lon_full<- reef_geog$lon #copy of the full longitude
reef_geog<-reef_geog%>% #separate out degrees and minutes
  separate(lat,into=c("lat_deg","lat_min"),sep=" ")%>%
  separate(lon,into=c("lon_deg","lon_min"),sep=" ")
col.num<-c("lat_deg","lat_min","lon_deg","lon_min") #create vec of vars to convert to numeric
reef_geog[col.num]<-sapply(reef_geog[col.num],as.numeric) #convert vars to numeric
reef_geog<-reef_geog[complete.cases(reef_geog),] #discard sites with no coordinates
reef_geog$lat_dd<- reef_geog$lat_deg+reef_geog$lat_min/60 #convert lat to decimal degrees
reef_geog$lon_dd<- reef_geog$lon_deg-reef_geog$lon_min/60 #convert long to decimal degrees
REEF<-REEF[-which(REEF$Date=="0000-00-00"),] #remove rows with no date reported
REEF$Date <- strptime(REEF$Date, format = "%m/%d/%Y") #convert date into format that is recognized by lubridate
REEF$Date<-ymd(REEF$Date)# change date variable to proper date format
REEF<-cbind(REEF,year=year(REEF$Date)) #create year var from date field
REEF<-cbind(REEF,month=month(REEF$Date)) #create month var from date field
REEF<-cbind(REEF,day=day(REEF$Date)) #create day var from date field
REEF<-filter(REEF, REEF$geogr %in% reef_geog$geogid) #remove sites without lat/long 
REEF<-REEF[as.numeric(REEF$btime)>20,] #exclude dives shorter than 20 mins
REEF<-REEF[as.numeric(REEF$btime)<120,] #exclude dives longer than 120 mins 
REEF<-REEF[as.numeric(REEF$start)>5,] #exclude dives before 5 am 
REEF<-REEF[as.numeric(REEF$start)<20,] #exclude dives after 8 pm 
REEF <- rename(REEF,c("formid"="Form","speciesid"="Species","abundance"="Abundance"))

## data filtering ## 

pf_occs <- vector(mode="list",length=length(regions))

for (i in 1:length(regions)){
  
  pf_occs[[i]] <- reef_filter_sp(REEF,sp=parrotfish,GZ=gr_codes[i])
}

save(pf_scientificnames,pf_commonnames,pf_speciesid, pf_occs,regions,gr_codes,file="./data_prep/regional_species_surveys_final.RData")



