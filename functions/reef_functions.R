reef_filter_sp = function(R,GZ,sp){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()

  
  TempDat<-R %>% subset(geogr4 %in% GZ) %>% select('formid','speciesid','abundance',everything())
  
  TempDat_sp<- TempDat %>% subset(speciesid%in%sp$speciesid)
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2)
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- complete(TempDat2,formid,nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
    site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid)) %>% subset(n>=5) #Calculate surveys per site
  
  TempDat4<- TempDat3 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat4 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat4<- TempDat4 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  
  occ_dat<- subset(TempDat4,speciesid%in%sp$speciesid) #Subset out each species in the provided dataframe
  return(occ_dat)
}


#estimate timeseries of occ, abund, n.survs
ts_reef = function(X){
  occ_by_year<- X %>% group_by(year) %>% dplyr::summarize(n.occ=sum(occ),n.surv=n_distinct(formid),p.occ=n.occ/n.surv)
  abun_by_year<- X %>% group_by(year,geogr) %>% dplyr::summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n_distinct(formid)) %>% group_by(year) %>% dplyr::summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
  total_sd<- X %>% group_by(year) %>% dplyr::summarize(sd=sd(abund_trans))
  
  comb<- left_join(occ_by_year,abun_by_year)
  comb2<- left_join(comb,total_sd)
  ts_dat<- comb2
  return(ts_dat)
}

#estimate timeseries of occ, abund, n.survs by species
ts_reef_sp = function(X){
  occ_by_year<- X %>% 
    group_by(year,speciesid) %>% 
    dplyr::summarize(n.occ=sum(occ),
                     n.surv=n_distinct(formid),
                     p.occ=n.occ/n.surv)
  abun_by_year<- X %>% 
    group_by(year,geogr,speciesid) %>% 
    dplyr::summarize(site_abun=mean(abund_trans),
                     sd_abund=sd(abund_trans),
                     n.surv=n_distinct(formid)) %>% 
    ungroup() %>%
    group_by(year,speciesid) %>% 
    dplyr::summarize(sp=unique(speciesid),
                     mean_abund=mean(site_abun),
                     n.survs=sum(n.surv),
                     n.sites=n())
  total_sd<- X %>% 
    group_by(year,speciesid) %>% 
    dplyr::summarize(sd=sd(abund_trans))
  
  comb<- left_join(occ_by_year,abun_by_year)
  comb2<- left_join(comb,total_sd)
  ts_dat<- comb2
  return(ts_dat)
}


#transform abund based on prob of landing in each category
abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

#convert ordered categories to number
ord_to_n<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
  if(ncol(c)==2){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=1-plogis(c[,2]-x)
    p[,4]=0
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==3){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=1-plogis(c[,3]-x)
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==4){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=plogis(c[,4]-x[,i])-plogis(c[,3]-x[,i])
    p[,5]=1-plogis(c[,4]-x[,i])
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  return(abund_x)
}

estimate_trends <- function(sp,region,grp,params1,TT,ts,n.groups){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params1$c)))
      
      if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==3){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
  } 
  return(lambda_mat)
  
  for(r in 1:n.groups){
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x)
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.95)
    }
    x_mat[[r]]$species <- sp[r]
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y)
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.95)
    }
  }
  return(lambda_mat)
}



