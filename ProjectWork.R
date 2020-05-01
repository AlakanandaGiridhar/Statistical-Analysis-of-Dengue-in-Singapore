setwd("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/R files")

#1. CALLING ALL REQUIRED LIBRARIES

library(dplyr)
library(reshape)
library(tidyr)
library(RColorBrewer)
library(trafo)
library(car)
library(MASS)
library(ggplot2)
library(multcompView)
library(stringr)
library(Hmisc)
library(corrplot)

#2. USER DEFINED FUNCTIONS
  #MODE FUNCTION
  
  getmode=function(v)
  {
    uniqv=unique(v)
    uniqv[which.max(tabulate(match(v,uniqv)))]
  }
  
  #TukeyPlot - Labels

  generate_label_df <- function(TUKEY, variable){
    
    # Extract labels and factor levels from Tukey post-hoc 
    Tukey.levels <- TUKEY[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    #to put the labels in the same order as in the boxplot :
    Tukey.labels$Year=rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$Year) , ]
    return(Tukey.labels)
  }
  
#3. READING DATASETS AND DATA WRANGLING
  #Sample Cluster-wise data transcribed from PDF to CSV manually; Data available for roughly 6-8 months, noted at 12 different dates 

  SampleDataLoc = read.csv("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/Trends-Normal+Epidemic+WB.csv")
  
  #Wrangling SampleData to obtain dataset for testing hypothesis about Reservoirs vs Normal locality
  SampleDataLocWB=SampleDataLoc %>% group_by(Date.Recorded,Region,Locality,WB) %>% summarise(Cluster.Cases=sum(Cluster.Cases))
  SampleDataLocWB=subset(SampleDataLocWB, Region %in% c("East","West"))
  SampleDataLocWB$WaterBodyExistence=NA
      for(i in 1:nrow(SampleDataLocWB))
      {
        if(SampleDataLocWB$WB[i]==1|SampleDataLocWB$WB[i]==2)
        {
          SampleDataLocWB$WaterBodyExistence[i]=1
        } else
        {
          SampleDataLocWB$WaterBodyExistence[i]=0
        }
      }
  
    #Further splitting SampleDataLocWB into East and West datasets
    SampleDataLocWBEast=subset(SampleDataLocWB,Region=="East")
    SampleDataLocWBWest=subset(SampleDataLocWB,Region=="West")
    
      #Further splitting SampleDataLocWBEast and SampleDataLocWBWest into separate datasets 
      #which contain localities near and away from reservoirs respectively
      SampleDataLocWBEast0 = subset(SampleDataLocWBEast,WaterBodyExistence ==0)
      SampleDataLocWBEast1 = subset(SampleDataLocWBEast,WaterBodyExistence ==1)
      SampleDataLocWBWest0 = subset(SampleDataLocWBWest,WaterBodyExistence ==0)
      SampleDataLocWBWest1 = subset(SampleDataLocWBWest,WaterBodyExistence ==1)
  
  #Wrangling SampleData for ANOVA
  SampleDataAnova=SampleDataLoc%>%group_by(Date.Recorded,Region)%>%summarise(no_cases=sum(Cluster.Cases))
  
  #Wrangling SampleData for calculating Population CI
  SampleData = cast(SampleDataLoc, Date.Recorded~Region, fun.aggregate = sum, value='Cluster.Cases')
  SampleData$TotalCases = SampleData$Central+SampleData$East+SampleData$North+SampleData$`North East`+SampleData$West
  
    #Further Splitting SampleData into Normal and Epidemic Periods
    SampleDataEpidemic=subset(SampleData, grepl("13",SampleData$Date.Recorded)|grepl("14",SampleData$Date.Recorded))
    SampleDataNormal=subset(SampleData, !(grepl("13",SampleData$Date.Recorded)|grepl("14",SampleData$Date.Recorded)))
    
  #Monthly mean Temperature data
  SgTemp1217=subset(read.csv("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/surface-air-temperature-monthly-mean.csv"),Year %in% c(2012,2013,2014,2015,2016,2017))
  
  #Monthly mean Temperature data
  SgRain1217=subset(read.csv("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/rainfall-monthly-number-of-rain-days.csv"),Year %in% c(2012,2013,2014,2015,2016,2017))
  
  #SGDengue Dataset obtained by subsetting data of all diseases
  SgDengue = subset(read.csv("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/weekly-infectious-disease-bulletin-cases.csv"), 
                    disease %in% c("Dengue Fever","Dengue Haemorrhagic Fever"))
  
  #Calculating proportion of cases for each zone based on Cluster-wise data - (CI and hypothesis testing done later on)
  cluster_proportion = c(sum(SampleData$Central)/sum(SampleData$TotalCases),sum(SampleData$East)/sum(SampleData$TotalCases),
                             sum(SampleData$North)/sum(SampleData$TotalCases), sum(SampleData$`North East`)/sum(SampleData$TotalCases),
                             sum(SampleData$West)/sum(SampleData$TotalCases))
  
  #Breaking down total sgdengue cases into clusters - Sg and Zone wise weekly Dengue cases for 6 yr period
  SgDengue$Central = round(SgDengue$no._of_cases*cluster_proportion[1],digits=0) 
  SgDengue$East = round(SgDengue$no._of_cases*cluster_proportion[2],digits=0) 
  SgDengue$North = round(SgDengue$no._of_cases*cluster_proportion[3],digits=0) 
  SgDengue$NorthEast = round(SgDengue$no._of_cases*cluster_proportion[4],digits=0)  
  SgDengue$West = round(SgDengue$no._of_cases*cluster_proportion[5],digits=0)  
  
        #Adding month column to SGDengue table
      for (i in 1:length(SgDengue$Week))
      {
        if(SgDengue$Week[i] %in% c(1,2,3,4,5))
        {
          SgDengue$Month[i] = 1
        }   else if(SgDengue$Week[i] %in% c(6,7,8,9))
        {
          SgDengue$Month[i] = 2
        }   else if(SgDengue$Week[i] %in% c(10,11,12,13))
        {
          SgDengue$Month[i] = 3
        }   else if(SgDengue$Week[i] %in% c(14,15,16,17,18))
        {
          SgDengue$Month[i] = 4
        }   else if(SgDengue$Week[i] %in% c(19,20,21,22))
        {
          SgDengue$Month[i] = 5
        }   else if(SgDengue$Week[i] %in% c(23,24,25,26))
        {
          SgDengue$Month[i] = 6
        }   else if(SgDengue$Week[i] %in% c(27,28,29,30,31))
        {
          SgDengue$Month[i] = 7
        }   else if(SgDengue$Week[i] %in% c(32,33,34,35))
        {
          SgDengue$Month[i] = 8
        }   else if(SgDengue$Week[i] %in% c(36,37,38,39))
        {
          SgDengue$Month[i] = 9
        }   else if(SgDengue$Week[i] %in% c(40,41,42,43,44))
        {
          SgDengue$Month[i] = 10
        }   else if(SgDengue$Week[i] %in% c(45,46,47,48))
        {
          SgDengue$Month[i] = 11
        }   else
        {
          SgDengue$Month[i] = 12
        }
      }
   
  #Further splitting the sgdengue data according to two types of dengue fevers
  SgDengueFever=subset(SgDengue, disease %in% c("Dengue Fever"))
  SgDengueHaemorrhagicFever=subset(SgDengue, disease %in% c("Dengue Haemorrhagic Fever"))
  
  #PopulationData - for population vs sample hypothesis testing and DescStats_SgDengue Calculation
  PopulationData=SgDengue%>%group_by(Year,Month,Week)%>%summarise(no_cases=sum(no._of_cases),Central=sum(Central),
                                                      East=sum(East), North=sum(North), NorthEast=sum(NorthEast),
                                                      West=sum(West))
      #Further splitting PopulationData into Normal and Epidemic Periods
      PopulationDataEpidemic=subset(PopulationData, Year %in% c(2013,2014))
      PopulationDataNormal=subset(PopulationData, !Year %in% c(2013,2014))
  
  #Joining Temp, Rainfall, Dengue Data
  
  JoinedTable=merge(merge(PopulationData,SgTemp1217,by.x= c("Year","Month"),by.y=c("Year","Month"),sort=FALSE),
                    SgRain1217,by.x= c("Year","Month"),by.y=c("Year","Month"),sort=FALSE)
  
  #LRTable - Aggregating JoinedTable into monthwise numbers to perform linear regression
  
  LRTable=JoinedTable%>%group_by(Year, Month)%>% 
    summarise(no_cases = sum(no_cases), Central = sum(Central), East = sum(East), 
              North = sum(North), NorthEast = sum(NorthEast), West = sum(West), 
              mean_temp= sprintf("%0.1f",mean(as.numeric(mean_temp))), 
              no_of_rainy_days= sprintf("%0.0f",mean(as.numeric(no_of_rainy_days))))
  sapply(LRTable,mode)
  LRTable=transform(LRTable, mean_temp=as.numeric(mean_temp), no_of_rainy_days=as.numeric(no_of_rainy_days)) 
  
  #Singapore Visitors Data
  SgVisitors=read.csv("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/SG_Visitors_by_air.csv")
  SgVisitors$Month=NA
  SgVisitors$Year=NA
  SgVisitors[,c(ncol(SgVisitors)-1,ncol(SgVisitors))]=str_split_fixed(SgVisitors$Month.Year,"-",2)
  sapply(SgVisitors,mode)
  SgVisitors=transform(SgVisitors, Month=as.numeric(Month), Year=as.numeric(Year)) 
  SgVisitors=subset(SgVisitors,Year %in% c(2012,2013,2014,2015,2016,2017))
  
  #LRTable with Visitors - merging LRTable with Visitors data
  LRTablewithVisitors= merge(LRTable,SgVisitors,by.x= c("Month","Year"),by.y=c("Month","Year"),sort=FALSE)
  sapply(LRTablewithVisitors,mode)
  
#4. DESCRIPTIVE STATISTICS
  
  #Descriptive statistics of sample data  
    Desc_Stats_SampleData = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                     "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
    Desc_Stats_SampleData$Region= c("Total Sg Cases","Central","East","North","`North East`","West")
    
    Desc_Stats_SampleData[1,2:8]=c(summary(SampleData$TotalCases),getmode(SampleData$TotalCases))
    Desc_Stats_SampleData[2,2:8]=c(summary(SampleData$Central),  getmode(SampleData$Central))
    Desc_Stats_SampleData[3,2:8]=c(summary(SampleData$East),getmode(SampleData$East))
    Desc_Stats_SampleData[4,2:8]=c(summary(SampleData$North),getmode(SampleData$North))
    Desc_Stats_SampleData[5,2:8]=c(summary(SampleData$`North East`),getmode(SampleData$`North East`))
    Desc_Stats_SampleData[6,2:8]=c(summary(SampleData$West),getmode(SampleData$West))
    
    
    # Adding SD
    
    Desc_Stats_SampleData$stddev=c(sd(SampleData$TotalCases),sd(SampleData$Central),sd(SampleData$East),
                                 sd(SampleData$North),sd(SampleData$`North East`),sd(SampleData$West))
    
    # Adding Sum
    Desc_Stats_SampleData$sum=c(sum(SampleData$TotalCases),sum(SampleData$Central),sum(SampleData$East),
                              sum(SampleData$North),sum(SampleData$`North East`),sum(SampleData$West))
  

  #Descriptive analytics of total singapore dengue cases (6 years worth of data)
    Desc_Stats_SgDengue = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                     "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
    Desc_Stats_SgDengue$Region= c("Total Sg Cases","Central","East","North","NorthEast","West")
    
    Desc_Stats_SgDengue[1,2:8]=c(summary(PopulationData$no_cases),getmode(PopulationData$no_cases))
    Desc_Stats_SgDengue[2,2:8]=c(summary(PopulationData$Central),  getmode(PopulationData$Central))
    Desc_Stats_SgDengue[3,2:8]=c(summary(PopulationData$East),getmode(PopulationData$East))
    Desc_Stats_SgDengue[4,2:8]=c(summary(PopulationData$North),getmode(PopulationData$North))
    Desc_Stats_SgDengue[5,2:8]=c(summary(PopulationData$NorthEast),getmode(PopulationData$NorthEast))
    Desc_Stats_SgDengue[6,2:8]=c(summary(PopulationData$West),getmode(PopulationData$West))


      # Adding SD
    
        Desc_Stats_SgDengue$stddev=c(sd(PopulationData$no_cases),sd(PopulationData$Central),sd(PopulationData$East),
                                     sd(PopulationData$North),sd(PopulationData$NorthEast),sd(PopulationData$West))
      
      # Adding Sum
        Desc_Stats_SgDengue$sum=c(sum(PopulationData$no_cases),sum(PopulationData$Central),sum(PopulationData$East),
                                      sum(PopulationData$North),sum(PopulationData$NorthEast),sum(PopulationData$West))
        
    #Descriptive analytics of total singapore dengue fever cases
        Desc_Stats_SgDengueFever = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                         "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
        Desc_Stats_SgDengueFever$Region= c("Total Sg Cases","Central","East","North","NorthEast","West")
        
        Desc_Stats_SgDengueFever[1,2:8]=c(summary(SgDengueFever$no._of_cases),getmode(SgDengueFever$no._of_cases))
        Desc_Stats_SgDengueFever[2,2:8]=c(summary(SgDengueFever$Central),  getmode(SgDengueFever$Central))
        Desc_Stats_SgDengueFever[3,2:8]=c(summary(SgDengueFever$East),getmode(SgDengueFever$East))
        Desc_Stats_SgDengueFever[4,2:8]=c(summary(SgDengueFever$North),getmode(SgDengueFever$North))
        Desc_Stats_SgDengueFever[5,2:8]=c(summary(SgDengueFever$NorthEast),getmode(SgDengueFever$NorthEast))
        Desc_Stats_SgDengueFever[6,2:8]=c(summary(SgDengueFever$West),getmode(SgDengueFever$West))
        
        
        # Adding SD
        
        Desc_Stats_SgDengueFever$stddev=c(sd(SgDengueFever$no._of_cases),sd(SgDengueFever$Central),sd(SgDengueFever$East),
                                     sd(SgDengueFever$North),sd(SgDengueFever$NorthEast),sd(SgDengueFever$West))
        
        # Adding Sum
        Desc_Stats_SgDengueFever$sum=c(sum(SgDengueFever$no._of_cases),sum(SgDengueFever$Central),sum(SgDengueFever$East),
                                  sum(SgDengueFever$North),sum(SgDengueFever$NorthEast),sum(SgDengueFever$West))
        
        
    #Descriptive analytics of total singapore dengue haemorrhagic fever cases
        Desc_Stats_SgDengueHaemorrhagicFever = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                              "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
        Desc_Stats_SgDengueHaemorrhagicFever$Region= c("Total Sg Cases","Central","East","North","NorthEast","West")
        
        Desc_Stats_SgDengueHaemorrhagicFever[1,2:8]=c(summary(SgDengueHaemorrhagicFever$no._of_cases),getmode(SgDengueHaemorrhagicFever$no._of_cases))
        Desc_Stats_SgDengueHaemorrhagicFever[2,2:8]=c(summary(SgDengueHaemorrhagicFever$Central),  getmode(SgDengueHaemorrhagicFever$Central))
        Desc_Stats_SgDengueHaemorrhagicFever[3,2:8]=c(summary(SgDengueHaemorrhagicFever$East),getmode(SgDengueHaemorrhagicFever$East))
        Desc_Stats_SgDengueHaemorrhagicFever[4,2:8]=c(summary(SgDengueHaemorrhagicFever$North),getmode(SgDengueHaemorrhagicFever$North))
        Desc_Stats_SgDengueHaemorrhagicFever[5,2:8]=c(summary(SgDengueHaemorrhagicFever$NorthEast),getmode(SgDengueHaemorrhagicFever$NorthEast))
        Desc_Stats_SgDengueHaemorrhagicFever[6,2:8]=c(summary(SgDengueHaemorrhagicFever$West),getmode(SgDengueHaemorrhagicFever$West))
        
        
        # Adding SD
        
        Desc_Stats_SgDengueHaemorrhagicFever$stddev=c(sd(SgDengueHaemorrhagicFever$no._of_cases),sd(SgDengueHaemorrhagicFever$Central),sd(SgDengueHaemorrhagicFever$East),
                                          sd(SgDengueHaemorrhagicFever$North),sd(SgDengueHaemorrhagicFever$NorthEast),sd(SgDengueHaemorrhagicFever$West))
        
        # Adding Sum
        Desc_Stats_SgDengueHaemorrhagicFever$sum=c(sum(SgDengueHaemorrhagicFever$no._of_cases),sum(SgDengueHaemorrhagicFever$Central),sum(SgDengueHaemorrhagicFever$East),
                                       sum(SgDengueHaemorrhagicFever$North),sum(SgDengueHaemorrhagicFever$NorthEast),sum(SgDengueHaemorrhagicFever$West))
        
        

      #Descriptive statistics of sample data Normal Period  
        Desc_Stats_SampleDataNormal = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                           "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
        Desc_Stats_SampleDataNormal$Region= c("Total Sg Cases","Central","East","North","`North East`","West")
        
        Desc_Stats_SampleDataNormal[1,2:8]=c(summary(SampleDataNormal$TotalCases),getmode(SampleDataNormal$TotalCases))
        Desc_Stats_SampleDataNormal[2,2:8]=c(summary(SampleDataNormal$Central),  getmode(SampleDataNormal$Central))
        Desc_Stats_SampleDataNormal[3,2:8]=c(summary(SampleDataNormal$East),getmode(SampleDataNormal$East))
        Desc_Stats_SampleDataNormal[4,2:8]=c(summary(SampleDataNormal$North),getmode(SampleDataNormal$North))
        Desc_Stats_SampleDataNormal[5,2:8]=c(summary(SampleDataNormal$`North East`),getmode(SampleDataNormal$`North East`))
        Desc_Stats_SampleDataNormal[6,2:8]=c(summary(SampleDataNormal$West),getmode(SampleDataNormal$West))
        
        
        # Adding SD
        
        Desc_Stats_SampleDataNormal$stddev=c(sd(SampleDataNormal$TotalCases),sd(SampleDataNormal$Central),sd(SampleDataNormal$East),
                                       sd(SampleDataNormal$North),sd(SampleDataNormal$`North East`),sd(SampleDataNormal$West))
        
        # Adding Sum
        Desc_Stats_SampleDataNormal$sum=c(sum(SampleDataNormal$TotalCases),sum(SampleDataNormal$Central),sum(SampleDataNormal$East),
                                    sum(SampleDataNormal$North),sum(SampleDataNormal$`North East`),sum(SampleDataNormal$West))
        
      #Descriptive statistics of sample data Epidemic Period  
        Desc_Stats_SampleDataEpidemic = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                                 "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
        Desc_Stats_SampleDataEpidemic$Region= c("Total Sg Cases","Central","East","North","`North East`","West")
        
        Desc_Stats_SampleDataEpidemic[1,2:8]=c(summary(SampleDataEpidemic$TotalCases),getmode(SampleDataEpidemic$TotalCases))
        Desc_Stats_SampleDataEpidemic[2,2:8]=c(summary(SampleDataEpidemic$Central),  getmode(SampleDataEpidemic$Central))
        Desc_Stats_SampleDataEpidemic[3,2:8]=c(summary(SampleDataEpidemic$East),getmode(SampleDataEpidemic$East))
        Desc_Stats_SampleDataEpidemic[4,2:8]=c(summary(SampleDataEpidemic$North),getmode(SampleDataEpidemic$North))
        Desc_Stats_SampleDataEpidemic[5,2:8]=c(summary(SampleDataEpidemic$`North East`),getmode(SampleDataEpidemic$`North East`))
        Desc_Stats_SampleDataEpidemic[6,2:8]=c(summary(SampleDataEpidemic$West),getmode(SampleDataEpidemic$West))
        
        
        # Adding SD
        
        Desc_Stats_SampleDataEpidemic$stddev=c(sd(SampleDataEpidemic$TotalCases),sd(SampleDataEpidemic$Central),sd(SampleDataEpidemic$East),
                                             sd(SampleDataEpidemic$North),sd(SampleDataEpidemic$`North East`),sd(SampleDataEpidemic$West))
        
        # Adding Sum
        Desc_Stats_SampleDataEpidemic$sum=c(sum(SampleDataEpidemic$TotalCases),sum(SampleDataEpidemic$Central),sum(SampleDataEpidemic$East),
                                          sum(SampleDataEpidemic$North),sum(SampleDataEpidemic$`North East`),sum(SampleDataEpidemic$West))
        
        
      #Descriptive analytics of total singapore dengue cases (6 years worth of data) - Normal Period
        Desc_Stats_SgDengueNormal = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                         "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
        Desc_Stats_SgDengueNormal$Region= c("Total Sg Cases","Central","East","North","NorthEast","West")
        
        Desc_Stats_SgDengueNormal[1,2:8]=c(summary(PopulationDataNormal$no_cases),getmode(PopulationDataNormal$no_cases))
        Desc_Stats_SgDengueNormal[2,2:8]=c(summary(PopulationDataNormal$Central),  getmode(PopulationDataNormal$Central))
        Desc_Stats_SgDengueNormal[3,2:8]=c(summary(PopulationDataNormal$East),getmode(PopulationDataNormal$East))
        Desc_Stats_SgDengueNormal[4,2:8]=c(summary(PopulationDataNormal$North),getmode(PopulationDataNormal$North))
        Desc_Stats_SgDengueNormal[5,2:8]=c(summary(PopulationDataNormal$NorthEast),getmode(PopulationDataNormal$NorthEast))
        Desc_Stats_SgDengueNormal[6,2:8]=c(summary(PopulationDataNormal$West),getmode(PopulationDataNormal$West))
        
        
        # Adding SD
        
        Desc_Stats_SgDengueNormal$stddev=c(sd(PopulationDataNormal$no_cases),sd(PopulationDataNormal$Central),sd(PopulationDataNormal$East),
                                     sd(PopulationDataNormal$North),sd(PopulationDataNormal$NorthEast),sd(PopulationDataNormal$West))
        
        # Adding Sum
        Desc_Stats_SgDengueNormal$sum=c(sum(PopulationDataNormal$no_cases),sum(PopulationDataNormal$Central),sum(PopulationDataNormal$East),
                                  sum(PopulationDataNormal$North),sum(PopulationDataNormal$NorthEast),sum(PopulationDataNormal$West))
        
      #Descriptive analytics of total singapore dengue cases (6 years worth of data) - Epidemic Period
        Desc_Stats_SgDengueEpidemic = data.frame("Region"=rep(NA,6),"Min"=rep(NA,6),"Quartile1"=rep(NA,6), "Median"=rep(NA,6), 
                                               "Mean"=rep(NA,6), "Quartile3"=rep(NA,6), "Max"=rep(NA,6),"Mode"=rep(NA,6))
        Desc_Stats_SgDengueEpidemic$Region= c("Total Sg Cases","Central","East","North","NorthEast","West")
        
        Desc_Stats_SgDengueEpidemic[1,2:8]=c(summary(PopulationDataEpidemic$no_cases),getmode(PopulationDataEpidemic$no_cases))
        Desc_Stats_SgDengueEpidemic[2,2:8]=c(summary(PopulationDataEpidemic$Central),  getmode(PopulationDataEpidemic$Central))
        Desc_Stats_SgDengueEpidemic[3,2:8]=c(summary(PopulationDataEpidemic$East),getmode(PopulationDataEpidemic$East))
        Desc_Stats_SgDengueEpidemic[4,2:8]=c(summary(PopulationDataEpidemic$North),getmode(PopulationDataEpidemic$North))
        Desc_Stats_SgDengueEpidemic[5,2:8]=c(summary(PopulationDataEpidemic$NorthEast),getmode(PopulationDataEpidemic$NorthEast))
        Desc_Stats_SgDengueEpidemic[6,2:8]=c(summary(PopulationDataEpidemic$West),getmode(PopulationDataEpidemic$West))
        
        
        # Adding SD
        
        Desc_Stats_SgDengueEpidemic$stddev=c(sd(PopulationDataEpidemic$no_cases),sd(PopulationDataEpidemic$Central),sd(PopulationDataEpidemic$East),
                                           sd(PopulationDataEpidemic$North),sd(PopulationDataEpidemic$NorthEast),sd(PopulationDataEpidemic$West))
        
        # Adding Sum
        Desc_Stats_SgDengueEpidemic$sum=c(sum(PopulationDataEpidemic$no_cases),sum(PopulationDataEpidemic$Central),sum(PopulationDataEpidemic$East),
                                        sum(PopulationDataEpidemic$North),sum(PopulationDataEpidemic$NorthEast),sum(PopulationDataEpidemic$West))
        
#5. INFERENTIAL STATISTICS
  
#Linear Regression with LRTable
  
  #Linear Regression stepwise without transformation
      
  fitstep1=lm(no_cases~Month+mean_temp+no_of_rainy_days,data=LRTable)
  step1=stepAIC(fitstep1,direction='both')
  fit1=lm(no_cases~mean_temp, data = LRTable)
  cor(LRTable$no_cases, LRTable$mean_temp)
  summary(fit1)
  anova(fit1)
  # plot(fit1)
  
  #Linear Regression stepwise with transformation of dependent variable using trafo package 
  #Best results
  
  # assumptions(fit1)
  fit_trafo=trafo_lm(fit1,trafo="sqrtshift")
  diagnostics(fit_trafo)
  summary(fit_trafo)
  # plot(fit_trafo)
  
  #Generalised Linear Regression stepwise using poisson distribution
  
  fit_glm=glm(no_cases~mean_temp+no_of_rainy_days, family=poisson, data = LRTable)
  step_glm=stepAIC(fit_glm,direction = 'both')
  summary(fit_glm)
  anova(fit_glm)
  # plot(fit_glm)
  cor(LRTable$no_of_rainy_days,LRTable$mean_temp)
  cor(LRTable$no_of_rainy_days,LRTable$no_cases)
  
  #Generalised Linear Regression stepwise using poisson distribution with transformation of dependent variable using trafo package
  
  # assumptions(step_glm)
  step_glm_trafo=trafo_lm(step_glm,trafo="sqrtshift")
  summary(step_glm_trafo)
  # plot(step_glm_trafo)
  
  #Transforming predictors using log
  
  LRTable$LNmean_temp=log(LRTable$mean_temp)
  LRTable$LNno_of_rainy_days=log(LRTable$no_of_rainy_days)
  
  #Linear Regression stepwise with log transformed predictors
  
  fitstep2=lm(no_cases~Month+LNmean_temp+LNno_of_rainy_days,data=LRTable)
  step2=stepAIC(fitstep2,direction='both')
  fit2=lm(no_cases~LNmean_temp, data = LRTable)
  cor(LRTable$no_cases, LRTable$LNmean_temp)
  summary(fit2)
  anova(fit2)
  # plot(fit2)
  
  #Linear Regression stepwise with log transformed predictors and transformation of dependent variable using trafo package
  #Second best
  
  # assumptions(fit2)
  fit_trafo2=trafo_lm(fit2,trafo="sqrtshift")
  diagnostics(fit_trafo2)
  summary(fit_trafo2)
  # plot(fit_trafo2)
  
  #Generalised Linear Regression stepwise using poisson distribution with log transformed predictors
  
  fit_glm2=glm(no_cases~LNmean_temp+LNno_of_rainy_days, family=poisson, data = LRTable)
  step_glm=stepAIC(fit_glm2,direction = 'both')
  summary(fit_glm2)
  anova(fit_glm2)
  # plot(fit_glm2)
  cor(LRTable$LNno_of_rainy_days,LRTable$LNmean_temp)
  cor(LRTable$LNno_of_rainy_days,LRTable$no_cases)
  cor(LRTable$LNmean_temp,LRTable$no_cases)
  

    #Residuals vs fitted plot shows that there's positive autocorrelation, requiring a time series analysis;
    #There might also be leakage of prediction in the following ways:
        # A missing variable
        # A missing higher-order term of a variable in the model to explain the curvature
        # A missing interaction between terms already in the model

#Linear Regression with LRTablewithVisitors - BEST OF ALL REGRESSION
  ScaledLRTablewithVisitors=data.frame(scale(LRTablewithVisitors[,c(3,9,10,14:ncol(LRTablewithVisitors))]))
  
  ccs <- as.matrix(ScaledLRTablewithVisitors)
  rcx=cor(ccs, method=c("pearson"))
  corrplot(rcx,method="circle", type="lower",diag = FALSE,tl.cex = 0.5)
 

  ScaledLRTablewithVisitors=ScaledLRTablewithVisitors[,-which(names(ScaledLRTablewithVisitors) %in% c("Spain","Germany","Vietnam","United.Kingdom",
                                                                                                      "France","Bangladesh","Italy","Netherlands",
                                                                                                      "Belgium...Luxembourg","Canada","Hong.Kong",
                                                                                                      "China","CIS","no_of_rainy_days","Philippines",
                                                                                                      "Republic.of.South.Africa","Taiwan","Brunei",
                                                                                                      "Greece","Finland","Kuwait","Malaysia",
                                                                                                      "Korea"))]
  
  xnam=paste("ScaledLRTablewithVisitors$",colnames(ScaledLRTablewithVisitors)[2:ncol(ScaledLRTablewithVisitors)],sep="",collapse=" + ")
  fml=as.formula(paste("ScaledLRTablewithVisitors$no_cases ~",xnam))
  
  fitstep3=lm(fml)
  step3=stepAIC(fitstep3,direction='both')
  summary(step3)
  anova(step3)
  # plot(step3)
  cor(ScaledLRTablewithVisitors$Myanmar,ScaledLRTablewithVisitors$no_cases)
  cor(ScaledLRTablewithVisitors$New.Zealand,ScaledLRTablewithVisitors$no_cases)
  cor(ScaledLRTablewithVisitors$Indonesia,ScaledLRTablewithVisitors$no_cases)
  cor(ScaledLRTablewithVisitors$Ireland,ScaledLRTablewithVisitors$no_cases)
  cor(ScaledLRTablewithVisitors$Switzerland,ScaledLRTablewithVisitors$no_cases)
  cor(ScaledLRTablewithVisitors$Australia,ScaledLRTablewithVisitors$no_cases)
  cor(ScaledLRTablewithVisitors$mean_temp,ScaledLRTablewithVisitors$no_cases)
  
  
  #Sample Vs Population - All 15 rows of sample data
    #n=15; Yet, assume normal distribution to derive CI, since bootstrapping "normal" algorithm gave similar results as other algorithms
      #Population Mean
        #Population mean - Point Estimate
        
        Mean=data.frame("Region"=c('Total','Central','East','North','NorthEast','West'), "SampleMean"=c(mean(SampleData$TotalCases),
                                                                                                       mean(SampleData$Central),
                                                                                                       mean(SampleData$East),
                                                                                                       mean(SampleData$North),
                                                                                                       mean(SampleData$`North East`), 
                                                                                                       mean(SampleData$West)))
        
        #95% & 90% CI for Population Mean
          xbar=Mean$SampleMean
          s=c(sd(SampleData$TotalCases),sd(SampleData$Central),sd(SampleData$East),sd(SampleData$North),sd(SampleData$`North East`),
              sd(SampleData$West))
          Mean$LowerCI5=NA
          Mean$UpperCI5=NA
          Mean$LowerCI10=NA
          Mean$UpperCI10=NA
          n= nrow(SampleData)
          alpha=0.05
          for(i in 1:6)
          {
            Mean$LowerCI5[i]=xbar[i]-qt(1-alpha/2,n-1)*s/sqrt(n)
            Mean$UpperCI5[i]=xbar[i]+qt(1-alpha/2,n-1)*s/sqrt(n)
          }
          alpha=0.1
          for(i in 1:6)
          {
            Mean$LowerCI10[i]=xbar[i]-qt(1-alpha/2,n-1)*s/sqrt(n)
            Mean$UpperCI10[i]=xbar[i]+qt(1-alpha/2,n-1)*s/sqrt(n)
          }
        
        
        #Hypothesis Testing for Population Mean: H0: Mu=PopMean[i], n<30,Single Population
          #alpha for two-tailed test is taken to be twice the alpha used for one-tailed tests   
          n= nrow(SampleData)
          xbar=Mean$SampleMean
          Mean$PopMean=c(mean(PopulationData$no_cases),mean(PopulationData$Central),mean(PopulationData$East),
                    mean(PopulationData$North),mean(PopulationData$NorthEast),mean(PopulationData$West))
          mu=Mean$PopMean
          s=c(sd(SampleData$TotalCases),sd(SampleData$Central),sd(SampleData$East),sd(SampleData$North),sd(SampleData$`North East`),
             sd(SampleData$West))
          Mean$pValueTwoTailed=NA
          Mean$pValueUpperTailed=NA
          Mean$pValueLowerTailed=NA
          Mean$NullHypTwoTailed=NA
          Mean$NullHypUpperTailed=NA
          Mean$NullHypLowerTailed=NA
          
          
          alpha=0.01
          
          for(i in 1:6)
          {
            t=(xbar[i]-mu[i])/(s[i]/sqrt(n))
            Mean$pValueTwoTailed[i]=2*(1-pt(t,n-1))
            Mean$pValueUpperTailed[i]=(1-pt(t,n-1))
            Mean$pValueLowerTailed[i]=pt(t,n-1)
            if(Mean$pValueTwoTailed[i]<2*alpha)
            {
                Mean$NullHypTwoTailed[i]="Reject"
                
            }else
              {
               Mean$NullHypTwoTailed[i]="Do Not Reject"
              }
            if(Mean$pValueUpperTailed[i]<alpha)
            {
              Mean$NullHypUpperTailed[i]="Reject"
            }else
              {
                Mean$NullHypUpperTailed[i]="Do Not Reject"
              }
            if(Mean$pValueLowerTailed[i]<alpha)
            {
              Mean$NullHypLowerTailed[i]="Reject"
            }else
            {
              Mean$NullHypLowerTailed[i]="Do Not Reject"
            }    
              
          }
      #Population Proportion
        #Population Proportion - Point Estimate
        
        Proportion = data.frame("Region"=c('Central','East','North','NorthEast','West'), 
                                "SampleProportion"=c(sum(SampleData$Central)/sum(SampleData$TotalCases),
                                                     sum(SampleData$East)/sum(SampleData$TotalCases),
                                                     sum(SampleData$North)/sum(SampleData$TotalCases), 
                                                     sum(SampleData$`North East`)/sum(SampleData$TotalCases),
                                                     sum(SampleData$West)/sum(SampleData$TotalCases)))
        
        #95% & 90% CI for Population Proportion
        p=Proportion$SampleProportion
        n= nrow(SampleData)
        Proportion$LowerCI5=NA
        Proportion$UpperCI5=NA
        Proportion$LowerCI10=NA
        Proportion$UpperCI10=NA
        alpha=0.05
        for(i in 1:5)
        {
          Proportion$LowerCI5[i]=p[i]-qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
          Proportion$UpperCI5[i]=p[i]+qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
        }
        alpha=0.1
        for(i in 1:5)
        {
          Proportion$LowerCI10[i]=p[i]-qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
          Proportion$UpperCI10[i]=p[i]+qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
        }
        
        
        #Hypothesis Testing for Population Proportion: H0: Mu=SampleProportion[i], n<30,Single Population
        n= nrow(SampleData)
        p=Proportion$SampleProportion
        pi=Proportion$SampleProportion
        Proportion$pValueTwoTailed=NA
        Proportion$pValueUpperTailed=NA
        Proportion$pValueLowerTailed=NA
        Proportion$NullHypTwoTailed=NA
        Proportion$NullHypUpperTailed=NA
        Proportion$NullHypLowerTailed=NA
        
        
        alpha=0.01
        
        for(i in 1:5)
        {
          z=(p[i]-pi[i])/sqrt(pi[i]*(1-pi[i])/n)
          Proportion$pValueTwoTailed[i]=2*(1-pnorm(z))
          Proportion$pValueUpperTailed[i]=1-pnorm(z)
          Proportion$pValueLowerTailed[i]=pnorm(z)
          if(Proportion$pValueTwoTailed[i]<2*alpha)
          {
            Proportion$NullHypTwoTailed[i]="Reject"
            
          }else
          {
            Proportion$NullHypTwoTailed[i]="Do Not Reject"
          }
          if(Proportion$pValueUpperTailed[i]<alpha)
          {
            Proportion$NullHypUpperTailed[i]="Reject"
          }else
          {
            Proportion$NullHypUpperTailed[i]="Do Not Reject"
          }
          if(Proportion$pValueLowerTailed[i]<alpha)
          {
            Proportion$NullHypLowerTailed[i]="Reject"
          }else
          {
            Proportion$NullHypLowerTailed[i]="Do Not Reject"
          }    
          
        }
        
  #SampleDataNormal - All  rows of sample data
    #n=8; Yet, assume normal distribution to derive CI, since bootstrapping "normal" algorithm gave similar results as other algorithms
      #Population Mean Normal
        #Population mean Normal - Point Estimate
        
        MeanNormal=data.frame("Region"=c('Total','Central','East','North','NorthEast','West'), "SampleMeanNormal"=c(mean(SampleDataNormal$TotalCases),
                                                                                                        mean(SampleDataNormal$Central),
                                                                                                        mean(SampleDataNormal$East),
                                                                                                        mean(SampleDataNormal$North),
                                                                                                        mean(SampleDataNormal$`North East`), 
                                                                                                        mean(SampleDataNormal$West)))
        
        #95% & 90% CI for Population Mean Normal
        xbar=MeanNormal$SampleMeanNormal
        s=c(sd(SampleDataNormal$TotalCases),sd(SampleDataNormal$Central),sd(SampleDataNormal$East),sd(SampleDataNormal$North),sd(SampleDataNormal$`North East`),
            sd(SampleDataNormal$West))
        MeanNormal$LowerCI5=NA
        MeanNormal$UpperCI5=NA
        MeanNormal$LowerCI10=NA
        MeanNormal$UpperCI10=NA
        n= nrow(SampleDataNormal)
        alpha=0.05
        for(i in 1:6)
        {
          MeanNormal$LowerCI5[i]=xbar[i]-qt(1-alpha/2,n-1)*s/sqrt(n)
          MeanNormal$UpperCI5[i]=xbar[i]+qt(1-alpha/2,n-1)*s/sqrt(n)
        }
        alpha=0.1
        for(i in 1:6)
        {
          MeanNormal$LowerCI10[i]=xbar[i]-qt(1-alpha/2,n-1)*s/sqrt(n)
          MeanNormal$UpperCI10[i]=xbar[i]+qt(1-alpha/2,n-1)*s/sqrt(n)
        }
        
        #Hypothesis Testing for Population Mean Normal: H0: Mu=PopMeanNormal[i], n<30,Single Population
        n=nrow(SampleDataNormal)
        xbar=MeanNormal$SampleMeanNormal
        MeanNormal$PopMeanNormal=c(mean(PopulationDataNormal$no_cases),mean(PopulationDataNormal$Central),mean(PopulationDataNormal$East),
                       mean(PopulationDataNormal$North),mean(PopulationDataNormal$NorthEast),mean(PopulationDataNormal$West))
        mu=MeanNormal$PopMeanNormal
        s=c(sd(SampleDataNormal$TotalCases),sd(SampleDataNormal$Central),sd(SampleDataNormal$East),sd(SampleDataNormal$North),sd(SampleDataNormal$`North East`),
            sd(SampleDataNormal$West))
        MeanNormal$pValueTwoTailed=NA
        MeanNormal$pValueUpperTailed=NA
        MeanNormal$pValueLowerTailed=NA
        MeanNormal$NullHypTwoTailed=NA
        MeanNormal$NullHypUpperTailed=NA
        MeanNormal$NullHypLowerTailed=NA
        
        
        alpha=0.01
        
        for(i in 1:6)
        {
          t=(xbar[i]-mu[i])/(s[i]/sqrt(n))
          MeanNormal$pValueTwoTailed[i]=2*(1-pt(t,n-1))
          MeanNormal$pValueUpperTailed[i]=(1-pt(t,n-1))
          MeanNormal$pValueLowerTailed[i]=pt(t,n-1)
          if(MeanNormal$pValueTwoTailed[i]<2*alpha)
          {
            MeanNormal$NullHypTwoTailed[i]="Reject"
            
          }else
          {
            MeanNormal$NullHypTwoTailed[i]="Do Not Reject"
          }
          if(MeanNormal$pValueUpperTailed[i]<alpha)
          {
            MeanNormal$NullHypUpperTailed[i]="Reject"
          }else
          {
            MeanNormal$NullHypUpperTailed[i]="Do Not Reject"
          }
          if(MeanNormal$pValueLowerTailed[i]<alpha)
          {
            MeanNormal$NullHypLowerTailed[i]="Reject"
          }else
          {
            MeanNormal$NullHypLowerTailed[i]="Do Not Reject"
          }    
          
        }
        
      #Population Proportion
        #Population Proportion Normal- Point Estimate
        
        ProportionNormal = data.frame("Region"=c('Central','East','North','NorthEast','West'), 
                                "SampleProportionNormal"=c(sum(SampleDataNormal$Central)/sum(SampleDataNormal$TotalCases),
                                                     sum(SampleDataNormal$East)/sum(SampleDataNormal$TotalCases),
                                                     sum(SampleDataNormal$North)/sum(SampleDataNormal$TotalCases), 
                                                     sum(SampleDataNormal$`North East`)/sum(SampleDataNormal$TotalCases),
                                                     sum(SampleDataNormal$West)/sum(SampleDataNormal$TotalCases)))
        
        #95% & 90% CI for Population Proportion
        p=ProportionNormal$SampleProportionNormal
        n= nrow(SampleDataNormal)
        ProportionNormal$LowerCI5=NA
        ProportionNormal$UpperCI5=NA
        ProportionNormal$LowerCI10=NA
        ProportionNormal$UpperCI10=NA
        alpha=0.05
        for(i in 1:5)
        {
          ProportionNormal$LowerCI5[i]=p[i]-qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
          ProportionNormal$UpperCI5[i]=p[i]+qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
        }
        alpha=0.1
        for(i in 1:5)
        {
          ProportionNormal$LowerCI10[i]=p[i]-qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
          ProportionNormal$UpperCI10[i]=p[i]+qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
        }
        
        #Hypothesis Testing for Population Proportion Normal: H0: Mu=SampleProportionNormal[i], n<30,Single Population
        n=nrow(SampleDataNormal)
        p=ProportionNormal$SampleProportionNormal
        pi=ProportionNormal$SampleProportionNormal
        ProportionNormal$pValueTwoTailed=NA
        ProportionNormal$pValueUpperTailed=NA
        ProportionNormal$pValueLowerTailed=NA
        ProportionNormal$NullHypTwoTailed=NA
        ProportionNormal$NullHypUpperTailed=NA
        ProportionNormal$NullHypLowerTailed=NA
        
        
        alpha=0.01
        
        for(i in 1:5)
        {
          z=(p[i]-pi[i])/sqrt(pi[i]*(1-pi[i])/n)
          ProportionNormal$pValueTwoTailed[i]=2*(1-pnorm(z))
          ProportionNormal$pValueUpperTailed[i]=1-pnorm(z)
          ProportionNormal$pValueLowerTailed[i]=pnorm(z)
          if(ProportionNormal$pValueTwoTailed[i]<2*alpha)
          {
            ProportionNormal$NullHypTwoTailed[i]="Reject"
            
          }else
          {
            ProportionNormal$NullHypTwoTailed[i]="Do Not Reject"
          }
          if(ProportionNormal$pValueUpperTailed[i]<alpha)
          {
            ProportionNormal$NullHypUpperTailed[i]="Reject"
          }else
          {
            ProportionNormal$NullHypUpperTailed[i]="Do Not Reject"
          }
          if(ProportionNormal$pValueLowerTailed[i]<alpha)
          {
            ProportionNormal$NullHypLowerTailed[i]="Reject"
          }else
          {
            ProportionNormal$NullHypLowerTailed[i]="Do Not Reject"
          }    
          
        }
        
  #SampleDataEpidemic - All  rows of sample data
    #n=7; Yet, assume normal distribution to derive CI, since bootstrapping "normal" algorithm gave similar results as other algorithms
      #Population Mean Epidemic
        #Population mean Epidemic - Point Estimate
        
        MeanEpidemic=data.frame("Region"=c('Total','Central','East','North','NorthEast','West'), "SampleMeanEpidemic"=c(mean(SampleDataEpidemic$TotalCases),
                                                                                                                    mean(SampleDataEpidemic$Central),
                                                                                                                    mean(SampleDataEpidemic$East),
                                                                                                                    mean(SampleDataEpidemic$North),
                                                                                                                    mean(SampleDataEpidemic$`North East`), 
                                                                                                                    mean(SampleDataEpidemic$West)))
        
        #95% & 90% CI for Population Mean Epidemic
        xbar=MeanEpidemic$SampleMeanEpidemic
        s=c(sd(SampleDataEpidemic$TotalCases),sd(SampleDataEpidemic$Central),sd(SampleDataEpidemic$East),sd(SampleDataEpidemic$North),sd(SampleDataEpidemic$`North East`),
            sd(SampleDataEpidemic$West))
        MeanEpidemic$LowerCI5=NA
        MeanEpidemic$UpperCI5=NA
        MeanEpidemic$LowerCI10=NA
        MeanEpidemic$UpperCI10=NA
        n= nrow(SampleDataEpidemic)
        alpha=0.05
        for(i in 1:6)
        {
          MeanEpidemic$LowerCI5[i]=xbar[i]-qt(1-alpha/2,n-1)*s/sqrt(n)
          MeanEpidemic$UpperCI5[i]=xbar[i]+qt(1-alpha/2,n-1)*s/sqrt(n)
        }
        alpha=0.1
        for(i in 1:6)
        {
          MeanEpidemic$LowerCI10[i]=xbar[i]-qt(1-alpha/2,n-1)*s/sqrt(n)
          MeanEpidemic$UpperCI10[i]=xbar[i]+qt(1-alpha/2,n-1)*s/sqrt(n)
        }
        
        #Hypothesis Testing for Population Mean Epidemic: H0: Mu=PopMeanEpidemic[i], n<30,Single Population
        n=nrow(SampleDataEpidemic)
        xbar=MeanEpidemic$SampleMeanEpidemic
        MeanEpidemic$PopMeanEpidemic=c(mean(PopulationDataEpidemic$no_cases),mean(PopulationDataEpidemic$Central),mean(PopulationDataEpidemic$East),
                                   mean(PopulationDataEpidemic$North),mean(PopulationDataEpidemic$NorthEast),mean(PopulationDataEpidemic$West))

        mu=MeanEpidemic$PopMeanEpidemic
        s=c(sd(SampleDataEpidemic$TotalCases),sd(SampleDataEpidemic$Central),sd(SampleDataEpidemic$East),sd(SampleDataEpidemic$North),sd(SampleDataEpidemic$`North East`),
            sd(SampleDataEpidemic$West))
        MeanEpidemic$pValueTwoTailed=NA
        MeanEpidemic$pValueUpperTailed=NA
        MeanEpidemic$pValueLowerTailed=NA
        MeanEpidemic$NullHypTwoTailed=NA
        MeanEpidemic$NullHypUpperTailed=NA
        MeanEpidemic$NullHypLowerTailed=NA
        
        
        alpha=0.01
        
        for(i in 1:6)
        {
          t=(xbar[i]-mu[i])/(s[i]/sqrt(n))
          MeanEpidemic$pValueTwoTailed[i]=2*(1-pt(t,n-1))
          MeanEpidemic$pValueUpperTailed[i]=(1-pt(t,n-1))
          MeanEpidemic$pValueLowerTailed[i]=pt(t,n-1)
          if(MeanEpidemic$pValueTwoTailed[i]<2*alpha)
          {
            MeanEpidemic$NullHypTwoTailed[i]="Reject"
            
          }else
          {
            MeanEpidemic$NullHypTwoTailed[i]="Do Not Reject"
          }
          if(MeanEpidemic$pValueUpperTailed[i]<alpha)
          {
            MeanEpidemic$NullHypUpperTailed[i]="Reject"
          }else
          {
            MeanEpidemic$NullHypUpperTailed[i]="Do Not Reject"
          }
          if(MeanEpidemic$pValueLowerTailed[i]<alpha)
          {
            MeanEpidemic$NullHypLowerTailed[i]="Reject"
          }else
          {
            MeanEpidemic$NullHypLowerTailed[i]="Do Not Reject"
          }    
          
        }
      #Population Proportion
        #Population Proportion Epidemic- Point Estimate
        
        ProportionEpidemic = data.frame("Region"=c('Central','East','North','NorthEast','West'), 
                                      "SampleProportionEpidemic"=c(sum(SampleDataEpidemic$Central)/sum(SampleDataEpidemic$TotalCases),
                                                                 sum(SampleDataEpidemic$East)/sum(SampleDataEpidemic$TotalCases),
                                                                 sum(SampleDataEpidemic$North)/sum(SampleDataEpidemic$TotalCases), 
                                                                 sum(SampleDataEpidemic$`North East`)/sum(SampleDataEpidemic$TotalCases),
                                                                 sum(SampleDataEpidemic$West)/sum(SampleDataEpidemic$TotalCases)))
        
        #95% & 90% CI for Population Proportion
        p=ProportionEpidemic$SampleProportionEpidemic
        n= nrow(SampleDataEpidemic)
        ProportionEpidemic$LowerCI5=NA
        ProportionEpidemic$UpperCI5=NA
        ProportionEpidemic$LowerCI10=NA
        ProportionEpidemic$UpperCI10=NA
        alpha=0.05
        for(i in 1:5)
        {
          ProportionEpidemic$LowerCI5[i]=p[i]-qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
          ProportionEpidemic$UpperCI5[i]=p[i]+qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
        }
        alpha=0.1
        for(i in 1:5)
        {
          ProportionEpidemic$LowerCI10[i]=p[i]-qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
          ProportionEpidemic$UpperCI10[i]=p[i]+qnorm(1-alpha/2)*sqrt(p[i]*(1-p[i])/n)
        }
        
        #Hypothesis Testing for Population Proportion Epidemic: H0: Mu=SampleProportionEpidemic[i], n<30,Single Population
        n=nrow(SampleDataEpidemic)
        p=ProportionEpidemic$SampleProportionEpidemic
        pi=ProportionEpidemic$SampleProportionEpidemic
        ProportionEpidemic$pValueTwoTailed=NA
        ProportionEpidemic$pValueUpperTailed=NA
        ProportionEpidemic$pValueLowerTailed=NA
        ProportionEpidemic$NullHypTwoTailed=NA
        ProportionEpidemic$NullHypUpperTailed=NA
        ProportionEpidemic$NullHypLowerTailed=NA
        
        
        alpha=0.01
        
        for(i in 1:5)
        {
          z=(p[i]-pi[i])/sqrt(pi[i]*(1-pi[i])/n)
          ProportionEpidemic$pValueTwoTailed[i]=2*(1-pnorm(z))
          ProportionEpidemic$pValueUpperTailed[i]=1-pnorm(z)
          ProportionEpidemic$pValueLowerTailed[i]=pnorm(z)
          if(ProportionEpidemic$pValueTwoTailed[i]<2*alpha)
          {
            ProportionEpidemic$NullHypTwoTailed[i]="Reject"
            
          }else
          {
            ProportionEpidemic$NullHypTwoTailed[i]="Do Not Reject"
          }
          if(ProportionEpidemic$pValueUpperTailed[i]<alpha)
          {
            ProportionEpidemic$NullHypUpperTailed[i]="Reject"
          }else
          {
            ProportionEpidemic$NullHypUpperTailed[i]="Do Not Reject"
          }
          if(ProportionEpidemic$pValueLowerTailed[i]<alpha)
          {
            ProportionEpidemic$NullHypLowerTailed[i]="Reject"
          }else
          {
            ProportionEpidemic$NullHypLowerTailed[i]="Do Not Reject"
          }    
          
        }
        
  #Outbreak Simulation Statistics
    OutbreakStats=data.frame("Region"=c('Total','Central','East','North','NorthEast','West'))
    AdjLCI=data.frame("Region"=c('Total','Central','East','North','NorthEast','West'))
    for(i in 1:6)
    {
      if(MeanNormal$LowerCI5[i]<0)
      {
        AdjLCI$MeanNormalLCI[i]=0
      }else
        {
          AdjLCI$MeanNormalLCI[i]=MeanNormal$LowerCI5[i]
        } 
      
      if(MeanEpidemic$LowerCI5[i]<0)
      {
        AdjLCI$MeanEpidemicLCI[i]=0
      }else
      {
        AdjLCI$MeanEpidemicLCI[i]=MeanEpidemic$LowerCI5[i]
      } 
    }
  
    
    OutbreakStats$PctIncEpMeanBlMean=((MeanEpidemic$SampleMeanEpidemic-MeanNormal$SampleMeanNormal)/MeanNormal$SampleMeanNormal)*100
    OutbreakStats$PctIncEpMeanBlLCI=((MeanEpidemic$SampleMeanEpidemic-AdjLCI$MeanNormalLCI)/(AdjLCI$MeanNormalLCI+1))*100
    OutbreakStats$PctIncEpMeanBlUCI=((MeanEpidemic$SampleMeanEpidemic-MeanNormal$UpperCI5)/MeanNormal$UpperCI5)*100   
    OutbreakStats$PctIncEpLCIBlMean=((AdjLCI$MeanEpidemicLCI-MeanNormal$SampleMeanNormal)/MeanNormal$SampleMeanNormal)*100
    OutbreakStats$PctIncEpLCIBlLCI=((AdjLCI$MeanEpidemicLCI-AdjLCI$MeanNormalLCI)/(AdjLCI$MeanNormalLCI+1))*100
    OutbreakStats$PctIncEpLCIBlUCI=((AdjLCI$MeanEpidemicLCI-MeanNormal$UpperCI5)/MeanNormal$UpperCI5)*100
    OutbreakStats$PctIncEpUCIBlMean=((MeanEpidemic$UpperCI5-MeanNormal$SampleMeanNormal)/MeanNormal$SampleMeanNormal)*100
    OutbreakStats$PctIncEpUCIBlLCI=((MeanEpidemic$UpperCI5-AdjLCI$MeanNormalLCI)/(AdjLCI$MeanNormalLCI+1))*100
    OutbreakStats$PctIncEpUCIBlUCI=((MeanEpidemic$UpperCI5-MeanNormal$UpperCI5)/MeanNormal$UpperCI5)*100
    
  #Hypothesis Testing - Difference in means between outbreak and normal period: 2 populations, n<30, unknown and unequal population variances
    #H0: Mu1-Mu2<=0
    #H1: Mu1-Mu2>0
    HypothesisEpidemicVsNormal=data.frame("Region"=c('Total','Central','East','North','NorthEast','West'))
    xbar=c(mean(SampleDataEpidemic$TotalCases),mean(SampleDataEpidemic$Central),
           mean(SampleDataEpidemic$East),mean(SampleDataEpidemic$North),
           mean(SampleDataEpidemic$`North East`),mean(SampleDataEpidemic$West))
    ybar=c(mean(SampleDataNormal$TotalCases),mean(SampleDataNormal$Central),
           mean(SampleDataNormal$East),mean(SampleDataNormal$North),
           mean(SampleDataNormal$`North East`),mean(SampleDataNormal$West))
    HypothesisEpidemicVsNormal$MeanDifferenceEpidemicNormal=xbar-ybar
    sx=c(sd(SampleDataEpidemic$TotalCases),sd(SampleDataEpidemic$Central),
           sd(SampleDataEpidemic$East),sd(SampleDataEpidemic$North),
           sd(SampleDataEpidemic$`North East`),sd(SampleDataEpidemic$West))
    sy=c(sd(SampleDataNormal$TotalCases),sd(SampleDataNormal$Central),
         sd(SampleDataNormal$East),sd(SampleDataNormal$North),
         sd(SampleDataNormal$`North East`),sd(SampleDataNormal$West))
    nx=nrow(SampleDataEpidemic)
    ny=nrow(SampleDataNormal)
    alpha=0.1
    HypothesisEpidemicVsNormal$t=NA
    taTwoTailed=NA
    taUpperTailed=NA
    taLowerTailed=NA
    HypothesisEpidemicVsNormal$pvalueTwoTailed=NA
    HypothesisEpidemicVsNormal$pvalueUpperTailed=NA
    HypothesisEpidemicVsNormal$pvalueLowerTailed=NA
    HypothesisEpidemicVsNormal$NullHypTwoTailed=NA
    HypothesisEpidemicVsNormal$NullHypUpperTailed=NA
    HypothesisEpidemicVsNormal$NullHypLowerTailed=NA
    
    for(i in 1:6)
    {
      HypothesisEpidemicVsNormal$t[i]=(xbar[i]-ybar[i])/sqrt(sx[i]^2/nx+sy[i]^2/ny)
      taTwoTailed[i]=qt(1-2*alpha/2,nx+ny-2)
      taUpperTailed[i]=qt(1-alpha,nx+ny-2)
      taLowerTailed[i]=qt(alpha,nx+ny-2)
      if(HypothesisEpidemicVsNormal$t[i]<0)
      {
        HypothesisEpidemicVsNormal$pvalueTwoTailed[i]= 2*(pt(HypothesisEpidemicVsNormal$t[i],nx+ny-2))
      }else
      {
        HypothesisEpidemicVsNormal$pvalueTwoTailed[i]= 2*(1-pt(HypothesisEpidemicVsNormal$t[i],nx+ny-2))
      }
      
      HypothesisEpidemicVsNormal$pvalueUpperTailed[i]=(1-pt(HypothesisEpidemicVsNormal$t[i],nx+ny-2))
      HypothesisEpidemicVsNormal$pvalueLowerTailed[i]=pt(HypothesisEpidemicVsNormal$t[i],nx+ny-2)
      
      
      if(HypothesisEpidemicVsNormal$pvalueTwoTailed[i]<(2*alpha))
      {
        HypothesisEpidemicVsNormal$NullHypTwoTailed[i]="Reject"
      }else
      {
        HypothesisEpidemicVsNormal$NullHypTwoTailed[i]="Do Not Reject"
      }
      
      if(HypothesisEpidemicVsNormal$pvalueUpperTailed[i]<alpha)
      {
        HypothesisEpidemicVsNormal$NullHypUpperTailed[i]="Reject"
      }else
      {
        HypothesisEpidemicVsNormal$NullHypUpperTailed[i]="Do Not Reject"
      }
      if(HypothesisEpidemicVsNormal$pvalueLowerTailed[i]<alpha)
      {
        HypothesisEpidemicVsNormal$NullHypLowerTailed[i]="Reject"
      }else
      {
        HypothesisEpidemicVsNormal$NullHypLowerTailed[i]="Do Not Reject"
      }
    }

  #Hypothesis Testing - ANOVA
    #H0:No difference in population means in 5 zones
    #H1: There is a difference in population means in 5 zones
    res.aov1=aov(SampleDataAnova$no_cases~SampleDataAnova$Region)
    summary(res.aov1)
    TUKEY1=TukeyHSD(res.aov1)
    
    LABELS=generate_label_df(TUKEY1 , 'SampleDataAnova$Region')
    
    
    my_colors=c( rgb(143,199,74,maxColorValue = 255),rgb(242,104,34,maxColorValue = 255), rgb(111,145,202,maxColorValue = 255),rgb(254,188,18,maxColorValue = 255) , rgb(74,132,54,maxColorValue = 255),rgb(236,33,39,maxColorValue = 255),rgb(165,103,40,maxColorValue = 255))
    
    # Draw the basic boxplot
    a=boxplot(SampleDataAnova$no_cases ~ SampleDataAnova$Region , ylim=c(min(SampleDataAnova$no_cases) , 1.1*max(SampleDataAnova$no_cases)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="Weekly Cases" , main="")
    
    # To write the letter over each box
    over=0.1*max( a$stats[nrow(a$stats),] )
    
    #Add the labels
    text( c(1:nlevels(SampleDataAnova$Region)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )
    
    
  #Hypothesis Testing - ANOVA
    #H0:No difference in population means across the 6 years
    #H1: There is a difference in population mean across the 6 years
    PopulationData$Year = as.factor(PopulationData$Year)
    res.aov2=aov(PopulationData$no_cases~PopulationData$Year)
    summary(res.aov2)
    TUKEY = TukeyHSD(res.aov2,  'PopulationData$Year',conf.level = 0.95)
    summary(PopulationData$no_cases)
    
    LABELS=generate_label_df(TUKEY, 'PopulationData$Year')
    
    
    my_colors=c( rgb(143,199,74,maxColorValue = 255),rgb(242,104,34,maxColorValue = 255), rgb(111,145,202,maxColorValue = 255),rgb(254,188,18,maxColorValue = 255) , rgb(74,132,54,maxColorValue = 255),rgb(236,33,39,maxColorValue = 255),rgb(165,103,40,maxColorValue = 255))
    
    # Draw the basic boxplot
    a=boxplot(PopulationData$no_cases ~ PopulationData$Year , ylim=c(min(PopulationData$no_cases) , 1.1*max(PopulationData$no_cases)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="WeeklyCount" , main="")
    
    # I want to write the letter over each box. Over is how high I want to write it.
    over=0.1*max( a$stats[nrow(a$stats),] )
    
    #Add the labels
    text( c(1:nlevels(PopulationData$Year)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )
    
    
  #Hypothesis Testing - Difference in mean between localities near reservoirs and other localities within 1 zone
    #East
      #H0: Mu1-Mu0<=0
      #H1: Mu1-Mu0>0
      
    HypothesisReservoirsVsNormal=data.frame("Region"=c("East","West"))
    HypothesisReservoirsVsNormal$xbar=c(mean(SampleDataLocWBEast1$Cluster.Cases),mean(SampleDataLocWBWest1$Cluster.Cases))
    HypothesisReservoirsVsNormal$ybar=c(mean(SampleDataLocWBEast0$Cluster.Cases),mean(SampleDataLocWBWest0$Cluster.Cases))
    HypothesisReservoirsVsNormal$MeanDifferenceReservoirNormal=HypothesisReservoirsVsNormal$xbar-HypothesisReservoirsVsNormal$ybar
    sx=c(sd(SampleDataLocWBEast1$Cluster.Cases),sd(SampleDataLocWBWest1$Cluster.Cases))
    sy=c(sd(SampleDataLocWBEast0$Cluster.Cases),sd(SampleDataLocWBWest0$Cluster.Cases))
    nx=c(nrow(SampleDataLocWBEast1),nrow(SampleDataLocWBWest1))
    ny=c(nrow(SampleDataLocWBEast0),nrow(SampleDataLocWBWest0))
    alpha=0.05
    HypothesisReservoirsVsNormal$t=NA
    taTwoTailed=NA
    taUpperTailed=NA
    taLowerTailed=NA
    HypothesisReservoirsVsNormal$pvalueTwoTailed=NA
    HypothesisReservoirsVsNormal$pvalueUpperTailed=NA
    HypothesisReservoirsVsNormal$pvalueLowerTailed=NA
    HypothesisReservoirsVsNormal$NullHypTwoTailed=NA
    HypothesisReservoirsVsNormal$NullHypUpperTailed=NA
    HypothesisReservoirsVsNormal$NullHypLowerTailed=NA
    
    for(i in 1:2)
    {
      HypothesisReservoirsVsNormal$t[i]=(HypothesisReservoirsVsNormal$xbar[i]-HypothesisReservoirsVsNormal$ybar[i])/sqrt(sx[i]^2/nx[i]+sy[i]^2/ny[i])
      taTwoTailed[i]=qt(1-2*alpha/2,nx[i]+ny[i]-2)
      taUpperTailed[i]=qt(1-alpha,nx[i]+ny[i]-2)
      taLowerTailed[i]=qt(alpha,nx[i]+ny[i]-2)
      HypothesisReservoirsVsNormal$pvalueTwoTailed[i]= 2*(1-pt(HypothesisReservoirsVsNormal$t[i],nx[i]+ny[i]-2))
      HypothesisReservoirsVsNormal$pvalueUpperTailed[i]=(1-pt(HypothesisReservoirsVsNormal$t[i],nx[i]+ny[i]-2))
      HypothesisReservoirsVsNormal$pvalueLowerTailed[i]=pt(HypothesisReservoirsVsNormal$t[i],nx[i]+ny[i]-2)
      
      if(HypothesisReservoirsVsNormal$pvalueTwoTailed[i]<(2*alpha))
      {
        HypothesisReservoirsVsNormal$NullHypTwoTailed[i]="Reject"
      }else
      {
        HypothesisReservoirsVsNormal$NullHypTwoTailed[i]="Do Not Reject"
      }
      
      if(HypothesisReservoirsVsNormal$pvalueUpperTailed[i]<alpha)
      {
        HypothesisReservoirsVsNormal$NullHypUpperTailed[i]="Reject"
      }else
      {
        HypothesisReservoirsVsNormal$NullHypUpperTailed[i]="Do Not Reject"
      }
      if(HypothesisReservoirsVsNormal$pvalueLowerTailed[i]<alpha)
      {
        HypothesisReservoirsVsNormal$NullHypLowerTailed[i]="Reject"
      }else
      {
        HypothesisReservoirsVsNormal$NullHypLowerTailed[i]="Do Not Reject"
      }
    }

#6. VISUALISATION

    # #Line chart of temperature
    # ggplot(LRTable,aes(x=c(1:nrow(LRTable)), y=mean_temp))+
    #   geom_point()+
    #   geom_line()+
    #   labs(x="Month (1 to 72)")
    # 
    # #Line chart of number of days of rainfall
    # ggplot(LRTable,aes(x=c(1:nrow(LRTable)), y=no_of_rainy_days))+
    #   geom_point()+
    #   geom_line()+
    #   labs(x="Month (1 to 72)")
    # 
    # #Plotting Year vs Sum of Dengue cases, and temperature vs Number of dengue cases
    # ggplot(JoinedTable,aes(x=mean_temp, y= no_cases)) +
    # geom_point()
    # 
    # ggplot(JoinedTable,aes(x=no_of_rainy_days, y= no_cases)) +
    # geom_point()
    # 
    # #Plot Year wise number of cases
    # yearlyCases = JoinedTable %>% group_by(Year) %>% summarise(Yearly_No_Cases=sum(no_cases))
    # ggplot(yearlyCases,aes(x=Year, y= Yearly_No_Cases)) +
    # geom_point()+
    # geom_line()
    # 
    # #Plot Monthwise temperature vs number of cases
    # ggplot(LRTable,aes(x=mean_temp, y= no_cases)) +
    # geom_point()+
    # geom_smooth(method='lm', formula=y~x)
    # 
    # #Plot Monthwise rainfall vs number of cases
    # ggplot(LRTable,aes(x=no_of_rainy_days, y= no_cases)) +
    #   geom_point()+
    #   geom_smooth(method='lm', formula=y~x)
    # 
    # #Plot Number of cases against Week, grouped by Year
    # ggplot(JoinedTable,aes(x=Week, y=no_cases, group=Year, colour=Year))+
    #   geom_point()+geom_line()+
    #   scale_color_gradientn(colours=rainbow(6))
    # 
    # #Plot Number of cases against Week
    # ggplot(PopulationData,aes(x=c(1:nrow(PopulationData)), y=no_cases))+
    #   geom_point()+
    #   geom_line()+
    #   labs(x="Week(1 to 313)")
    # 
    # #Plotting histogram
    # hist(JoinedTable$no_cases, breaks=25)
    

#7. EXPORTING DATA    
# #Descriptive Stats
# setwd("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/Data/Resultant Tables from R/DescriptiveStats")
# write.csv(Desc_Stats_SampleData,"SampleDataOverall_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SampleDataEpidemic,"SampleDataEpidemic_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SampleDataNormal,"SampleDataNormal_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SgDengue,"MainDataOverall_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SgDengueEpidemic,"MainDataEpidemic_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SgDengueNormal,"MainDataNormal_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SgDengueFever,"MainData_DF_Zonewise_desc_stats.csv")
# write.csv(Desc_Stats_SgDengueHaemorrhagicFever,"MainData_DHF_Zonewise_desc_stats.csv")
# 
# #Inferential Stats
# setwd("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/Data/Resultant Tables from R/InferentialStats")
# write.csv(Mean,"OverallMean - Estimate-CI-HypTest.csv")
# write.csv(MeanNormal,"Normal Period Mean - Estimate-CI-HypTest.csv")
# write.csv(MeanEpidemic,"Epidemic Period Mean - Estimate-CI-HypTest.csv")
# write.csv(Proportion,"Overall Proportion - Estimate-CI-HypTest.csv")
# write.csv(ProportionNormal,"Normal Period Proportion - Estimate-CI-HypTest.csv")
# write.csv(ProportionEpidemic,"Epidemic Period Proportion - Estimate-CI-HypTest.csv")
# write.csv(OutbreakStats,"Outbreak Percentage Increase Statistics.csv")
# write.csv(HypothesisEpidemicVsNormal,"HypothesisEpidemicVsNormal.csv")
# write.csv(HypothesisReservoirsVsNormal,"HypothesisReservoirsVsNormal.csv")
# 
# #Wrangled Tables
# setwd("E:/MITB/Term 1/ISSS616 Applied Statistical Analysis with R/ASAR Project/Work/Data/Resultant Tables from R/WrangledTables")
# write.csv(JoinedTable,"JoinedTable.csv")
# write.csv(LRTable,"LRTable.csv")
# write.csv(PopulationData,"PopulationData.csv")
# write.csv(PopulationDataEpidemic,"PopulationDataEpidemic.csv")
# write.csv(PopulationDataNormal,"PopulationDataNormal.csv")
# write.csv(SampleData,"SampleData.csv")
# write.csv(SampleDataEpidemic,"SampleDataEpidemic.csv")
# write.csv(SampleDataNormal,"SampleDataNormal.csv")
# write.csv(SampleDataAnova,"SampleDataAnova.csv")
# write.csv(SampleDataLoc,"SampleDataLoc.csv")
# write.csv(SampleDataLocWB,"SampleDataLocWB.csv")
# write.csv(SampleDataLocWBEast,"SampleDataLocWBEast.csv")
# write.csv(SampleDataLocWBEast0,"SampleDataLocWBEast0.csv")
# write.csv(SampleDataLocWBEast1,"SampleDataLocWBEast1.csv")
# write.csv(SampleDataLocWBWest,"SampleDataLocWBWest.csv")
# write.csv(SampleDataLocWBWest0,"SampleDataLocWBWest0.csv")
# write.csv(SampleDataLocWBWest1,"SampleDataLocWBWest1.csv")
# write.csv(SgDengue,"SgDengue.csv")
# write.csv(SgDengueFever,"SgDengueFever.csv")
# write.csv(SgDengueHaemorrhagicFever,"SgDengueHaemorrhagicFever.csv")
# write.csv(SgRain1217,"SgRain1217.csv")
# write.csv(SgTemp1217,"SgTemp1217.csv")
# write.csv(SgVisitors,"SgVisitors.csv")
# write.csv(LRTablewithVisitors,"LRTablewithVisitors.csv")
# write.csv(ScaledLRTablewithVisitors,"ScaledLRTablewithVisitors.csv")
#     
# #SHINY - DESCRIPTIVE STATS
#     
#     
#     library(shiny)
#     library(ggplot2)
#     library(stringr)
#     
#     f = data.frame(Region = Desc_Stats_SampleData$Region,x = Desc_Stats_SampleData$Min)
#     class(f$Region)
#     class(f$x)
#     
#     plot<-ggplot(data = LRTable,aes(x=as.factor(Month),y= mean_temp, group = Year, colour=Year))
#     three<-plot+geom_point(stat = "identity")+geom_line(stat="identity")+ggtitle("Temperature Vs Month")
#     print(three)
#     barplot(f$x,names.arg = f$Region)
#     
#     plot(f$Region,f$x)
#     ui <- fluidPage(
#       
#       titlePanel("Dengue Analysis"),
#       tabsetPanel(
#         tabPanel( "Stats",
#                   sidebarLayout(
#                     sidebarPanel (
#                       selectInput("Data","Data:",c("Sample Data"= "sampledata", "Population Data" = "populationdata")),
#                       selectInput("Period","Period::",c("Total"= "Total", "Outbreak" = "Outbreak","Normal" = "Normal")),
#                       selectInput("Measure","Measure::",c("Mean"= "Mean", "Median" = "Median","Mode" = "Mode"))
#                       
#                       
#                     )
#                     ,
#                     mainPanel(
#                       plotOutput(outputId = "barplot"))))))
#     
#     server <- function(input, output, session) {
#       
#       
#       
#       output$barplot = renderPlot({ 
#         
#         if(input$Data == "sampledata")
#         {
#           
#           if(input$Period == "Total")
#           {
#             if(input$Measure == "Mean")
#             {
#               a = barplot(Desc_Stats_SampleData$Mean[2:6],names.arg = Desc_Stats_SampleData$Region[2:6])
#             }
#             
#             if(input$Measure == "Median")
#             {
#               a = barplot(Desc_Stats_SampleData$Median[2:6],names.arg = Desc_Stats_SampleData$Region[2:6] ) 
#             }
#             
#             if(input$Measure == "Mode")
#             {
#               a = barplot(Desc_Stats_SampleData$Mode[2:6],names.arg = Desc_Stats_SampleData$Region[2:6] )
#             }
#             
#           } else if(input$Period == "Outbreak")
#           {
#             if(input$Measure == "Mean")
#             {
#               a = barplot(Desc_Stats_SampleDataEpidemic$Mean[2:6],names.arg = Desc_Stats_SampleDataEpidemic$Region[2:6])
#             } 
#             
#             if(input$Measure == "Median")
#             {
#               a = barplot(Desc_Stats_SampleDataEpidemic$Median[2:6],names.arg = Desc_Stats_SampleDataEpidemic$Region[2:6] ) 
#             }
#             
#             if(input$Measure == "Mode")
#             {
#               a = barplot(Desc_Stats_SampleDataEpidemic$Mode[2:6],names.arg = Desc_Stats_SampleDataEpidemic$Region[2:6] )
#             }
#           } else {
#             
#             if(input$Measure == "Mean")
#             {
#               a = barplot(Desc_Stats_SampleDataNormal$Mean[2:6],names.arg = Desc_Stats_SampleDataNormal$Region[2:6])
#             }
#             
#             if(input$Measure == "Median")
#             {
#               a = barplot(Desc_Stats_SampleDataNormal$Median[2:6],names.arg = Desc_Stats_SampleDataNormal$Region[2:6] ) 
#             }
#             
#             if(input$Measure == "Mode")
#             {
#               a = barplot(Desc_Stats_SampleDataNormal$Mode[2:6],names.arg = Desc_Stats_SampleDataNormal$Region[2:6] )
#             }
#           }
#         }else 
#         {
#           if(input$Period == "Total")
#           {
#             if(input$Measure == "Mean")
#             {
#               a = barplot(Desc_Stats_SgDengue$Mean[2:6],names.arg = Desc_Stats_SgDengue$Region[2:6])
#             }
#             
#             if(input$Measure == "Median")
#             {
#               a = barplot(Desc_Stats_SgDengue$Median[2:6],names.arg = Desc_Stats_SgDengue$Region[2:6] ) 
#             }
#             
#             if(input$Measure == "Mode")
#             {
#               a = barplot(Desc_Stats_SgDengue$Mode[2:6],names.arg = Desc_Stats_SgDengue$Region[2:6] )
#             }
#             
#           } else if(input$Period == "Outbreak")
#           {
#             if(input$Measure == "Mean")
#             {
#               a = barplot(Desc_Stats_SgDengueEpidemic$Mean[2:6],names.arg = Desc_Stats_SgDengueEpidemic$Region[2:6])
#             } 
#             
#             if(input$Measure == "Median")
#             {
#               a = barplot(Desc_Stats_SgDengueEpidemic$Median[2:6],names.arg = Desc_Stats_SgDengueEpidemic$Region[2:6] ) 
#             }
#             
#             if(input$Measure == "Mode")
#             {
#               a = barplot(Desc_Stats_SgDengueEpidemic$Mode[2:6],names.arg = Desc_Stats_SgDengueEpidemic$Region[2:6] )
#             }
#           } else {
#             
#             if(input$Measure == "Mean")
#             {
#               a = barplot(Desc_Stats_SgDengueNormal$Mean[2:6],names.arg = Desc_Stats_SgDengueNormal$Region[2:6])
#             }
#             
#             if(input$Measure == "Median")
#             {
#               a = barplot(Desc_Stats_SgDengueNormal$Median[2:6],names.arg = Desc_Stats_SgDengueNormal$Region[2:6] ) 
#             }
#             
#             if(input$Measure == "Mode")
#             {
#               a = barplot(Desc_Stats_SgDengueNormal$Mode[2:6],names.arg = Desc_Stats_SgDengueNormal$Region[2:6] )
#             }
#           }
#         }
#         
#         
#         print(a)
#         
#       })
#     }
#     shinyApp(ui, server)
    
    
    
    
    
    
    