# Script used to visualize and analyze the Palau Temperature Data from CRRF and Cohen Lab Loggers
## Loading libraries ####
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(xts)
library(zoo)
library(TTR)
library(scales)
library(ggpubr)
library(signal)
library(ggridges)

setwd("/Users/hannyrivera/Documents/MIT/Research_Papers/Palau_Microsat_Paper/Drafts/Comms_Bio_submission/Palau_porites/Site_Temperatures/")

## Importing data for weekly temps ####
# Import data from R friendly (tidy) version (data reformatted in excel from HOBO logger output) 
# Using timezone = GMT make it so that it just reads the date time values as is and doesn't try to convert it as these were logggin in Palau local time
DO<-xts(zoo(read.table("DropOff_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2000,3,2,0,0,0),ISOdate(2017,8,16,23,30,0), "30 min", tz="GMT")))
Ngerdiluches<-xts(zoo(read.table("Ngerdiluches_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2000,6,7,0,0,0),ISOdate(2017,6,8,23,30,0), "30 min", tz="GMT")))
Helen<-xts(zoo(read.table("Helen_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2014,8,30,0,0,0),ISOdate(2017,2,15,23,30,0), "30 min", tz="GMT")))   
Kayan<-xts(zoo(read.table("Kayangel_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2010,4,1,0,0,0),ISOdate(2017,8,16,23,30,0), "30 min", tz="GMT")))
Mecherchar<-xts(zoo(read.table("Mecherchar_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2011,4,25,0,0,0),ISOdate(2017,6,12,23,30,0), "30 min", tz="GMT")))
Ngerchelong<-xts(zoo(read.table("Ngerchelong_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2000,6,27,0,0,0),ISOdate(2017,6,29,23,30,0), "30 min", tz="GMT")))
Ngermid<-xts(zoo(read.table("Ngermid_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2007,5,25,0,0,0),ISOdate(2017,6,27,23,30,0), "30 min", tz="GMT")))
Risong<-xts(zoo(read.table("Risong_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2011,9,19,0,0,0),ISOdate(2013,5,4,23,45,0), "15 min", tz="GMT")))
Taoch<-xts(zoo(read.table("Taoch_tidy.txt", header=TRUE)$Temperature,seq.POSIXt(ISOdate(2012,3,28,0,0,0),ISOdate(2012,12,5,23,45,0), "15 min", tz="GMT")))

## Weekly average temps analysis/plot (Fig 1A) ####
# mean raw temps for outer reefs
mean(c(mean(as.numeric(DO), na.rm=TRUE),
       mean(as.numeric(Ngerdiluches), na.rm=TRUE),
       mean(as.numeric(Helen), na.rm=TRUE),
       mean(as.numeric(Kayan), na.rm=TRUE),
       mean(as.numeric(Ngerchelong), na.rm=TRUE)))
# 29.11208

#average standard dev
mean(c(sd(as.numeric(DO), na.rm=TRUE),
       sd(as.numeric(Ngerdiluches), na.rm=TRUE),
       sd(as.numeric(Helen), na.rm=TRUE),
       sd(as.numeric(Kayan), na.rm=TRUE),
       sd(as.numeric(Ngerchelong), na.rm=TRUE)))
#0.6954924

# mean raw temps for rock islands
mean(c(mean(as.numeric(Mecherchar), na.rm=TRUE),
       mean(as.numeric(Ngermid), na.rm=TRUE),
       mean(as.numeric(Risong), na.rm=TRUE),
       mean(as.numeric(Taoch), na.rm=TRUE)))
#30.28624

#average standard dev
mean(c(sd(as.numeric(Mecherchar), na.rm=TRUE),
       sd(as.numeric(Ngermid), na.rm=TRUE),
       sd(as.numeric(Risong), na.rm=TRUE),
       sd(as.numeric(Taoch), na.rm=TRUE)))
# 0.6191633

#Calculate yearly top 90 percentile temperatures 
#Outer reefs
DO_perct<-apply.yearly(DO,  FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))
Ngerdiluches_perct<-apply.yearly(Ngerdiluches, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))  
Helen_perct<-apply.yearly(Helen, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))
Kayan_perct<-apply.yearly(Kayan, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))
Ngerchelong_perct<-apply.yearly(Ngerchelong, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))

mean(c(mean(as.numeric(DO_perct), na.rm=TRUE),
       mean(as.numeric(Ngerdiluches_perct), na.rm=TRUE),
       mean(as.numeric(Helen_perct), na.rm=TRUE),
       mean(as.numeric(Kayan_perct), na.rm=TRUE),
       mean(as.numeric(Ngerchelong_perct), na.rm=TRUE)))
#  29.78032

#Rock Islands
Mecherchar_perct<-apply.yearly(Mecherchar, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))
Ngermid_perct<-apply.yearly(Ngermid, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))
Taoch_perct<-apply.yearly(Taoch, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))
Risong_perct<-apply.yearly(Risong, FUN=function(x) quantile(x,c(0.9), na.rm=TRUE))

mean(c(mean(as.numeric(Mecherchar_perct), na.rm=TRUE),
       mean(as.numeric(Ngermid_perct), na.rm=TRUE),
       mean(as.numeric(Risong_perct), na.rm=TRUE),
       mean(as.numeric(Taoch_perct), na.rm=TRUE)))
# 31.08199

# Merge all the data together taking weekly averages for plotting 
AllData<-merge(apply.weekly(DO["2010/"],FUN = mean), apply.weekly(Ngerdiluches["2010/"], FUN=mean), join="left", fill=NA)
AllData<-merge(AllData, apply.weekly(Helen["2010/"],FUN=mean), join="left", fill=NA)
AllData<-merge(AllData, apply.weekly(Ngerchelong["2010/"],FUN=mean), join="left", fill=NA)
AllData<-merge(AllData, apply.weekly(Kayan["2010/"],FUN=mean), join="left", fill=NA)
AllData<-merge(AllData, apply.weekly(Mecherchar["2010/"],FUN=mean), join="left", fill=NA)
AllData<-merge(AllData, apply.weekly(Ngermid["2010/"],FUN=mean), join="left", fill=NA)
## For the Cohen loggers, need to combine indexes to be on same scale
Risongweekly<-apply.weekly(Risong["2011/"], FUN=mean)
Taochweekly<-apply.weekly(Taoch["2012/"], FUN=mean)
index(Risongweekly)<-index(Risongweekly)-15*60
index(Taochweekly)<-index(Taochweekly)-15*60

AllData<-merge(AllData, Taochweekly, join="left", fill=NA)
AllData<-merge(AllData, Risongweekly, join="left", fill=NA)

# Name columns
colnames(AllData)<-c("Drop_Off", "Ngerdiluches", "Helen", "Ngerchelong", "Kayangel","Mecherchar", "Ngermid", "Taoch", "Risong")
# Convert to data frame for ggplot plotting
AllData_gg<-as.data.frame(AllData)
#Readjust columns to the order in which they appear in the rest of the paper (Rock Islands first)
AllData_gg<-AllData_gg[, c("Mecherchar","Taoch", "Risong" ,"Ngermid", "Helen", "Drop_Off", "Ngerdiluches", "Ngerchelong", "Kayangel")]
# renames rows as increasing numbers
row.names(AllData_gg)<-as.character(seq(1,length(AllData_gg[,1])))
# add in date time values as variable
AllData_gg<-cbind(time(apply.weekly(DO["2010/"],FUN = mean)),AllData_gg)
# fix column header 
colnames(AllData_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
AllData_gg<-melt(AllData_gg, id.vars = "DateTime", measure.vars = c("Mecherchar", "Taoch", "Risong", "Ngermid", "Helen", "Drop_Off", "Ngerdiluches", "Ngerchelong","Kayangel"))
# fix column headers again
colnames(AllData_gg)<-c("DateTime", "Site", "Temp")

# Plot weekly time series
(A<-ggplot(AllData_gg,aes(x=DateTime, y=Temp, colour=Site))+geom_line(size=.6)+theme_minimal()+
  scale_y_continuous(breaks= c(28, 30, 32), labels=c("28°C", "30°C","32°C"),limits = c(26.8,32.35))+
  scale_x_datetime(date_labels="%Y", date_breaks="year", date_minor_breaks = "months", expand=c(0.02,0.0))+ggtitle("")+
  ylab("")+xlab("")+
    theme(axis.text.x = element_text(size=10, face="bold"),
          plot.title = element_text(hjust=-0.07, size = 7),
          axis.text.y = element_text(size=10, face="bold"),
          panel.grid.major.y = element_line(size=.75),
          panel.grid.major.x = element_line(size=.75), 
          panel.grid.minor.x = element_line(size=.3), 
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom", 
          legend.title = element_blank(), 
          legend.text = element_text(size=10, face="bold"),
          plot.margin = margin(15,15,0,0,"pt"))+
  scale_color_manual(values=c("coral2", "darkgoldenrod3", "brown4", "maroon3","springgreen1", "skyblue2","turquoise3", "springgreen4", "royalblue3"),
                     labels=c("Mecherchar","Taoch", "Risong", "Ngermid", "Helen", "Drop Off", "Ngerdiluches", "Ngerchelong", "Kayangel"))+
  guides(colour=FALSE)+theme(plot.margin = margin(0,0,0,0,"cm"))+
  guides(color = guide_legend(overrride.aes = list (size=10))))

  
## Import data for daily variability ####
# For this, bring in the data that was frequency filtered in MATLAB to remove high and low frequency variation (see temp_filtering.m)
setwd("Freq_filtered/")
DO_filt<-read.table("DO_filtered_temp.txt")
Ngerd_filt<-read.table("Ngerdiluches_filtered_temp.txt")
Ngerchelong_filt<-read.table("Ngerchelong_filtered_temp.txt")
Ngermid_filt<-read.table("Ngermid_filtered_temp.txt")
Mecherchar_filt<-read.table("Mecherchar_filtered_temp.txt")
Risong_filt<-read.table("Risong_filtered_temp.txt")
Taoch_filt<-read.table("Toach_filtered_temp.txt")
Helen_filt<-read.table("Helen_filtered_temp.txt")
Kayan_filt<-read.table("Kayangel_filtered_temp.txt")

colnames(DO_filt)<-c("Temp")
colnames(Ngerd_filt)<-c("Temp")
colnames(Ngerchelong_filt)<-c("Temp")
colnames(Ngermid_filt)<-c("Temp")
colnames(Mecherchar_filt)<-c("Temp")
colnames(Risong_filt)<-c("Temp")
colnames(Taoch_filt)<-c("Temp")
colnames(Helen_filt)<-c("Temp")
colnames(Kayan_filt)<-c("Temp")


### 
# Adding the mean back to the values 
# The mean values were calculated from the time series in MATLAB since that's were I did the filtering 
# I have to add back the mean because MATLAB outputs the filtered data as an anomaly from the mean so the scale wouldn't make as much sense
# in the context of the other plots 

DO_filt$Temp_adj<-DO_filt$Temp+29.1878
Helen_filt$Temp_adj<-Helen_filt$Temp+28.9612
Kayan_filt$Temp_adj<-Kayan_filt$Temp+29.0864
Mecherchar_filt$Temp_adj<-Mecherchar_filt$Temp+30.7557
Ngerd_filt$Temp_adj<-Ngerd_filt$Temp+29.3951
Ngerchelong_filt$Temp_adj<-Ngerchelong_filt$Temp+28.9298
Ngermid_filt$Temp_adj<-Ngermid_filt$Temp+30.2928
Risong_filt$Temp_adj<-Risong_filt$Temp+30.3357
Taoch_filt$Temp_adj<-Taoch_filt$Temp+29.7604

#Add in the sampling intervals
DO_filt$Sample<-rep(seq(1,48), length(DO_filt$Temp)/48)
Helen_filt$Sample<-rep(seq(1,48), length(Helen_filt$Temp)/48)
Kayan_filt$Sample<-rep(seq(1,48), length(Kayan_filt$Temp)/48)
Mecherchar_filt$Sample<-rep(seq(1,48), length(Mecherchar_filt$Temp)/48)
Ngerd_filt$Sample<-rep(seq(1,48), length(Ngerd_filt$Temp)/48)
Ngerchelong_filt$Sample<-rep(seq(1,48), length(Ngerchelong_filt$Temp)/48)
#Risong, Ngermid, and Taoch are at 15 min sampling intervals
Risong_filt$Day<-rep(seq(1,96), length(Risong_filt$Temp)/96)
Taoch_filt$Day<-rep(seq(1,96), length(Taoch_filt$Temp)/96)
Ngermid_filt$Day<-rep(seq(1,96), length(Ngermid_filt$Temp)/96)

## Daily variability analyses and plot (Fig 1B) ####
# Calculate average daily range at each site and plot density plots of ranges
DO_var<-setNames(data.frame(matrix(ncol=5, nrow = length(DO_filt$Temp)/48)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
DO_var$Site<-"Drop_Off"
Helen_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Helen_filt$Temp)/48)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
Helen_var$Site<-"Helen"
Kayan_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Kayan_filt$Temp)/48)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
Kayan_var$Site<-"Kayangel"
Mecherchar_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Mecherchar_filt$Temp)/48)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
Mecherchar_var$Site<-"Mecherchar"
Ngerd_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Ngerd_filt$Temp)/48)), c("MaxTemp", "MinTemp", "Range","Variance",  "Site"))
Ngerd_var$Site<-"Ngerdiluches"
Ngermid_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Ngermid_filt$Temp)/96)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
Ngermid_var$Site<-"Ngermid"
Ngerchelong_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Ngerchelong_filt$Temp)/48)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
Ngerchelong_var$Site<-"Ngerchelong"
Risong_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Risong_filt$Temp)/96)), c("MaxTemp", "MinTemp", "Range","Variance",  "Site"))
Risong_var$Site<-"Risong"
Taoch_var<-setNames(data.frame(matrix(ncol=5, nrow = length(Taoch_filt$Temp)/96)), c("MaxTemp", "MinTemp", "Range", "Variance", "Site"))
Taoch_var$Site<-"Taoch"

i=1
for (x in seq(1,length(DO_filt$Temp_adj)/48)){
  DO_var$MaxTemp[x]<-max(DO_filt$Temp_adj[i:(i+47)])
  DO_var$MinTemp[x]<-min(DO_filt$Temp_adj[i:(i+47)])
  DO_var$Range[x]<-max(DO_filt$Temp_adj[i:(i+47)])-min(DO_filt$Temp_adj[i:(i+47)])
  DO_var$Variance[x]<-var(DO_filt$Temp_adj[i:(i+47)])
  i=i+48
}

i=1
for (x in seq(1,length(Helen_filt$Temp_adj)/48)){
  Helen_var$MaxTemp[x]<-max(Helen_filt$Temp_adj[i:(i+47)])
  Helen_var$MinTemp[x]<-min(Helen_filt$Temp_adj[i:(i+47)])
  Helen_var$Range[x]<-max(Helen_filt$Temp_adj[i:(i+47)])-min(Helen_filt$Temp_adj[i:(i+47)])
  Helen_var$Variance[x]<-var(Helen_filt$Temp_adj[i:(i+47)])
  i=i+48
}


i=1
for (x in seq(1,(length(Kayan_filt$Temp_adj)/48)-1)){
  Kayan_var$MaxTemp[x]<-max(Kayan_filt$Temp_adj[i:(i+47)])
  Kayan_var$MinTemp[x]<-min(Kayan_filt$Temp_adj[i:(i+47)])
  Kayan_var$Range[x]<-max(Kayan_filt$Temp_adj[i:(i+47)])-min(Kayan_filt$Temp_adj[i:(i+47)])
  Kayan_var$Variance[x]<-var(Kayan_filt$Temp_adj[i:(i+47)])
  i=i+48
}


i=1
for (x in seq(1,length(Mecherchar_filt$Temp_adj)/48)){
  Mecherchar_var$MaxTemp[x]<-max(Mecherchar_filt$Temp_adj[i:(i+47)])
  Mecherchar_var$MinTemp[x]<-min(Mecherchar_filt$Temp_adj[i:(i+47)])
  Mecherchar_var$Range[x]<-max(Mecherchar_filt$Temp_adj[i:(i+47)])-min(Mecherchar_filt$Temp_adj[i:(i+47)])
  Mecherchar_var$Variance[x]<-var(Mecherchar_filt$Temp_adj[i:(i+47)])
  i=i+48
}

i=1
for (x in seq(1,length(Ngerd_filt$Temp_adj)/48)){
  Ngerd_var$MaxTemp[x]<-max(Ngerd_filt$Temp_adj[i:(i+47)])
  Ngerd_var$MinTemp[x]<-min(Ngerd_filt$Temp_adj[i:(i+47)])
  Ngerd_var$Range[x]<-max(Ngerd_filt$Temp_adj[i:(i+47)])-min(Ngerd_filt$Temp_adj[i:(i+47)])
  Ngerd_var$Variance[x]<-var(Ngerd_filt$Temp_adj[i:(i+47)])
  i=i+48
}

i=1
for (x in seq(1,length(Ngermid_filt$Temp_adj)/96)){
  Ngermid_var$MaxTemp[x]<-max(Ngermid_filt$Temp_adj[i:(i+95)])
  Ngermid_var$MinTemp[x]<-min(Ngermid_filt$Temp_adj[i:(i+95)])
  Ngermid_var$Range[x]<-max(Ngermid_filt$Temp_adj[i:(i+95)])-min(Ngermid_filt$Temp_adj[i:(i+95)])
  Ngermid_var$Variance[x]<-var(Ngermid_filt$Temp_adj[i:(i+95)])
  i=i+96
}

i=1
for (x in seq(1,length(Ngerchelong_filt$Temp_adj)/48)){
  Ngerchelong_var$MaxTemp[x]<-max(Ngerchelong_filt$Temp_adj[i:(i+47)])
  Ngerchelong_var$MinTemp[x]<-min(Ngerchelong_filt$Temp_adj[i:(i+47)])
  Ngerchelong_var$Range[x]<-max(Ngerchelong_filt$Temp_adj[i:(i+47)])-min(Ngerchelong_filt$Temp_adj[i:(i+47)])
  Ngerchelong_var$Variance[x]<-var(Ngerchelong_filt$Temp_adj[i:(i+47)])
  i=i+48
}

i=1
for (x in seq(1,length(Risong_filt$Temp_adj)/96)){
  Risong_var$MaxTemp[x]<-max(Risong_filt$Temp_adj[i:(i+95)])
  Risong_var$MinTemp[x]<-min(Risong_filt$Temp_adj[i:(i+95)])
  Risong_var$Range[x]<-max(Risong_filt$Temp_adj[i:(i+95)])-min(Risong_filt$Temp_adj[i:(i+95)])
  Risong_var$Variance[x]<-var(Risong_filt$Temp_adj[i:(i+95)])
  i=i+96
}

i=1
for (x in seq(1,length(Taoch_filt$Temp_adj)/96)){
  Taoch_var$MaxTemp[x]<-max(Taoch_filt$Temp_adj[i:(i+95)])
  Taoch_var$MinTemp[x]<-min(Taoch_filt$Temp_adj[i:(i+95)])
  Taoch_var$Range[x]<-max(Taoch_filt$Temp_adj[i:(i+95)])-min(Taoch_filt$Temp_adj[i:(i+95)])
  Taoch_var$Variance[x]<-var(Taoch_filt$Temp_adj[i:(i+95)])
  i=i+96
}

DO_var$Region="OR"
Helen_var$Region="OR"
Kayan_var$Region="OR"
Ngerchelong_var$Region="OR"
Ngerd_var$Region="OR"
Mecherchar_var$Region="RI"
Ngermid_var$Region="RI"
Risong_var$Region="RI"
Taoch_var$Region="RI"

All_var<-rbind(DO_var, Helen_var, Kayan_var, Mecherchar_var,Ngerd_var, Ngermid_var, Ngerchelong_var, Risong_var, Taoch_var)

OR_var<-subset(All_var, Site!="Mecherchar" & Site!="Taoch" & Site!="Risong" & Site!="Ngermid")
mean(OR_var$Range, na.rm=TRUE)
# 0.3167784
RI_var<-subset(All_var, Site=="Mecherchar" | Site=="Taoch" | Site=="Risong" | Site=="Ngermid")
mean(RI_var$Range, na.rm=TRUE)
# 0.5925325


## ggridges variability plot 
(B<-ggplot(All_var, aes(x=Range, y=Site, fill=Site))+geom_density_ridges(alpha=0.7, quantile_lines=TRUE, quantiles=2)+
  theme_bw()+guides(fill="none")+ xlab("Diurnal Range (°C)")+ylab("")+
  scale_x_continuous(limits=c(-0.03,2))+
  scale_fill_manual(values=c("coral2", "darkgoldenrod3", "brown4", "maroon3","springgreen1", "skyblue2","turquoise3", "springgreen4", "royalblue3"), 
                     limits=c("Mecherchar","Taoch", "Risong","Ngermid",  "Helen", "Drop_Off", "Ngerdiluches", "Ngerchelong", "Kayangel"),
                     labels=c("Mecherchar","Taoch" , "Risong","Ngermid",  "Helen", "Drop Off", "Ngerdiluches", "Ngerchelong", "Kayangel"))+
  scale_y_discrete(limits=c("Mecherchar", "Risong","Ngermid","Taoch", "Ngerchelong", "Helen", "Drop_Off", "Kayangel", "Ngerdiluches"), 
                   labels=c("Mecherchar", "Risong","Ngermid","Taoch", "Ngerchelong", "Helen", "Drop Off", "Kayangel", "Ngerdiluches"))+
  theme(axis.text = element_text(face="bold"),
        axis.title =element_text(face="bold")))

## Final Figure 1 #### 
ggarrange(A, B, common.legend = TRUE, ncol = 2, nrow = 1, labels=c("A","B"), legend = "bottom",widths=c(2,1))
ggsave("../figure1.png", units="in", height = 4, width=8, dpi=300)

          