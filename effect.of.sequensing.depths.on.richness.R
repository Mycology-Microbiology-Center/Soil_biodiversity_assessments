#Richness
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
library(MuMIn)
library("parameters")
#bacteria
library(scales)
library(effectsize)
bacteria.normal.m <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal.m$method<-gsub("40A","_40A",bacteria.normal.m$method)
bacteria.normal.m$method<-gsub("40B","_40B",bacteria.normal.m$method)
bacteria.normal.m$method[bacteria.normal.m$method=="GSMc"]<-"GSMc_62"
##bacteria soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.soil.m<-read.csv("bacteria.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.soil.m<-bacteria.soil.m[grepl("soil",row.names(bacteria.soil.m)),]
row.names(bacteria.soil.m)<-gsub("soil","",row.names(bacteria.soil.m))
bacteria.soil.m$site<-row.names(bacteria.soil.m)
bacteria.soil.m<-merge(bacteria.soil.m,soil.metadata,by="site",all=T)
bacteria.soil.m$site2<-substr(bacteria.soil.m$site,1,2)
bacteria.soil.m$estimate<-round(bacteria.soil.m$estimate)
##visualize
bacteria.soilLV<-bacteria.soil.m[grepl("LV",bacteria.soil.m$site2),]
bacteria.soilLV<-bacteria.soilLV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="mean"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,2,3)]
bacteria.soil<-merge(bacteria.soilLV,bacteria.normalLV,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LV"
bacteria.soil$pooling.effect<-round(bacteria.soil$estimate)/round(bacteria.soil$normalestimate)
bacteria.soilLV<-bacteria.soil
##LW
bacteria.soilLW<-bacteria.soil.m[grepl("LW",bacteria.soil.m$site2),]
bacteria.soilLW<-bacteria.soilLW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="mean"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,2,3)]
bacteria.soil<-merge(bacteria.soilLW,bacteria.normalLW,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LW"
bacteria.soil$pooling.effect<-round((bacteria.soil$estimate))/round(bacteria.soil$normalestimate)
bacteria.soilLW<-bacteria.soil
##LZ
bacteria.soilLZ<-bacteria.soil.m[grepl("LZ",bacteria.soil.m$site2),]
bacteria.soilLZ<-bacteria.soilLZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="mean"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,2,3)]
#
bacteria.soil<-merge(bacteria.soilLZ,bacteria.normalLZ,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$pooling.effect<-(round(bacteria.soil$estimate))/round(bacteria.soil$normalestimate)
bacteria.soilLZ<-bacteria.soil
bacteria.all<-rbind(bacteria.soilLZ,bacteria.soilLV,bacteria.soilLW)
##bacteria DNA
bacteria.normal.m <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal.m$method<-gsub("40A","_40A",bacteria.normal.m$method)
bacteria.normal.m$method<-gsub("40B","_40B",bacteria.normal.m$method)
bacteria.normal.m$method[bacteria.normal.m$method=="GSMc"]<-"GSMc_62"
##bacteria DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.DNA.m<-read.csv("bacteria.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.DNA.m<-bacteria.DNA.m[grepl("DNA",row.names(bacteria.DNA.m)),]
row.names(bacteria.DNA.m)<-gsub("DNA","",row.names(bacteria.DNA.m))
bacteria.DNA.m$site<-row.names(bacteria.DNA.m)

bacteria.DNA.m<-merge(bacteria.DNA.m,DNA.metadata,by="site",all=T)
bacteria.DNA.m$site2<-substr(bacteria.DNA.m$site,1,2)
bacteria.DNA.m$estimate<-round(bacteria.DNA.m$estimate)
##visualize
bacteria.DNALV<-bacteria.DNA.m[grepl("LV",bacteria.DNA.m$site2),]
bacteria.DNALV<-bacteria.DNALV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="mean"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,2,3)]
bacteria.DNA<-merge(bacteria.DNALV,bacteria.normalLV,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LV"
bacteria.DNA$pooling.effect<-round(bacteria.DNA$estimate)/round(bacteria.DNA$normalestimate)
bacteria.DNALV<-bacteria.DNA
##LW
bacteria.DNALW<-bacteria.DNA.m[grepl("LW",bacteria.DNA.m$site2),]
bacteria.DNALW<-bacteria.DNALW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="mean"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,2,3)]
bacteria.DNA<-merge(bacteria.DNALW,bacteria.normalLW,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LW"
bacteria.DNA$pooling.effect<-round((bacteria.DNA$estimate))/round(bacteria.DNA$normalestimate)
bacteria.DNALW<-bacteria.DNA
##LZ
bacteria.DNALZ<-bacteria.DNA.m[grepl("LZ",bacteria.DNA.m$site2),]
bacteria.DNALZ<-bacteria.DNALZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="mean"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,2,3)]
#
bacteria.DNA<-merge(bacteria.DNALZ,bacteria.normalLZ,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$pooling.effect<-(round(bacteria.DNA$estimate))/round(bacteria.DNA$normalestimate)
bacteria.DNALZ<-bacteria.DNA
bacteria.all2<-rbind(bacteria.DNALZ,bacteria.DNALV,bacteria.DNALW)
##separate the site
bacteria.all$pooling<-"soil"
bacteria.all2$pooling<-"DNA"
all<-rbind(bacteria.all,bacteria.all2)
all<-all[,c(1,4,5:7)]
all2<-aggregate(all[,4],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
out<-ggplot(all2, aes(x = method, y = x,fill=pooling)) +
  geom_bar(stat = "identity",position = 'dodge') +
  facet_wrap(~site2)+
  scale_x_discrete(limits = as.character(c( "GSMc_62", "GSMc_40A", "GSMc_40B","Zobel","DarkDiv", "LUCAS", "SUCC", "deep_SUCC", "deep"))) +
  labs(x = "Method", y = "Ratio of pooled methods and separated samples' richness") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )+
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey")

cbPalette <- rev(c("grey","#55AAFF"))
out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
bacteria.pool.ratio<-all2
bacteria.pool.ratio$community<-"bacteria"
###
#bacteria
#divide each method
bacteria.normal.m <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal.m$method<-gsub("40A","_40A",bacteria.normal.m$method)
bacteria.normal.m$method<-gsub("40B","_40B",bacteria.normal.m$method)
bacteria.normal.m$method[bacteria.normal.m$method=="GSMc"]<-"GSMc_62"
##bacteria soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.soil.m<-read.csv("bacteria.richness.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.soil.m<-bacteria.soil.m[grepl("soil",row.names(bacteria.soil.m)),]
row.names(bacteria.soil.m)<-gsub("soil","",row.names(bacteria.soil.m))
bacteria.soil.m$site<-row.names(bacteria.soil.m)
bacteria.soil.m<-merge(bacteria.soil.m,soil.metadata,by="site",all=T)
bacteria.soil.m$site2<-substr(bacteria.soil.m$site,1,2)
bacteria.soil.m$estimate<-round(bacteria.soil.m$estimate)
##visualize
bacteria.soilLV<-bacteria.soil.m[grepl("LV",bacteria.soil.m$site2),]
bacteria.soilLV<-bacteria.soilLV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="mean"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,2,3)]
bacteria.soil<-merge(bacteria.soilLV,bacteria.normalLV,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LV"
bacteria.soil$pooling.effect<-round(bacteria.soil$estimate)/round(bacteria.soil$normalestimate)
bacteria.soilLV<-bacteria.soil
##LW
bacteria.soilLW<-bacteria.soil.m[grepl("LW",bacteria.soil.m$site2),]
bacteria.soilLW<-bacteria.soilLW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="mean"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,2,3)]
bacteria.soil<-merge(bacteria.soilLW,bacteria.normalLW,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LW"
bacteria.soil$pooling.effect<-round((bacteria.soil$estimate))/round(bacteria.soil$normalestimate)
bacteria.soilLW<-bacteria.soil
##LZ
bacteria.soilLZ<-bacteria.soil.m[grepl("LZ",bacteria.soil.m$site2),]
bacteria.soilLZ<-bacteria.soilLZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="mean"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,2,3)]
#
bacteria.soil<-merge(bacteria.soilLZ,bacteria.normalLZ,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$pooling.effect<-(round(bacteria.soil$estimate))/round(bacteria.soil$normalestimate)
bacteria.soilLZ<-bacteria.soil
bacteria.all<-rbind(bacteria.soilLZ,bacteria.soilLV,bacteria.soilLW)
##bacteria DNA
bacteria.normal.m <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal.m$method<-gsub("40A","_40A",bacteria.normal.m$method)
bacteria.normal.m$method<-gsub("40B","_40B",bacteria.normal.m$method)
bacteria.normal.m$method[bacteria.normal.m$method=="GSMc"]<-"GSMc_62"
##bacteria DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.DNA.m<-read.csv("bacteria.richness.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.DNA.m<-bacteria.DNA.m[grepl("DNA",row.names(bacteria.DNA.m)),]
row.names(bacteria.DNA.m)<-gsub("DNA","",row.names(bacteria.DNA.m))
bacteria.DNA.m$site<-row.names(bacteria.DNA.m)
bacteria.DNA.m<-merge(bacteria.DNA.m,DNA.metadata,by="site",all=T)
bacteria.DNA.m$site2<-substr(bacteria.DNA.m$site,1,2)
bacteria.DNA.m$estimate<-round(bacteria.DNA.m$estimate)
##visualize
bacteria.DNALV<-bacteria.DNA.m[grepl("LV",bacteria.DNA.m$site2),]
bacteria.DNALV<-bacteria.DNALV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="mean"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,2,3)]
bacteria.DNA<-merge(bacteria.DNALV,bacteria.normalLV,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LV"
bacteria.DNA$pooling.effect<-round(bacteria.DNA$estimate)/round(bacteria.DNA$normalestimate)
bacteria.DNALV<-bacteria.DNA
##LW
bacteria.DNALW<-bacteria.DNA.m[grepl("LW",bacteria.DNA.m$site2),]
bacteria.DNALW<-bacteria.DNALW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="mean"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,2,3)]
bacteria.DNA<-merge(bacteria.DNALW,bacteria.normalLW,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LW"
bacteria.DNA$pooling.effect<-round((bacteria.DNA$estimate))/round(bacteria.DNA$normalestimate)
bacteria.DNALW<-bacteria.DNA
##LZ
bacteria.DNALZ<-bacteria.DNA.m[grepl("LZ",bacteria.DNA.m$site2),]
bacteria.DNALZ<-bacteria.DNALZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="mean"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,2,3)]
#
bacteria.DNA<-merge(bacteria.DNALZ,bacteria.normalLZ,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$pooling.effect<-(round(bacteria.DNA$estimate))/round(bacteria.DNA$normalestimate)
bacteria.DNALZ<-bacteria.DNA
bacteria.all2<-rbind(bacteria.DNALZ,bacteria.DNALV,bacteria.DNALW)
##separate the site
bacteria.all$pooling<-"soil"
bacteria.all2$pooling<-"DNA"
all<-rbind(bacteria.all,bacteria.all2)
all<-all[,c(1,4,5:7)]
all2<-aggregate(all[,4],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
bacteria.normal.ratio<-all2
bacteria.normal.ratio$type<-"low"
bacteria.normal.ratio$community<-"bacteria"
bacteria.pool.ratio$type<-"high"
##
bacteria.summary.m<-read.csv("proportion.bacteria.richness.summarydepthforpooled.csv",header = TRUE,row.names = 1,sep = ",")
infor<-bacteria.summary.m[,c(1,8,11)]
infor$mer<-paste0(infor$site,infor$method)
infor<-infor[!duplicated(infor),]
bacteria.normal.ratio$mer<-paste0(bacteria.normal.ratio$site2,bacteria.normal.ratio$method)
bacteria.normal.ratio<-merge(infor[,c(2,4)],bacteria.normal.ratio,by="mer")
bacteria.normal.ratio$pooled.seqs<-1067
bacteria.normal.ratio$proportion<-percent(bacteria.normal.ratio$pooled.seqs/bacteria.normal.ratio$unpooledn_seqs)
bacteria.pool.ratio$mer<-paste0(bacteria.pool.ratio$site2,bacteria.pool.ratio$method)
bacteria.pool.ratio<-merge(infor[,c(2,4)],bacteria.pool.ratio,by="mer")
bacteria.pool.ratio$pooled.seqs<-39290
bacteria.pool.ratio$proportion<-percent(bacteria.pool.ratio$pooled.seqs/bacteria.pool.ratio$unpooledn_seqs)
##
bacteria.summary.m<-bacteria.summary.m[c(1,4,6,8:11)]
bacteria.pool.ratio<-bacteria.pool.ratio[,c(2:6,9,10)]
bacteria.normal.ratio<-bacteria.normal.ratio[,c(2:6,9,10)]
names(bacteria.normal.ratio)[c(3,5)]<-c("site","ratio")
names(bacteria.pool.ratio)[c(3,5)]<-c("site","ratio")
names(bacteria.summary.m)[c(2,3)]<-c("pooled.seqs","pooling")
bacteria.normal.ratio$type<-"low"
bacteria.pool.ratio$type<-"high"
bacteria.summary.m$type<-"summary"
bacteria<-rbind(bacteria.normal.ratio,bacteria.pool.ratio,bacteria.summary.m)
bacteria$community<-"bacteria"

#animal
#divide each method
animal.normal.m <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
animal.normal.m$method<-gsub("40A","_40A",animal.normal.m$method)
animal.normal.m$method<-gsub("40B","_40B",animal.normal.m$method)
animal.normal.m$method[animal.normal.m$method=="GSMc"]<-"GSMc_62"
##animal soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.soil.m<-read.csv("animal.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
animal.soil.m<-animal.soil.m[grepl("soil",row.names(animal.soil.m)),]
row.names(animal.soil.m)<-gsub("soil","",row.names(animal.soil.m))
animal.soil.m$site<-row.names(animal.soil.m)

animal.soil.m<-merge(animal.soil.m,soil.metadata,by="site",all=T)
animal.soil.m$site2<-substr(animal.soil.m$site,1,2)
animal.soil.m$estimate<-round(animal.soil.m$estimate)
##visualize
animal.soilLV<-animal.soil.m[grepl("LV",animal.soil.m$site2),]
animal.soilLV<-animal.soilLV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="mean"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,2,3)]
animal.soil<-merge(animal.soilLV,animal.normalLV,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LV"
animal.soil$pooling.effect<-round(animal.soil$estimate)/round(animal.soil$normalestimate)
animal.soilLV<-animal.soil
##LW
animal.soilLW<-animal.soil.m[grepl("LW",animal.soil.m$site2),]
animal.soilLW<-animal.soilLW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="mean"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,2,3)]
animal.soil<-merge(animal.soilLW,animal.normalLW,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LW"
animal.soil$pooling.effect<-round((animal.soil$estimate))/round(animal.soil$normalestimate)
animal.soilLW<-animal.soil
##LZ
animal.soilLZ<-animal.soil.m[grepl("LZ",animal.soil.m$site2),]
animal.soilLZ<-animal.soilLZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="mean"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,2,3)]
#
animal.soil<-merge(animal.soilLZ,animal.normalLZ,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$pooling.effect<-(round(animal.soil$estimate))/round(animal.soil$normalestimate)
animal.soilLZ<-animal.soil
animal.all<-rbind(animal.soilLZ,animal.soilLV,animal.soilLW)
##animal DNA
animal.normal.m <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
animal.normal.m$method<-gsub("40A","_40A",animal.normal.m$method)
animal.normal.m$method<-gsub("40B","_40B",animal.normal.m$method)
animal.normal.m$method[animal.normal.m$method=="GSMc"]<-"GSMc_62"
##animal DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.DNA.m<-read.csv("animal.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
animal.DNA.m<-animal.DNA.m[grepl("DNA",row.names(animal.DNA.m)),]
row.names(animal.DNA.m)<-gsub("DNA","",row.names(animal.DNA.m))
animal.DNA.m$site<-row.names(animal.DNA.m)

animal.DNA.m<-merge(animal.DNA.m,DNA.metadata,by="site",all=T)
animal.DNA.m$site2<-substr(animal.DNA.m$site,1,2)
animal.DNA.m$estimate<-round(animal.DNA.m$estimate)
##visualize
animal.DNALV<-animal.DNA.m[grepl("LV",animal.DNA.m$site2),]
animal.DNALV<-animal.DNALV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="mean"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,2,3)]
animal.DNA<-merge(animal.DNALV,animal.normalLV,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LV"
animal.DNA$pooling.effect<-round(animal.DNA$estimate)/round(animal.DNA$normalestimate)
animal.DNALV<-animal.DNA
##LW
animal.DNALW<-animal.DNA.m[grepl("LW",animal.DNA.m$site2),]
animal.DNALW<-animal.DNALW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="mean"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,2,3)]
animal.DNA<-merge(animal.DNALW,animal.normalLW,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LW"
animal.DNA$pooling.effect<-round((animal.DNA$estimate))/round(animal.DNA$normalestimate)
animal.DNALW<-animal.DNA
##LZ
animal.DNALZ<-animal.DNA.m[grepl("LZ",animal.DNA.m$site2),]
animal.DNALZ<-animal.DNALZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="mean"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,2,3)]
#
animal.DNA<-merge(animal.DNALZ,animal.normalLZ,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$pooling.effect<-(round(animal.DNA$estimate))/round(animal.DNA$normalestimate)
animal.DNALZ<-animal.DNA
animal.all2<-rbind(animal.DNALZ,animal.DNALV,animal.DNALW)
##separate the site
animal.all$pooling<-"soil"
animal.all2$pooling<-"DNA"
all<-rbind(animal.all,animal.all2)
all<-all[,c(1,4,5:7)]
all2<-aggregate(all[,4],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
animal.pool.ratio<-all2
animal.pool.ratio$community<-"animal"
animal.pool.ratio$type<-"high"
##low sequencing
#divide each method
animal.normal.m <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
animal.normal.m$method<-gsub("40A","_40A",animal.normal.m$method)
animal.normal.m$method<-gsub("40B","_40B",animal.normal.m$method)
animal.normal.m$method[animal.normal.m$method=="GSMc"]<-"GSMc_62"
##animal soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.soil.m<-read.csv("animal.richness.csv",header = TRUE,row.names = 1,sep = ",")
animal.soil.m<-animal.soil.m[grepl("soil",row.names(animal.soil.m)),]
row.names(animal.soil.m)<-gsub("soil","",row.names(animal.soil.m))
animal.soil.m$site<-row.names(animal.soil.m)

animal.soil.m<-merge(animal.soil.m,soil.metadata,by="site",all=T)
animal.soil.m$site2<-substr(animal.soil.m$site,1,2)
animal.soil.m$estimate<-round(animal.soil.m$estimate)
##visualize
animal.soilLV<-animal.soil.m[grepl("LV",animal.soil.m$site2),]
animal.soilLV<-animal.soilLV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="mean"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,2,3)]
animal.soil<-merge(animal.soilLV,animal.normalLV,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LV"
animal.soil$pooling.effect<-round(animal.soil$estimate)/round(animal.soil$normalestimate)
animal.soilLV<-animal.soil
##LW
animal.soilLW<-animal.soil.m[grepl("LW",animal.soil.m$site2),]
animal.soilLW<-animal.soilLW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="mean"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,2,3)]
animal.soil<-merge(animal.soilLW,animal.normalLW,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LW"
animal.soil$pooling.effect<-round((animal.soil$estimate))/round(animal.soil$normalestimate)
animal.soilLW<-animal.soil
##LZ
animal.soilLZ<-animal.soil.m[grepl("LZ",animal.soil.m$site2),]
animal.soilLZ<-animal.soilLZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="mean"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,2,3)]
#
animal.soil<-merge(animal.soilLZ,animal.normalLZ,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$pooling.effect<-(round(animal.soil$estimate))/round(animal.soil$normalestimate)
animal.soilLZ<-animal.soil
animal.all<-rbind(animal.soilLZ,animal.soilLV,animal.soilLW)
##animal DNA
animal.normal.m <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
animal.normal.m$method<-gsub("40A","_40A",animal.normal.m$method)
animal.normal.m$method<-gsub("40B","_40B",animal.normal.m$method)
animal.normal.m$method[animal.normal.m$method=="GSMc"]<-"GSMc_62"
##animal DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.DNA.m<-read.csv("animal.richness.csv",header = TRUE,row.names = 1,sep = ",")
animal.DNA.m<-animal.DNA.m[grepl("DNA",row.names(animal.DNA.m)),]
row.names(animal.DNA.m)<-gsub("DNA","",row.names(animal.DNA.m))
animal.DNA.m$site<-row.names(animal.DNA.m)

animal.DNA.m<-merge(animal.DNA.m,DNA.metadata,by="site",all=T)
animal.DNA.m$site2<-substr(animal.DNA.m$site,1,2)
animal.DNA.m$estimate<-round(animal.DNA.m$estimate)
##visualize
animal.DNALV<-animal.DNA.m[grepl("LV",animal.DNA.m$site2),]
animal.DNALV<-animal.DNALV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="mean"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,2,3)]
animal.DNA<-merge(animal.DNALV,animal.normalLV,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LV"
animal.DNA$pooling.effect<-round(animal.DNA$estimate)/round(animal.DNA$normalestimate)
animal.DNALV<-animal.DNA
##LW
animal.DNALW<-animal.DNA.m[grepl("LW",animal.DNA.m$site2),]
animal.DNALW<-animal.DNALW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="mean"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,2,3)]
animal.DNA<-merge(animal.DNALW,animal.normalLW,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LW"
animal.DNA$pooling.effect<-round((animal.DNA$estimate))/round(animal.DNA$normalestimate)
animal.DNALW<-animal.DNA
##LZ
animal.DNALZ<-animal.DNA.m[grepl("LZ",animal.DNA.m$site2),]
animal.DNALZ<-animal.DNALZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="mean"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,2,3)]
#
animal.DNA<-merge(animal.DNALZ,animal.normalLZ,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$pooling.effect<-(round(animal.DNA$estimate))/round(animal.DNA$normalestimate)
animal.DNALZ<-animal.DNA
animal.all2<-rbind(animal.DNALZ,animal.DNALV,animal.DNALW)
##separate the site
animal.all$pooling<-"soil"
animal.all2$pooling<-"DNA"
all<-rbind(animal.all,animal.all2)
all<-all[,c(1,4,5:7)]
all2<-aggregate(all[,4],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
animal.normal.ratio<-all2
animal.normal.ratio$community<-"animal"
animal.normal.ratio$type<-"low"
animal.pool.ratio$type<-"high"
#summary
animal.summary.m<-read.csv("proportion.animal.richness.summarydepthforpooled.csv",header = TRUE,row.names = 1,sep = ",")
infor<-animal.summary.m[,c(1,8,11)]
infor$mer<-paste0(infor$site,infor$method)
infor<-infor[!duplicated(infor),]
animal.normal.ratio$mer<-paste0(animal.normal.ratio$site2,animal.normal.ratio$method)
animal.normal.ratio<-merge(infor[,c(2,4)],animal.normal.ratio,by="mer")
animal.normal.ratio$pooled.seqs<-880
animal.normal.ratio$proportion<-percent(animal.normal.ratio$pooled.seqs/animal.normal.ratio$unpooledn_seqs)
animal.pool.ratio$mer<-paste0(animal.pool.ratio$site2,animal.pool.ratio$method)
animal.pool.ratio<-merge(infor[,c(2,4)],animal.pool.ratio,by="mer")
animal.pool.ratio$pooled.seqs<-15962
animal.pool.ratio$proportion<-percent(animal.pool.ratio$pooled.seqs/animal.pool.ratio$unpooledn_seqs)
##
animal.summary.m<-animal.summary.m[c(1,4,6,8:11)]
animal.pool.ratio<-animal.pool.ratio[,c(2:6,9,10)]
animal.normal.ratio<-animal.normal.ratio[,c(2:6,9,10)]
names(animal.normal.ratio)[c(3,5)]<-c("site","ratio")
names(animal.pool.ratio)[c(3,5)]<-c("site","ratio")
names(animal.summary.m)[c(2,3)]<-c("pooled.seqs","pooling")
animal.normal.ratio$type<-"low"
animal.pool.ratio$type<-"high"
animal.summary.m$type<-"summary"
animal<-rbind(animal.normal.ratio,animal.pool.ratio,animal.summary.m)
animal$community<-"animal"

#fungi
#divide each method
fungi.normal.m <- read.csv("fungi.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal.m$method<-gsub("40A","_40A",fungi.normal.m$method)
fungi.normal.m$method<-gsub("40B","_40B",fungi.normal.m$method)
fungi.normal.m$method[fungi.normal.m$method=="GSMc"]<-"GSMc_62"
##fungi soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.soil.m<-read.csv("fungi.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
fungi.soil.m<-fungi.soil.m[grepl("soil",row.names(fungi.soil.m)),]
row.names(fungi.soil.m)<-gsub("soil","",row.names(fungi.soil.m))
fungi.soil.m$site<-row.names(fungi.soil.m)

fungi.soil.m<-merge(fungi.soil.m,soil.metadata,by="site",all=T)
fungi.soil.m$site2<-substr(fungi.soil.m$site,1,2)
fungi.soil.m$estimate<-round(fungi.soil.m$estimate)
##visualize
fungi.soilLV<-fungi.soil.m[grepl("LV",fungi.soil.m$site2),]
fungi.soilLV<-fungi.soilLV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="mean"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,2,3)]
fungi.soil<-merge(fungi.soilLV,fungi.normalLV,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LV"
fungi.soil$pooling.effect<-round(fungi.soil$estimate)/round(fungi.soil$normalestimate)
fungi.soilLV<-fungi.soil
##LW
fungi.soilLW<-fungi.soil.m[grepl("LW",fungi.soil.m$site2),]
fungi.soilLW<-fungi.soilLW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="mean"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,2,3)]
fungi.soil<-merge(fungi.soilLW,fungi.normalLW,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LW"
fungi.soil$pooling.effect<-round((fungi.soil$estimate))/round(fungi.soil$normalestimate)
fungi.soilLW<-fungi.soil
##LZ
fungi.soilLZ<-fungi.soil.m[grepl("LZ",fungi.soil.m$site2),]
fungi.soilLZ<-fungi.soilLZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="mean"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,2,3)]
#
fungi.soil<-merge(fungi.soilLZ,fungi.normalLZ,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$pooling.effect<-(round(fungi.soil$estimate))/round(fungi.soil$normalestimate)
fungi.soilLZ<-fungi.soil
fungi.all<-rbind(fungi.soilLZ,fungi.soilLV,fungi.soilLW)
##fungi DNA
fungi.normal.m <- read.csv("fungi.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal.m$method<-gsub("40A","_40A",fungi.normal.m$method)
fungi.normal.m$method<-gsub("40B","_40B",fungi.normal.m$method)
fungi.normal.m$method[fungi.normal.m$method=="GSMc"]<-"GSMc_62"
##fungi DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.DNA.m<-read.csv("fungi.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
fungi.DNA.m<-fungi.DNA.m[grepl("DNA",row.names(fungi.DNA.m)),]
row.names(fungi.DNA.m)<-gsub("DNA","",row.names(fungi.DNA.m))
fungi.DNA.m$site<-row.names(fungi.DNA.m)

fungi.DNA.m<-merge(fungi.DNA.m,DNA.metadata,by="site",all=T)
fungi.DNA.m$site2<-substr(fungi.DNA.m$site,1,2)
fungi.DNA.m$estimate<-round(fungi.DNA.m$estimate)
##visualize
fungi.DNALV<-fungi.DNA.m[grepl("LV",fungi.DNA.m$site2),]
fungi.DNALV<-fungi.DNALV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="mean"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,2,3)]
fungi.DNA<-merge(fungi.DNALV,fungi.normalLV,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LV"
fungi.DNA$pooling.effect<-round(fungi.DNA$estimate)/round(fungi.DNA$normalestimate)
fungi.DNALV<-fungi.DNA
##LW
fungi.DNALW<-fungi.DNA.m[grepl("LW",fungi.DNA.m$site2),]
fungi.DNALW<-fungi.DNALW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="mean"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,2,3)]
fungi.DNA<-merge(fungi.DNALW,fungi.normalLW,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LW"
fungi.DNA$pooling.effect<-round((fungi.DNA$estimate))/round(fungi.DNA$normalestimate)
fungi.DNALW<-fungi.DNA
##LZ
fungi.DNALZ<-fungi.DNA.m[grepl("LZ",fungi.DNA.m$site2),]
fungi.DNALZ<-fungi.DNALZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="mean"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,2,3)]
#
fungi.DNA<-merge(fungi.DNALZ,fungi.normalLZ,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$pooling.effect<-(round(fungi.DNA$estimate))/round(fungi.DNA$normalestimate)
fungi.DNALZ<-fungi.DNA
fungi.all2<-rbind(fungi.DNALZ,fungi.DNALV,fungi.DNALW)
##separate the site
fungi.all$pooling<-"soil"
fungi.all2$pooling<-"DNA"
all<-rbind(fungi.all,fungi.all2)
all<-all[,c(1,4,5:7)]
all2<-aggregate(all[,4],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
fungi.pool.ratio<-all2
fungi.pool.ratio$community<-"fungi"
fungi.pool.ratio$type<-"high"
##low.sequencing
#divide each method
fungi.normal.m <- read.csv("fungi.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal.m$method<-gsub("40A","_40A",fungi.normal.m$method)
fungi.normal.m$method<-gsub("40B","_40B",fungi.normal.m$method)
fungi.normal.m$method[fungi.normal.m$method=="GSMc"]<-"GSMc_62"
##fungi soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.soil.m<-read.csv("fungi.richness.csv",header = TRUE,row.names = 1,sep = ",")
fungi.soil.m<-fungi.soil.m[grepl("soil",row.names(fungi.soil.m)),]
row.names(fungi.soil.m)<-gsub("soil","",row.names(fungi.soil.m))
fungi.soil.m$site<-row.names(fungi.soil.m)

fungi.soil.m<-merge(fungi.soil.m,soil.metadata,by="site",all=T)
fungi.soil.m$site2<-substr(fungi.soil.m$site,1,2)
fungi.soil.m$estimate<-round(fungi.soil.m$estimate)
##visualize
fungi.soilLV<-fungi.soil.m[grepl("LV",fungi.soil.m$site2),]
fungi.soilLV<-fungi.soilLV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="mean"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,2,3)]
fungi.soil<-merge(fungi.soilLV,fungi.normalLV,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LV"
fungi.soil$pooling.effect<-round(fungi.soil$estimate)/round(fungi.soil$normalestimate)
fungi.soilLV<-fungi.soil
##LW
fungi.soilLW<-fungi.soil.m[grepl("LW",fungi.soil.m$site2),]
fungi.soilLW<-fungi.soilLW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="mean"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,2,3)]
fungi.soil<-merge(fungi.soilLW,fungi.normalLW,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LW"
fungi.soil$pooling.effect<-round((fungi.soil$estimate))/round(fungi.soil$normalestimate)
fungi.soilLW<-fungi.soil
##LZ
fungi.soilLZ<-fungi.soil.m[grepl("LZ",fungi.soil.m$site2),]
fungi.soilLZ<-fungi.soilLZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="mean"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,2,3)]
#
fungi.soil<-merge(fungi.soilLZ,fungi.normalLZ,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$pooling.effect<-(round(fungi.soil$estimate))/round(fungi.soil$normalestimate)
fungi.soilLZ<-fungi.soil
fungi.all<-rbind(fungi.soilLZ,fungi.soilLV,fungi.soilLW)
##fungi DNA
fungi.normal.m <- read.csv("fungi.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal.m$method<-gsub("40A","_40A",fungi.normal.m$method)
fungi.normal.m$method<-gsub("40B","_40B",fungi.normal.m$method)
fungi.normal.m$method[fungi.normal.m$method=="GSMc"]<-"GSMc_62"
##fungi DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.DNA.m<-read.csv("fungi.richness.csv",header = TRUE,row.names = 1,sep = ",")
fungi.DNA.m<-fungi.DNA.m[grepl("DNA",row.names(fungi.DNA.m)),]
row.names(fungi.DNA.m)<-gsub("DNA","",row.names(fungi.DNA.m))
fungi.DNA.m$site<-row.names(fungi.DNA.m)

fungi.DNA.m<-merge(fungi.DNA.m,DNA.metadata,by="site",all=T)
fungi.DNA.m$site2<-substr(fungi.DNA.m$site,1,2)
fungi.DNA.m$estimate<-round(fungi.DNA.m$estimate)
##visualize
fungi.DNALV<-fungi.DNA.m[grepl("LV",fungi.DNA.m$site2),]
fungi.DNALV<-fungi.DNALV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="mean"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,2,3)]
fungi.DNA<-merge(fungi.DNALV,fungi.normalLV,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LV"
fungi.DNA$pooling.effect<-round(fungi.DNA$estimate)/round(fungi.DNA$normalestimate)
fungi.DNALV<-fungi.DNA
##LW
fungi.DNALW<-fungi.DNA.m[grepl("LW",fungi.DNA.m$site2),]
fungi.DNALW<-fungi.DNALW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="mean"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,2,3)]
fungi.DNA<-merge(fungi.DNALW,fungi.normalLW,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LW"
fungi.DNA$pooling.effect<-round((fungi.DNA$estimate))/round(fungi.DNA$normalestimate)
fungi.DNALW<-fungi.DNA
##LZ
fungi.DNALZ<-fungi.DNA.m[grepl("LZ",fungi.DNA.m$site2),]
fungi.DNALZ<-fungi.DNALZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="mean"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,2,3)]
#
fungi.DNA<-merge(fungi.DNALZ,fungi.normalLZ,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$pooling.effect<-(round(fungi.DNA$estimate))/round(fungi.DNA$normalestimate)
fungi.DNALZ<-fungi.DNA
fungi.all2<-rbind(fungi.DNALZ,fungi.DNALV,fungi.DNALW)
##separate the site
fungi.all$pooling<-"soil"
fungi.all2$pooling<-"DNA"
all<-rbind(fungi.all,fungi.all2)
all<-all[,c(1,4,5:7)]
all2<-aggregate(all[,4],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
fungi.normal.ratio<-all2
fungi.normal.ratio$community<-"fungi"
fungi.normal.ratio$type<-"low"
fungi.pool.ratio$type<-"high"
#summary
fungi.summary.m<-read.csv("proportion.fungi.richness.summarydepthforpooled.csv",header = TRUE,row.names = 1,sep = ",")
infor<-fungi.summary.m[,c(1,8,11)]
infor$mer<-paste0(infor$site,infor$method)
infor<-infor[!duplicated(infor),]
fungi.normal.ratio$mer<-paste0(fungi.normal.ratio$site2,fungi.normal.ratio$method)
fungi.normal.ratio<-merge(infor[,c(2,4)],fungi.normal.ratio,by="mer")
fungi.normal.ratio$pooled.seqs<-814
fungi.normal.ratio$proportion<-percent(fungi.normal.ratio$pooled.seqs/fungi.normal.ratio$unpooledn_seqs)
fungi.pool.ratio$mer<-paste0(fungi.pool.ratio$site2,fungi.pool.ratio$method)
fungi.pool.ratio<-merge(infor[,c(2,4)],fungi.pool.ratio,by="mer")
fungi.pool.ratio$pooled.seqs<-14409
fungi.pool.ratio$proportion<-percent(fungi.pool.ratio$pooled.seqs/fungi.pool.ratio$unpooledn_seqs)
##
fungi.summary.m<-fungi.summary.m[c(1,4,6,8:11)]
fungi.pool.ratio<-fungi.pool.ratio[,c(2:6,9,10)]
fungi.normal.ratio<-fungi.normal.ratio[,c(2:6,9,10)]
names(fungi.normal.ratio)[c(3,5)]<-c("site","ratio")
names(fungi.pool.ratio)[c(3,5)]<-c("site","ratio")
names(fungi.summary.m)[c(2,3)]<-c("pooled.seqs","pooling")
fungi.normal.ratio$type<-"low"
fungi.pool.ratio$type<-"high"
fungi.summary.m$type<-"summary"
fungi<-rbind(fungi.normal.ratio,fungi.pool.ratio,fungi.summary.m)
fungi$community<-"fungi"
##
all<-rbind(animal,bacteria,fungi)
write.csv(all,"proportion.all.csv")
##animal
animal$number<-as.numeric(sub("%", "", animal$proportion)) / 100
animal$try<-NA
animal[animal$number<0.25,]$try<-"<25%"
animal[animal$number>=0.25 & animal$number< 0.5,]$try<-"=25%~<50%"
animal[animal$number>=0.50&animal$number<0.75,]$try<-"=50%~<75%"
animal[animal$number==0.75,]$try<-"75%"
animal[animal$number==1,]$try<-"100%"
animal[animal$number>1,]$try<-">100%"
soil<-animal[grepl("soil",animal$pooling),]
soil<-soil[!is.na(soil$ratio),]
soil <- soil %>%
  filter(!(method == "deep"))
soil <- soil %>%
  filter(!(method == "deep_SUCC"))
mod01<-lm(ratio~ try+method+ site,data = soil)
#broom::tidy(mod01)
options(scipen = 999)
result.soil<-parameters::model_parameters(mod01)
result.soil$pooling<-'soil'
result.soil$organism<-'animal'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
letter$Group[3]<-"=25%~<50%"
letter$Group[4]<-"=50%~<75%"
letter$Group[5]<-">100%"
difference <- soil%>%
  group_by(try) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="try")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "try")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "<25%","=25%~<50%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
#DNA
DNA<-animal[grepl("DNA",animal$pooling),]
DNA<-DNA[!is.na(DNA$ratio),]
DNA <- DNA %>%
  filter(!(method == "deep"))
DNA <- DNA %>%
  filter(!(method == "deep_SUCC"))
mod01<-lm(ratio~ try+method+ site,data = DNA)
performance::performance(mod01)
options(scipen = 999)
result.DNA<-parameters::model_parameters(mod01)
result.DNA$pooling<-'DNA'
result.DNA$organism<-'animal'
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
letter$Group[3]<-"=25%~<50%"
letter$Group[4]<-"=50%~<75%"
letter$Group[5]<-">100%"
difference <- DNA%>%
  group_by(try) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="try")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "try")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "<25%","=25%~<50%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##compare the pooling type
pooling<-rbind(soil,DNA)
mod01<-lmer(ratio~ pooling+method+(1|try)+ (1|site),data = pooling)
performance::performance(mod01)
options(scipen = 999)
result.pooling<-parameters::model_parameters(mod01)
result.pooling$organism<-'animal'
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
difference <- DNA%>%
  group_by(pooling) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="pooling")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "pooling")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##fungi
fungi$proportion[fungi$proportion=="1 770.1%"]<-c("1770.1%","1770.1%")
fungi$number<-as.numeric(sub("%", "", fungi$proportion)) / 100
fungi$try<-NA
fungi[fungi$number<0.25,]$try<-"<25%"
fungi[fungi$number>=0.25 & fungi$number< 0.5,]$try<-"=25%~<50%"
fungi[fungi$number>=0.50&fungi$number<0.75,]$try<-"=50%~<75%"
fungi[fungi$number==0.75,]$try<-"75%"
fungi[fungi$number==1,]$try<-"100%"
fungi[fungi$number>1,]$try<-">100%"
fungi <- fungi %>%
  filter(!(site == "LV" & method == "LUCAS"))
fungi <- fungi %>%
  filter(!(site == "LV" & method == "deep"))
fungi <- fungi %>%
  filter(!(site == "LV" & method == "deep_SUCC"))
soil<-fungi[grepl("soil",fungi$pooling),]
soil<-soil[!is.na(soil$ratio),]
soil<-soil[!is.na(soil$try),]
mod01<-lm(ratio~ try+method+site,data = soil)
performance::performance(mod01)
options(scipen = 999)
fungi.result.soil<-parameters::model_parameters(mod01)
fungi.result.soil$pooling<-'soil'
fungi.result.soil$organism<-'fungi'
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
letter$Group[3]<-"=25%~<50%"
letter$Group[4]<-"=50%~<75%"
letter$Group[5]<-">100%"
difference <- soil%>%
  group_by(try) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="try")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "try")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "<25%","=25%~<50%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  scale_y_continuous(breaks = c(0.5,1,1.5,2))+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
#DNA
DNA<-fungi[grepl("DNA",fungi$pooling),]
DNA<-DNA[!is.na(DNA$ratio),]
DNA<-DNA[!is.na(DNA$try),]
mod01<-lm(ratio~ try+method+site,data = DNA)
performance::performance(mod01)
options(scipen = 999)
fungi.result.DNA<-parameters::model_parameters(mod01)
fungi.result.DNA$pooling<-'DNA'
fungi.result.DNA$organism<-'fungi'
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
letter$Group[3]<-"=25%~<50%"
letter$Group[4]<-"=50%~<75%"
letter$Group[5]<-">100%"
difference <- DNA%>%
  group_by(try) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="try")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "try")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "<25%","=25%~<50%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  scale_y_continuous(breaks = c(0.5,1,1.5,2))+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##compare the pooling type
pooling<-rbind(soil,DNA)
mod01<-lmer(ratio~ pooling+method+(1|try)+ (1|site),data = pooling)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
performance::performance(mod01)
options(scipen = 999)
fungi.result.pooling<-parameters::model_parameters(mod01)
fungi.result.pooling$organism<-'fungi'
eta_squared(mod01)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
difference <- DNA%>%
  group_by(pooling) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="pooling")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "pooling")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##bacteria
bacteria$number<-as.numeric(sub("%", "", bacteria$proportion)) / 100
bacteria$try<-NA
bacteria[bacteria$number<0.25,]$try<-"<25%"
bacteria[bacteria$number>=0.25 & bacteria$number< 0.5,]$try<-"=25%~<50%"
bacteria[bacteria$number>=0.50&bacteria$number<0.75,]$try<-"=50%~<75%"
bacteria[bacteria$number==0.75,]$try<-"75%"
bacteria[bacteria$number<0.95&bacteria$number>0.75,]$try<-"<75%~95%"
bacteria[bacteria$number==1,]$try<-"100%"
bacteria[bacteria$number>1,]$try<-">100%"
soil<-bacteria[grepl("soil",bacteria$pooling),]
soil<-soil[!is.na(soil$ratio),]
soil<-soil[!is.na(soil$try),]
mod01<-lm(ratio~ try+method+site,data = soil)
performance::performance(mod01)
options(scipen = 999)
bacteria.result.soil<-parameters::model_parameters(mod01)
bacteria.result.soil$pooling<-'soil'
bacteria.result.soil$organism<-'bacteria'
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
letter$Group[4]<-"=25%~<50%"
letter$Group[5]<-"=50%~<75%"
letter$Group[6]<-">100%"
difference <- soil%>%
  group_by(try) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="try")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "try")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "<25%","=25%~<50%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5))+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
#DNA
DNA<-bacteria[grepl("DNA",bacteria$pooling),]
DNA<-DNA[!is.na(DNA$ratio),]
mod01<-lm(ratio~ try+method+ site,data = DNA)
performance::performance(mod01)
options(scipen = 999)
bacteria.result.DNA<-parameters::model_parameters(mod01)
bacteria.result.DNA$pooling<-'DNA'
bacteria.result.DNA$organism<-'bacteria'
DNA<-DNA[!is.na(DNA$try),]
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
letter$Group[4]<-"=25%~<50%"
letter$Group[5]<-"=50%~<75%"
letter$Group[6]<-">100%"
difference <- DNA%>%
  group_by(try) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="try")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "try")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "<25%","=25%~<50%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  scale_y_continuous(breaks = c(0.5,1,1.5,2,2.5))+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##compare the pooling type
pooling<-rbind(soil,DNA)
mod01<-lmer(ratio~ pooling+method+(1|try)+ (1|site),data = pooling)
performance::performance(mod01)
options(scipen = 999)
bacteria.result.pooling<-parameters::model_parameters(mod01)
bacteria.result.pooling$organism<-'bacteria'
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
difference <- DNA%>%
  group_by(pooling) %>%
  summarise(max=mean(ratio))
y.site<-merge(letter,difference,by.x="Group",by.y="pooling")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)
plot_predictions(mod01,condition = "pooling")+
  labs(x="proportion",y="Ratio")+
  theme_light() +
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
result<-rbind(bacteria.result.DNA,bacteria.result.soil,
      fungi.result.DNA,fungi.result.soil,
      result.DNA,result.soil)
write.csv(result,"richness.parameters.csv")
result<-rbind(bacteria.result.pooling,
              fungi.result.pooling,
              result.pooling)
write.csv(result,"pooling.compare.parameters.csv")

