##bacteria
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
bacteria.normal.m<-read.csv("bacteria.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m$site<-row.names(bacteria.normal.m)
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata2<-metadata
metadata$site<-paste0(metadata$site,"DNA")
metadata2$site<-paste0(metadata2$site,"soil")
metadata<-rbind(metadata,metadata2)
bacteria.normal.m<-merge(bacteria.normal.m,metadata,by="site",all=T)
bacteria.normal.m <- bacteria.normal.m[!is.na(bacteria.normal.m$estimate),]
bacteria.normal.m <- bacteria.normal.m[!is.na(bacteria.normal.m$method),]
bacteria.normal.m$site2<-substr(bacteria.normal.m$site,1,2)

bacteria.normal.m<-bacteria.normal.m[,c(1,2,8,9)]
###DNA
DNA<-bacteria.normal.m[grepl("DNA",bacteria.normal.m$site),]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
bacteria.result.DNA<-parameters::model_parameters(mod01)
bacteria.result.DNA$organism<-'bacteria'
bacteria.result.DNA$pooling<-'DNA'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- DNA%>%
  group_by(method.x) %>%
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###soil
soil<-bacteria.normal.m[grepl("soil",bacteria.normal.m$site),]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
bacteria.result.soil<-parameters::model_parameters(mod01)
bacteria.result.soil$organism<-'bacteria'
bacteria.result.soil$pooling<-'soil'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- soil%>%
  group_by(method.x) %>%
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
###compare soil and DNA
DNA$type<-"DNA"
soil$type<-"soil"
bacteria2<-rbind(DNA,soil)
mod01<-lmer(estimate.x~  type+(1|method.x) + (1|site2),data = bacteria2)
bacteria.result.DNA<-parameters::model_parameters(mod01)
bacteria.result.DNA$organism<-'bacteria'
bacteria.result.DNA$pooling<-'DNA'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(type = "pairwise")) 
plot_predictions(mod01,condition = "type", transform = exp )+
  labs(x="method",y="shannon")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

##fungi #use the mean value
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
fungi.normal.m<-read.csv("fungi.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m$site<-row.names(fungi.normal.m)
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata2<-metadata
metadata$site<-paste0(metadata$site,"DNA")
metadata2$site<-paste0(metadata2$site,"soil")
metadata<-rbind(metadata,metadata2)
fungi.normal.m<-merge(fungi.normal.m,metadata,by="site",all=T)
fungi.normal.m <- fungi.normal.m[!is.na(fungi.normal.m$estimate),]
fungi.normal.m <- fungi.normal.m[!is.na(fungi.normal.m$method),]
fungi.normal.m$site2<-substr(fungi.normal.m$site,1,2)

fungi.normal.m<-fungi.normal.m[,c(1,2,8,9)]
###DNA
DNA<-fungi.normal.m[grepl("DNA",fungi.normal.m$site),]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x +site2,data = DNA)
fungi.result.DNA<-parameters::model_parameters(mod01)
fungi.result.DNA$organism<-'fungi'
fungi.result.DNA$pooling<-'DNA'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- DNA%>%
  group_by(method.x) %>%
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(limits=c(3,6.2))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###soil
soil<-fungi.normal.m[grepl("soil",fungi.normal.m$site),]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
fungi.result.soil<-parameters::model_parameters(mod01)
fungi.result.soil$organism<-'fungi'
fungi.result.soil$pooling<-'soil'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- soil%>%
  group_by(method.x) %>%
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###compare soil and DNA
DNA$type<-"DNA"
soil$type<-"soil"
fungi2<-rbind(DNA,soil)
mod01<-lmer(log(estimate.x)~  type+(1|method.x) + (1|site2),data = fungi2)

b<-avg_comparisons(mod01, variables = list(type = "pairwise")) 
plot_predictions(mod01,condition = "type", transform = exp )+
  labs(x="method",y="shannon")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
##animal #use the mean value
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
animal.normal.m<-read.csv("animal.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m$site<-row.names(animal.normal.m)
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata2<-metadata
metadata$site<-paste0(metadata$site,"DNA")
metadata2$site<-paste0(metadata2$site,"soil")
metadata<-rbind(metadata,metadata2)
animal.normal.m<-merge(animal.normal.m,metadata,by="site",all=T)
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$estimate),]
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$method),]
animal.normal.m$site2<-substr(animal.normal.m$site,1,2)

animal.normal.m<-animal.normal.m[,c(1,2,8,9)]
###DNA
DNA<-animal.normal.m[grepl("DNA",animal.normal.m$site),]
DNA<-DNA[DNA$method!="deep",]
#DNA<-DNA[DNA$method!="LUCAS",]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
options(scipen = 999)
animal.result.DNA<-parameters::model_parameters(mod01)
animal.result.DNA$organism<-'animal'
animal.result.DNA$pooling<-'DNA'
performance::performance(mod01)
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- DNA%>%
  group_by(method.x) %>%
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC")))+ 
  scale_y_continuous(limits = c(0,6))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###soil
soil<-animal.normal.m[grepl("soil",animal.normal.m$site),]
soil<-soil[soil$method!="deep",]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
animal.result.soil<-parameters::model_parameters(mod01)
animal.result.soil$organism<-'animal'
animal.result.soil$pooling<-'soil'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- soil%>%
  group_by(method.x) %>%
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_y_continuous(limits = c(0,6))+
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###compare soil and DNA
DNA$type<-"DNA"
soil$type<-"soil"
animal2<-rbind(DNA,soil)
mod01<-lmer(log(estimate.x)~  type+(1|method.x) + (1|site2),data = animal2)

b<-avg_comparisons(mod01, variables = list(type = "pairwise")) 
plot_predictions(mod01,condition = "type", transform = exp )+
  labs(x="method",y="shannon")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
result<-rbind(bacteria.result.DNA,bacteria.result.soil,
              fungi.result.DNA,fungi.result.soil,
              animal.result.DNA,animal.result.soil)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
write.csv(result,"a.pooled.shannon.parameters.csv")