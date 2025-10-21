##animal
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
animal.normal.m<-read.csv("unpooled.animal.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m<-animal.normal.m[,c(1,6,7)]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lm(log(shannon)~  method + site,data = animal.normal.m)
options(scipen = 999)
animal.result<-parameters::model_parameters(mod01)
animal.result$organism<-'animal'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- animal.normal.m%>%
  group_by(method) %>%
  summarise(max=mean(log(shannon)))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks = c(1,2.5,5,5.5))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##fungi
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
fungi.normal.m<-read.csv("unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m<-fungi.normal.m[,c(1,6,7)]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lm(log(shannon)~  method + site,data = fungi.normal.m)
options(scipen = 999)
fungi.result<-parameters::model_parameters(mod01)
fungi.result$organism<-'fungi'
performance::performance(mod01)

b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- fungi.normal.m%>%
  group_by(method) %>%
  summarise(max=mean(log(shannon)))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks = c(1,2.5,5,5.5,6))+
  geom_hline(yintercept=1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##bacteria
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
bacteria.normal.m<-read.csv("unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m<-bacteria.normal.m[,c(1,6,7)]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lm(log(shannon)~  method + site,data = bacteria.normal.m)
bacteria.result<-parameters::model_parameters(mod01)
bacteria.result$organism<-'bacteria'
performance::performance(mod01)

b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- bacteria.normal.m%>%
  group_by(method) %>%
  summarise(max=mean(log(shannon)))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks = c(1,3,4,5,5.5,6))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

result<-rbind(bacteria.result,
              fungi.result,
              animal.result)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
write.csv(result,"a.unpooled.totalshannon.parameters.csv")

