##animal
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
animal.normal.m<-read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
library("parameters")
names(animal.normal.m)[names(animal.normal.m)=="method"]<-"method.x"
mod01<-lm(mean~  method.x + site,data = animal.normal.m)
options(scipen = 999)
animal.result<-parameters::model_parameters(mod01)
animal.result$organism<-'animal'
performance::performance(mod01)

b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- animal.normal.m%>%
  group_by(method.x) %>%
  summarise(max=mean(mean))
y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x")+
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc","GSMc40A","GSMc40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##fungi
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
fungi.normal.m<-read.csv("fungi.totalrichness.csv",header = TRUE,row.names = 1,sep = ",")
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(fungi.normal.m)[names(fungi.normal.m)=="method"]<-"method.x"
mod01<-lm(mean~  method.x + site,data = fungi.normal.m)
options(scipen = 999)
fungi.result<-parameters::model_parameters(mod01)
fungi.result$organism<-'fungi'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- fungi.normal.m%>%
  group_by(method.x) %>%
  summarise(max=mean(mean))
y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x")+
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc","GSMc40A","GSMc40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##bacteria
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
bacteria.normal.m<-read.csv("bacteria.totalrichness.csv",header = TRUE,row.names = 1,sep = ",")
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(bacteria.normal.m)[names(bacteria.normal.m)=="method"]<-"method.x"
mod01<-lm(mean~  method.x + site,data = bacteria.normal.m)
options(scipen = 999)
bacteria.result<-parameters::model_parameters(mod01)
bacteria.result$organism<-'bacteria'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method.x = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- bacteria.normal.m%>%
  group_by(method.x) %>%
  summarise(max=mean(mean))
y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x")+
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc","GSMc40A","GSMc40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
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
write.csv(result,"a.unpooled.totalrichness.parameters.csv")

