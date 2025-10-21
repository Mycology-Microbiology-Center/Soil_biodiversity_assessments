##animal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
#DNA
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/14.sequencing depth/correct1")
all<-read.csv("proportion.shannon.all.csv",header = TRUE,row.names = 1,sep = ",")
all<-all[grepl("100%",all$proportion),]
all<-all[,c(2:5,9)]
bacteria<-all[grepl("bacteria",all$community),]
fungi<-all[grepl("fungi",all$community),]
animal<-all[grepl("animal",all$community),]
DNA<-animal[animal$pooling=="DNA",]
DNA<-DNA[DNA$method!="deep",]
DNA<-DNA[DNA$method!="deep_SUCC",]
mod01<-lm(ratio ~ method + site,data = DNA)
options(scipen = 999)
animal.result.DNA<-parameters::model_parameters(mod01)
animal.result.DNA$organism<-'animal'
animal.result.DNA$pooling<-'DNA'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- DNA%>%
  group_by(method) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="Richness ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
#soil
soil<-animal[animal$pooling=="soil",]
soil<-soil[soil$method!="deep",]
soil<-soil[soil$method!="deep_SUCC",]
mod01<-lm(ratio ~ method + site,data = soil)
animal.result.soil<-parameters::model_parameters(mod01)
animal.result.soil$organism<-'animal'
animal.result.soil$pooling<-'soil'
performance::performance(mod01)

b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- soil%>%
  group_by(method) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="Richness ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3))+ 
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##fungi
##pooling effect
#DNA
DNA<-fungi[fungi$pooling=="DNA",]
DNA <- DNA %>%
  filter(!(site == "LV" & method == "LUCAS")) ###because of the unbalance sample size (eg., LUCAS: 1 vs.5)
DNA <- DNA %>%
  filter(!(site == "LV" & method == "deep"))
DNA <- DNA %>%
  filter(!(site == "LV" & method == "deep_SUCC"))
mod01<-lm(ratio ~ method + site,data = DNA)
fungi.result.DNA<-parameters::model_parameters(mod01)
fungi.result.DNA$organism<-'fungi'
fungi.result.DNA$pooling<-'DNA'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- DNA%>%
  group_by(method) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="Richness ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks=c(-1,0,0.7,1,1.5,2,3))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
#soil
soil<-fungi[fungi$pooling=="soil",]
soil <- soil %>%
  filter(!(site == "LV" & method == "LUCAS"))
soil <- soil %>%
  filter(!(site == "LV" & method == "deep"))
soil <- soil %>%
  filter(!(site == "LV" & method == "deep_SUCC"))
mod01<-lm(ratio ~ method + site,data = soil)
fungi.result.soil<-parameters::model_parameters(mod01)
fungi.result.soil$organism<-'fungi'
fungi.result.soil$pooling<-'soil'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- soil%>%
  group_by(method) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="Richness ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  scale_y_continuous(breaks=c(-1,0,0.5,1,1.5,2,3))+ 
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##bacteria
DNA<-bacteria[bacteria$pooling=="DNA",]
mod01<-lm(ratio ~ method + site,data = DNA)
bacteria.result.DNA<-parameters::model_parameters(mod01)
bacteria.result.DNA$organism<-'bacteria'
bacteria.result.DNA$pooling<-'DNA'
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- DNA%>%
  group_by(method) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="Richness ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks=c(0,1,2,3))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
#soil
soil<-bacteria[bacteria$pooling=="soil",]
mod01<-lm(ratio ~ method + site,data = soil)
bacteria.result.soil<-parameters::model_parameters(mod01)
bacteria.result.soil$organism<-'bacteria'
bacteria.result.soil$pooling<-'soil'
performance::performance(mod01)

b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group<-gsub("4A","40A",letter$Group)
letter$Group<-gsub("4B","40B",letter$Group)

difference <- soil%>%
  group_by(method) %>%
  summarise(max=mean(ratio))

y.site<-merge(letter,difference,by.x="Group",by.y="method")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method" )+
  labs(x="method",y="Richness ratio")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  scale_y_continuous(breaks=c(0,1,2))+ 
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
result<-rbind(bacteria.result.DNA,bacteria.result.soil,
              fungi.result.DNA,fungi.result.soil,
              animal.result.DNA,animal.result.soil)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/14.sequencing depth/correct1/revise")
write.csv(result,"PE.methods.shannon.csv")
