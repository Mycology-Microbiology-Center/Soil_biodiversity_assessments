##animal #use the mean value
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
animal.normal.m<-read.csv("animal.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m$site<-row.names(animal.normal.m)
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/normal_metadata.csv",header = TRUE,sep = ",")
animal.normal.m<-merge(animal.normal.m,metadata,by="site",all=T)
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$estimate),]
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$method),]
animal.normal.m$site2<-substr(animal.normal.m$site,1,2)

animal.normal.m<-animal.normal.m[,c(1,2,5,10)]
animal.normal.m<-animal.normal.m[!grepl("B",animal.normal.m$site),]
###GSMc40
select<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/1.accumulative/sheet2_for_pooled.csv")
##separate the site
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-animal.normal.m[grepl(GSMc40A,animal.normal.m$site),]
A40$site<-paste0(A40$site,"40A")
A40$method<-paste0(A40$method,"40A")

B40<-animal.normal.m[grepl(GSMc40B,animal.normal.m$site),]
B40$site<-paste0(B40$site,"40B")
B40$method<-paste0(B40$method,"40B")
animal.normal.m<-rbind(A40,B40,animal.normal.m)
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(animal.normal.m)[names(animal.normal.m)=="method"]<-"method.x"
names(animal.normal.m)[names(animal.normal.m)=="estimate"]<-"estimate.x"
animal.normal.m<-animal.normal.m[!grepl("deep_SUCC",animal.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = animal.normal.m)
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
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "Zobel","GSMc","GSMc40A","GSMc40B","DarkDiv","SUCC","LUCAS","Deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##bacteria #use the mean value
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
bacteria.normal.m<-read.csv("bacteria.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m$site<-row.names(bacteria.normal.m)
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/normal_metadata.csv",header = TRUE,sep = ",")
bacteria.normal.m<-merge(bacteria.normal.m,metadata,by="site",all=T)
bacteria.normal.m <- bacteria.normal.m[!is.na(bacteria.normal.m$estimate),]
bacteria.normal.m <- bacteria.normal.m[!is.na(bacteria.normal.m$method),]
bacteria.normal.m$site2<-substr(bacteria.normal.m$site,1,2)

bacteria.normal.m<-bacteria.normal.m[,c(1,2,5,10)]
bacteria.normal.m<-bacteria.normal.m[!grepl("B",bacteria.normal.m$site),]
###GSMc40
##separate the site
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-bacteria.normal.m[grepl(GSMc40A,bacteria.normal.m$site),]
A40$site<-paste0(A40$site,"40A")
A40$method<-paste0(A40$method,"40A")

B40<-bacteria.normal.m[grepl(GSMc40B,bacteria.normal.m$site),]
B40$site<-paste0(B40$site,"40B")
B40$method<-paste0(B40$method,"40B")
bacteria.normal.m<-rbind(A40,B40,bacteria.normal.m)
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(bacteria.normal.m)[names(bacteria.normal.m)=="method"]<-"method.x"
names(bacteria.normal.m)[names(bacteria.normal.m)=="estimate"]<-"estimate.x"
bacteria.normal.m<-bacteria.normal.m[!grepl("deep_SUCC",bacteria.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = bacteria.normal.m)
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
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_y_continuous(breaks=c(4.1,4.5,4.8))+ 
  scale_x_discrete(limits=as.character(c( "Zobel","GSMc","GSMc40A","GSMc40B","DarkDiv","SUCC","LUCAS","Deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##fungi #use the mean value
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
fungi.normal.m<-read.csv("fungi.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m$site<-row.names(fungi.normal.m)
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/normal_metadata.csv",header = TRUE,sep = ",")
fungi.normal.m<-merge(fungi.normal.m,metadata,by="site",all=T)
fungi.normal.m <- fungi.normal.m[!is.na(fungi.normal.m$estimate),]
fungi.normal.m <- fungi.normal.m[!is.na(fungi.normal.m$method),]
fungi.normal.m$site2<-substr(fungi.normal.m$site,1,2)

fungi.normal.m<-fungi.normal.m[,c(1,2,5,10)]
fungi.normal.m<-fungi.normal.m[!grepl("B",fungi.normal.m$site),]
###GSMc40
##separate the site
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-fungi.normal.m[grepl(GSMc40A,fungi.normal.m$site),]
A40$site<-paste0(A40$site,"40A")
A40$method<-paste0(A40$method,"40A")

B40<-fungi.normal.m[grepl(GSMc40B,fungi.normal.m$site),]
B40$site<-paste0(B40$site,"40B")
B40$method<-paste0(B40$method,"40B")
fungi.normal.m<-rbind(A40,B40,fungi.normal.m)
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(fungi.normal.m)[names(fungi.normal.m)=="method"]<-"method.x"
names(fungi.normal.m)[names(fungi.normal.m)=="estimate"]<-"estimate.x"
fungi.normal.m<-fungi.normal.m[!grepl("deep_SUCC",fungi.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = fungi.normal.m)
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
  summarise(max=mean(estimate.x))

y.site<-merge(letter,difference,by.x="Group",by.y="method.x")
y.site<-data.frame(area=factor(y.site$Group),ymax=y.site$max,letter=y.site$Letter)

plot_predictions(mod01,condition = "method.x" )+
  labs(x="method",y="shannon")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "Zobel","GSMc","GSMc40A","GSMc40B","DarkDiv","SUCC","LUCAS","Deep")))+ 
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
write.csv(result,"a.unpooled.shannon.parameters.csv")
