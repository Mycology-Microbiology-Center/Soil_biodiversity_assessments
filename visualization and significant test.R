library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
library(vegan)
library(scales)
library(stringr)
library(metagMisc)
###accumulation curve visualization
##read the table
fungi<- read.table("fungi.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
bacteria<- read.table("bacteria.rarefy.table.csv",header = T,row.names = 1,sep = ",")
animal<- read.table("animal.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
##animal accumulation curve
animal2<-as.data.frame(t(animal))
animal2<-animal2[!grepl("DNA|soil",row.names(animal2)),]
animal.LZ<-animal2[grepl("LZ",row.names(animal2)),]
animal.LZ<-animal.LZ[!grepl("LZ7",row.names(animal.LZ)),]
animal.LZ<-animal.LZ[!grepl("B",row.names(animal.LZ)),]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
rownames(metadata) <- metadata$site
metadata<-metadata[,1:2]
metadata.LZ<-metadata[match(row.names(animal.LZ),metadata$site),]
select<-read.csv("sheet2_for_pooled.csv")
##separate the site
animal.LZ<-animal.LZ[,colSums(animal.LZ)!=0]
animal.LZ<-as.data.frame(t(animal.LZ))
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
A40<-animal.LZ[,grepl(GSMc40A,names(animal.LZ))]
names(A40)<-paste0(names(A40),"40A")

B40<-animal.LZ[,grepl(GSMc40B,names(animal.LZ))]
names(B40)<-paste0(names(B40),"40B")

metadata.LZ<-rbind(
data.frame(method="GSMc_40A",site=names(A40)),
data.frame(method="GSMc_40B",site=names(B40)),
metadata.LZ)

try<-function(x,y){
  merge(x,y,by="OTU",all=T)
}
A40$OTU<-row.names(A40)
B40$OTU<-row.names(B40)

animal.LZ$OTU<-row.names(animal.LZ)
a<-list(A40,B40,animal.LZ)
animal.LZ<-Reduce(try,a)
row.names(animal.LZ)<-animal.LZ$OTU
animal.LZ<-animal.LZ[,-1]
animal.LZ<-as.data.frame(t(animal.LZ))
methods<-unique(metadata.LZ$method)
accu.LZ<-list()
for (i in 1:length(methods)) {
  c<-methods[i]
  accu.LZ[[c]] <- vegan::specaccum(animal.LZ[metadata.LZ$method == 
                                              methods[i], ])
}
##plot
tidy_specaccum <- function(x) {
  table1<-data.frame()
  method<-names(x)

  for (i in 1:length(names(x))) {
    richness = x[[i]][["richness"]]
    sd = x[[i]]$sd
    table<-data.frame(
      site = x[[i]]$sites,
      richness = richness,
      sd = sd,
      ymin=(richness - 2*sd),
      ymax =( richness + 2*sd),
      methods=method[i])
    table1<-rbind(table1,table)
  }
  return(table1)
}

accu_.LZ2<-tidy_specaccum(accu.LZ)
library(ggplot2)
p<-ggplot(data = accu_.LZ2[!grepl("all",accu_.LZ2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
  geom_line(size=0.7)  +
  labs(x="Number of plots",y="Richness")+
  scale_color_discrete(name = "methods") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3",
                   "#EE8110","#E00B00" ,"#111FDD","#A99EFF"))
names(cbPalette)<-methods

p + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) 
#LV
animal.LV<-animal2[grepl("LV",row.names(animal2)),]
animal.LV<-animal.LV[!grepl("LV7",row.names(animal.LV)),]
animal.LV<-animal.LV[!grepl("B",row.names(animal.LV)),]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
rownames(metadata) <- metadata$site
metadata<-metadata[,1:2]
metadata.LV<-metadata[match(row.names(animal.LV),metadata$site),]
##separate the site
animal.LV<-animal.LV[,colSums(animal.LV)!=0]

animal.LV<-as.data.frame(t(animal.LV))
GSMc40A<-gsub("LZ","LV",GSMc40A)
GSMc40B<-gsub("LZ","LV",GSMc40B)

A40<-animal.LV[,grepl(GSMc40A,names(animal.LV))]
names(A40)<-paste0(names(A40),"40A")

B40<-animal.LV[,grepl(GSMc40B,names(animal.LV))]
names(B40)<-paste0(names(B40),"40B")

metadata.LV<-rbind(
  data.frame(method="GSMc_40A",site=names(A40)),
  data.frame(method="GSMc_40B",site=names(B40)),
  metadata.LV)

try<-function(x,y){
  merge(x,y,by="OTU",all=T)
}
A40$OTU<-row.names(A40)
B40$OTU<-row.names(B40)

animal.LV$OTU<-row.names(animal.LV)
a<-list(A40,B40,animal.LV)
animal.LV<-Reduce(try,a)
row.names(animal.LV)<-animal.LV$OTU
animal.LV<-animal.LV[,-1]
animal.LV<-as.data.frame(t(animal.LV))

methods<-unique(metadata.LV$method)
accu.LV<-list()
for (i in 1:length(methods)) {
  c<-methods[i]
  accu.LV[[c]] <- vegan::specaccum(animal.LV[metadata.LV$method == 
                                               methods[i], ])
}

##plot
accu_.LV2<-tidy_specaccum(  accu.LV)
library(ggplot2)
p<-ggplot(data = accu_.LV2, aes(x = site, y = richness, group = methods,col = methods)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
  geom_line(size=0.7)  +
  labs(x="Number of plots",y="Richness")+
  scale_color_discrete(name = "methods") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

p + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) 
#LW
animal.LW<-animal2[grepl("LW",row.names(animal2)),]
animal.LW<-animal.LW[!grepl("LW7",row.names(animal.LW)),]
animal.LW<-animal.LW[!grepl("B",row.names(animal.LW)),]

metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
rownames(metadata) <- metadata$site
metadata<-metadata[,1:2]
metadata.LW<-metadata[match(row.names(animal.LW),metadata$site),]
##separate the site
animal.LW<-animal.LW[,colSums(animal.LW)!=0]

animal.LW<-as.data.frame(t(animal.LW))
GSMc40A<-gsub("LV","LW",GSMc40A)
GSMc40B<-gsub("LV","LW",GSMc40B)

A40<-animal.LW[,grepl(GSMc40A,names(animal.LW))]
names(A40)<-paste0(names(A40),"40A")

B40<-animal.LW[,grepl(GSMc40B,names(animal.LW))]
names(B40)<-paste0(names(B40),"40B")

metadata.LW<-rbind(
  data.frame(method="GSMc_40A",site=names(A40)),
  data.frame(method="GSMc_40B",site=names(B40)),
  metadata.LW)

try<-function(x,y){
  merge(x,y,by="OTU",all=T)
}
A40$OTU<-row.names(A40)
B40$OTU<-row.names(B40)

animal.LW$OTU<-row.names(animal.LW)
a<-list(A40,B40,animal.LW)
animal.LW<-Reduce(try,a)
row.names(animal.LW)<-animal.LW$OTU
animal.LW<-animal.LW[,-1]
animal.LW<-as.data.frame(t(animal.LW))

methods<-unique(metadata.LW$method)
accu.LW<-list()
for (i in 1:length(methods)) {
  c<-methods[i]
  accu.LW[[c]] <- vegan::specaccum(animal.LW[metadata.LW$method == 
                                               methods[i], ])
}

accu.LW[["all"]]<-vegan::specaccum(animal.LW)
##plot
accu_.LW2<-tidy_specaccum(accu.LW)
library(ggplot2)
p<-ggplot(data = accu_.LW2[!grepl("all",accu_.LW2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
  geom_line(size=0.7)  +
  labs(x="Number of plots",y="Richness")+
  scale_color_discrete(name = "methods") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

label<-data.frame(xlab=tapply(accu_.LW2$site, accu_.LW2$methods,max),
                  ylab=tapply(accu_.LW2$richness,accu_.LW2$methods,  max),
                  methods=names(tapply(accu_.LW2$richness,accu_.LW2$methods,  max)))

p + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) 

##fungal accumulation curve
  fungi2<-as.data.frame(t(fungi))
  fungi2<-fungi2[!grepl("DNA|soil",row.names(fungi2)),]
  library(vegan)
  fungi.LZ<-fungi2[grepl("LZ",row.names(fungi2)),]
  fungi.LZ<-fungi.LZ[!grepl("LZ7",row.names(fungi.LZ)),]
  fungi.LZ<-fungi.LZ[!grepl("B",row.names(fungi.LZ)),]
  metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
  rownames(metadata) <- metadata$site
  metadata<-metadata[,1:2]
  metadata.LZ<-metadata[match(row.names(fungi.LZ),metadata$site),]
  ##separate the site
  fungi.LZ<-fungi.LZ[,colSums(fungi.LZ)!=0]
  
  fungi.LZ<-as.data.frame(t(fungi.LZ))
  GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
  GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
  A40<-fungi.LZ[,grepl(GSMc40A,names(fungi.LZ))]
  names(A40)<-paste0(names(A40),"40A")
  
  B40<-fungi.LZ[,grepl(GSMc40B,names(fungi.LZ))]
  names(B40)<-paste0(names(B40),"40B")
  
  metadata.LZ<-rbind(
    data.frame(method="GSMc_40A",site=names(A40)),
    data.frame(method="GSMc_40B",site=names(B40)),
    metadata.LZ)
  
  try<-function(x,y){
    merge(x,y,by="OTU",all=T)
  }
  A40$OTU<-row.names(A40)
  B40$OTU<-row.names(B40)
  
  fungi.LZ$OTU<-row.names(fungi.LZ)
  a<-list(A40,B40,fungi.LZ)
  fungi.LZ<-Reduce(try,a)
  row.names(fungi.LZ)<-fungi.LZ$OTU
  fungi.LZ<-fungi.LZ[,-1]
  fungi.LZ<-as.data.frame(t(fungi.LZ))
  
  methods<-unique(metadata.LZ$method)
  accu.LZ<-list()
  for (i in 1:length(methods)) {
    c<-methods[i]
    accu.LZ[[c]] <- vegan::specaccum(fungi.LZ[metadata.LZ$method == 
                                                methods[i], ])
  }
  ##plot
  accu_.LZ2<-tidy_specaccum(  accu.LZ)
  library(ggplot2)
  p<-ggplot(data = accu_.LZ2[!grepl("all",accu_.LZ2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
    geom_line(size=0.7)  +
    labs(x="Number of plots",y="Richness")+
    scale_color_discrete(name = "methods") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  p + scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) 

  #LV
  fungi.LV<-fungi2[grepl("LV",row.names(fungi2)),]
  fungi.LV<-fungi.LV[!grepl("LV7",row.names(fungi.LV)),]
  fungi.LV<-fungi.LV[!grepl("B",row.names(fungi.LV)),]
  metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
  rownames(metadata) <- metadata$site
  metadata<-metadata[,1:2]
  metadata.LV<-metadata[match(row.names(fungi.LV),metadata$site),]
  ##separate the site
  fungi.LV<-fungi.LV[,colSums(fungi.LV)!=0]
  
  fungi.LV<-as.data.frame(t(fungi.LV))
  GSMc40A<-gsub("LZ","LV",GSMc40A)
  GSMc40B<-gsub("LZ","LV",GSMc40B)
  
  A40<-fungi.LV[,grepl(GSMc40A,names(fungi.LV))]
  names(A40)<-paste0(names(A40),"40A")
  
  B40<-fungi.LV[,grepl(GSMc40B,names(fungi.LV))]
  names(B40)<-paste0(names(B40),"40B")
  
  metadata.LV<-rbind(
    data.frame(method="GSMc_40A",site=names(A40)),
    data.frame(method="GSMc_40B",site=names(B40)),
    metadata.LV)
  
  try<-function(x,y){
    merge(x,y,by="OTU",all=T)
  }
  A40$OTU<-row.names(A40)
  B40$OTU<-row.names(B40)
  
  fungi.LV$OTU<-row.names(fungi.LV)
  a<-list(A40,B40,fungi.LV)
  fungi.LV<-Reduce(try,a)
  row.names(fungi.LV)<-fungi.LV$OTU
  fungi.LV<-fungi.LV[,-1]
  fungi.LV<-as.data.frame(t(fungi.LV))
  methods<-unique(metadata.LV$method)
  accu.LV<-list()
  for (i in 1:length(methods)) {
    c<-methods[i]
    accu.LV[[c]] <- vegan::specaccum(fungi.LV[metadata.LV$method == 
                                                methods[i], ])
  }
  ##plot
  accu_.LV2<-tidy_specaccum(  accu.LV)
  p<-ggplot(data = accu_.LV2[!grepl("all",accu_.LV2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
    geom_line(size=0.7)  +
    labs(x="Number of plots",y="Richness")+
    scale_color_discrete(name = "methods") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
p + scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) 
  #LW
  fungi.LW<-fungi2[grepl("LW",row.names(fungi2)),]
  fungi.LW<-fungi.LW[!grepl("LW7",row.names(fungi.LW)),]
  fungi.LW<-fungi.LW[!grepl("B",row.names(fungi.LW)),]
  metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
  rownames(metadata) <- metadata$site
  metadata<-metadata[,1:2]
  metadata.LW<-metadata[match(row.names(fungi.LW),metadata$site),]
  ##separate the site
  fungi.LW<-fungi.LW[,colSums(fungi.LW)!=0]
  fungi.LW<-as.data.frame(t(fungi.LW))
  GSMc40A<-gsub("LV","LW",GSMc40A)
  GSMc40B<-gsub("LV","LW",GSMc40B)
  
  A40<-fungi.LW[,grepl(GSMc40A,names(fungi.LW))]
  names(A40)<-paste0(names(A40),"40A")
  
  B40<-fungi.LW[,grepl(GSMc40B,names(fungi.LW))]
  names(B40)<-paste0(names(B40),"40B")
  
  metadata.LW<-rbind(
    data.frame(method="GSMc_40A",site=names(A40)),
    data.frame(method="GSMc_40B",site=names(B40)),
    metadata.LW)
  
  try<-function(x,y){
    merge(x,y,by="OTU",all=T)
  }
  A40$OTU<-row.names(A40)
  B40$OTU<-row.names(B40)
  
  fungi.LW$OTU<-row.names(fungi.LW)
  a<-list(A40,B40,fungi.LW)
  fungi.LW<-Reduce(try,a)
  row.names(fungi.LW)<-fungi.LW$OTU
  fungi.LW<-fungi.LW[,-1]
  fungi.LW<-as.data.frame(t(fungi.LW))
  methods<-unique(metadata.LW$method)
  accu.LW<-list()
  for (i in 1:length(methods)) {
    c<-methods[i]
    accu.LW[[c]] <- vegan::specaccum(fungi.LW[metadata.LW$method == 
                                                methods[i], ])
  }
  
  ##plot
  accu_.LW2<-tidy_specaccum(  accu.LW)
  library(ggplot2)
  p<-ggplot(data = accu_.LW2[!grepl("all",accu_.LW2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
    geom_line(size=0.7)  +
    labs(x="Number of plots",y="Richness")+
    scale_color_discrete(name = "methods") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  p + scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) 
 ##bacteria accumulation curve
  bacteria2<-as.data.frame(t(bacteria))
  bacteria2<-bacteria2[!grepl("DNA|soil",row.names(bacteria2)),]
  bacteria.LZ<-bacteria2[grepl("LZ",row.names(bacteria2)),]
  bacteria.LZ<-bacteria.LZ[!grepl("LZ7",row.names(bacteria.LZ)),]
  bacteria.LZ<-bacteria.LZ[!grepl("B",row.names(bacteria.LZ)),]
  
  metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
  rownames(metadata) <- metadata$site
  metadata<-metadata[,1:2]
  metadata.LZ<-metadata[match(row.names(bacteria.LZ),metadata$site),]
  ##separate the site
  bacteria.LZ<-bacteria.LZ[,colSums(bacteria.LZ)!=0]
  bacteria.LZ<-as.data.frame(t(bacteria.LZ))
  GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
  GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
  A40<-bacteria.LZ[,grepl(GSMc40A,names(bacteria.LZ))]
  names(A40)<-paste0(names(A40),"40A")
  
  B40<-bacteria.LZ[,grepl(GSMc40B,names(bacteria.LZ))]
  names(B40)<-paste0(names(B40),"40B")
  
  metadata.LZ<-rbind(
    data.frame(method="GSMc_40A",site=names(A40)),
    data.frame(method="GSMc_40B",site=names(B40)),
    metadata.LZ)
  
  try<-function(x,y){
    merge(x,y,by="OTU",all=T)
  }
  A40$OTU<-row.names(A40)
  B40$OTU<-row.names(B40)
  
  bacteria.LZ$OTU<-row.names(bacteria.LZ)
  a<-list(A40,B40,bacteria.LZ)
  bacteria.LZ<-Reduce(try,a)
  row.names(bacteria.LZ)<-bacteria.LZ$OTU
  bacteria.LZ<-bacteria.LZ[,-1]
  bacteria.LZ<-as.data.frame(t(bacteria.LZ))
  methods<-unique(metadata.LZ$method)
  accu.LZ<-list()
  for (i in 1:length(methods)) {
    c<-methods[i]
    accu.LZ[[c]] <- vegan::specaccum(bacteria.LZ[metadata.LZ$method ==methods[i], ])
  }
  
  ##plot
  accu_.LZ2<-tidy_specaccum(  accu.LZ)
  library(ggplot2)
  p<-ggplot(data = accu_.LZ2[!grepl("all",accu_.LZ2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
    geom_line(size=0.7)  +
    labs(x="Number of plots",y="Richness")+
    scale_color_discrete(name = "methods") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
    
 p + scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) 
   #LV
  library(vegan)
  bacteria.LV<-bacteria2[grepl("LV",row.names(bacteria2)),]
  bacteria.LV<-bacteria.LV[!grepl("LV7",row.names(bacteria.LV)),]
  bacteria.LV<-bacteria.LV[!grepl("LVB5",row.names(bacteria.LV)),]
  bacteria.LV<-bacteria.LV[!grepl("B",row.names(bacteria.LV)),]
  metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
  rownames(metadata) <- metadata$site
  metadata<-metadata[,1:2]
  metadata.LV<-metadata[match(row.names(bacteria.LV),metadata$site),]
  ##separate the site
  bacteria.LV<-as.data.frame(t(bacteria.LV))
  GSMc40A<-gsub("LZ","LV",GSMc40A)
  GSMc40B<-gsub("LZ","LV",GSMc40B)
  
  A40<-bacteria.LV[,grepl(GSMc40A,names(bacteria.LV))]
  names(A40)<-paste0(names(A40),"40A")
  
  B40<-bacteria.LV[,grepl(GSMc40B,names(bacteria.LV))]
  names(B40)<-paste0(names(B40),"40B")
  
  metadata.LV<-rbind(
    data.frame(method="GSMc_40A",site=names(A40)),
    data.frame(method="GSMc_40B",site=names(B40)),
    metadata.LV)
  
  try<-function(x,y){
    merge(x,y,by="OTU",all=T)
  }
  A40$OTU<-row.names(A40)
  B40$OTU<-row.names(B40)
  
  bacteria.LV$OTU<-row.names(bacteria.LV)
  a<-list(A40,B40,bacteria.LV)
  bacteria.LV<-Reduce(try,a)
  row.names(bacteria.LV)<-bacteria.LV$OTU
  bacteria.LV<-bacteria.LV[,-1]
  bacteria.LV<-as.data.frame(t(bacteria.LV))
  
  methods<-unique(metadata.LV$method)
  
  accu.LV<-list()
  for (i in 1:length(methods)) {
    c<-methods[i]
    accu.LV[[c]] <- vegan::specaccum(bacteria.LV[metadata.LV$method == methods[i], ])
  }

  library(ggplot2)
  accu_.LV2<-tidy_specaccum(  accu.LV)
  p<-ggplot(data = accu_.LV2[!grepl("all",accu_.LV2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
    geom_line(size=0.7)  +
    labs(x="Number of plots",y="Richness")+
    scale_color_discrete(name = "methods") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
p + scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) 
  ##LW
  bacteria.LW<-bacteria2[grepl("LW",row.names(bacteria2)),]
  bacteria.LW<-bacteria.LW[!grepl("LW7",row.names(bacteria.LW)),]
  bacteria.LW<-bacteria.LW[!grepl("B",row.names(bacteria.LW)),]
  
  metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
  rownames(metadata) <- metadata$site
  metadata<-metadata[,1:2]
  metadata.LW<-metadata[match(row.names(bacteria.LW),metadata$site),]
  ##separate the site
  bacteria.LW<-as.data.frame(t(bacteria.LW))
  GSMc40A<-gsub("LV","LW",GSMc40A)
  GSMc40B<-gsub("LV","LW",GSMc40B)
  
  A40<-bacteria.LW[,grepl(GSMc40A,names(bacteria.LW))]
  names(A40)<-paste0(names(A40),"40A")
  
  B40<-bacteria.LW[,grepl(GSMc40B,names(bacteria.LW))]
  names(B40)<-paste0(names(B40),"40B")
  
  metadata.LW<-rbind(
    data.frame(method="GSMc_40A",site=names(A40)),
    data.frame(method="GSMc_40B",site=names(B40)),
    metadata.LW)
  
  try<-function(x,y){
    merge(x,y,by="OTU",all=T)
  }
  A40$OTU<-row.names(A40)
  B40$OTU<-row.names(B40)
  
  bacteria.LW$OTU<-row.names(bacteria.LW)
  a<-list(A40,B40,bacteria.LW)
  bacteria.LW<-Reduce(try,a)
  row.names(bacteria.LW)<-bacteria.LW$OTU
  bacteria.LW<-bacteria.LW[,-1]
  bacteria.LW<-as.data.frame(t(bacteria.LW))
  methods<-unique(metadata.LW$method)
  
  accu.LW<-list()
  for (i in 1:length(methods)) {
    c<-methods[i]
    accu.LW[[c]] <- vegan::specaccum(bacteria.LW[metadata.LW$method == methods[i], ])
  }
  

  accu_.LW2<-tidy_specaccum(  accu.LW)
  p<-ggplot(data = accu_.LW2[!grepl("all",accu_.LW2$methods),], aes(x = site, y = richness, group = methods,col = methods)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,col = NULL, fill = methods), alpha = 0.05) +
    geom_line(size=0.7)  +
    labs(x="Number of plots",y="Richness")+
    scale_color_discrete(name = "methods") +
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
p + scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) 

#####Marginal estimates
###mean within-site richness and Shannon diversity index
###richness
##animal
animal.normal.m<-read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
##normal
names(animal.normal.m)[names(animal.normal.m)=="method"]<-"method.x"
mod01<-lmer(mean~  method.x + (1|site),data = animal.normal.m)
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
fungi.normal.m<-read.csv("fungi.totalrichness.csv",header = TRUE,row.names = 1,sep = ",")
##normal
names(fungi.normal.m)[names(fungi.normal.m)=="method"]<-"method.x"
mod01<-lmer(mean~  method.x + (1|site),data = fungi.normal.m)
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
bacteria.normal.m<-read.csv("bacteria.totalrichness.csv",header = TRUE,row.names = 1,sep = ",")
##normal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
names(bacteria.normal.m)[names(bacteria.normal.m)=="method"]<-"method.x"
mod01<-lmer(mean~  method.x + (1|site),data = bacteria.normal.m)
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

###shannon diversity index
##animal
animal.normal.m<-read.csv("unpooled.animal.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m<-animal.normal.m[,c(1,6,7)]
mod01<-lm(log(shannon)~  method + site,data = animal.normal.m)
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
fungi.normal.m<-read.csv("unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m<-fungi.normal.m[,c(1,6,7)]
mod01<-lm(log(shannon)~  method + site,data = fungi.normal.m)
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
bacteria.normal.m<-read.csv("unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m<-bacteria.normal.m[,c(1,6,7)]
mod01<-lm(log(shannon)~  method + site,data = bacteria.normal.m)
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

###individual samples richness comparison
##animal
animal.normal.m<-read.csv("animal.richness.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m$site<-row.names(animal.normal.m)
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
animal.normal.m<-merge(animal.normal.m,metadata,by="site",all=T)
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$estimate),]
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$method),]
animal.normal.m$site2<-substr(animal.normal.m$site,1,2)

animal.normal.m<-animal.normal.m[,c(1,2,5,10)]
animal.normal.m<-animal.normal.m[!grepl("B",animal.normal.m$site),]
###GSMc40
select<-read.csv("sheet2_for_pooled.csv")
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
names(animal.normal.m)[names(animal.normal.m)=="method"]<-"method.x"
names(animal.normal.m)[names(animal.normal.m)=="estimate"]<-"estimate.x"
animal.normal.m<-animal.normal.m[!grepl("deep_SUCC",animal.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = animal.normal.m)
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
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "Zobel","GSMc","GSMc40A","GSMc40B","DarkDiv","SUCC","LUCAS","Deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##bacteria
bacteria.normal.m<-read.csv("bacteria.richness.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m$site<-row.names(bacteria.normal.m)
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
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
names(bacteria.normal.m)[names(bacteria.normal.m)=="method"]<-"method.x"
names(bacteria.normal.m)[names(bacteria.normal.m)=="estimate"]<-"estimate.x"
bacteria.normal.m<-bacteria.normal.m[!grepl("deep_SUCC",bacteria.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = bacteria.normal.m)
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
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "Zobel","GSMc","GSMc40A","GSMc40B","DarkDiv","SUCC","LUCAS","Deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

##fungi
fungi.normal.m<-read.csv("fungi.richness.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m$site<-row.names(fungi.normal.m)
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
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
names(fungi.normal.m)[names(fungi.normal.m)=="method"]<-"method.x"
names(fungi.normal.m)[names(fungi.normal.m)=="estimate"]<-"estimate.x"
fungi.normal.m<-fungi.normal.m[!grepl("deep_SUCC",fungi.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = fungi.normal.m)
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
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "Zobel","GSMc","GSMc40A","GSMc40B","DarkDiv","SUCC","LUCAS","Deep")))+ 
  scale_y_continuous(breaks=c(0,25,50,75,100,150))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###shannon diversity index
##animal
animal.normal.m<-read.csv("animal.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m$site<-row.names(animal.normal.m)
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
animal.normal.m<-merge(animal.normal.m,metadata,by="site",all=T)
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$estimate),]
animal.normal.m <- animal.normal.m[!is.na(animal.normal.m$method),]
animal.normal.m$site2<-substr(animal.normal.m$site,1,2)

animal.normal.m<-animal.normal.m[,c(1,2,5,10)]
animal.normal.m<-animal.normal.m[!grepl("B",animal.normal.m$site),]
###GSMc40
select<-read.csv("sheet2_for_pooled.csv")
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
names(animal.normal.m)[names(animal.normal.m)=="method"]<-"method.x"
names(animal.normal.m)[names(animal.normal.m)=="estimate"]<-"estimate.x"
animal.normal.m<-animal.normal.m[!grepl("deep_SUCC",animal.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = animal.normal.m)
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
##bacteria
bacteria.normal.m<-read.csv("bacteria.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m$site<-row.names(bacteria.normal.m)
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
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
names(bacteria.normal.m)[names(bacteria.normal.m)=="method"]<-"method.x"
names(bacteria.normal.m)[names(bacteria.normal.m)=="estimate"]<-"estimate.x"
bacteria.normal.m<-bacteria.normal.m[!grepl("deep_SUCC",bacteria.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = bacteria.normal.m)
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

##fungi
fungi.normal.m<-read.csv("fungi.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m$site<-row.names(fungi.normal.m)
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
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

names(fungi.normal.m)[names(fungi.normal.m)=="method"]<-"method.x"
names(fungi.normal.m)[names(fungi.normal.m)=="estimate"]<-"estimate.x"
fungi.normal.m<-fungi.normal.m[!grepl("deep_SUCC",fungi.normal.m$method.x),]
mod01<-lmer(estimate.x~  method.x + (1|site2),data = fungi.normal.m)
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

###ranking accumulation curve
#animal
normal<-read.table("animal.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata2<-metadata
metadata$site<-paste0(metadata$site,"DNA")
metadata2$site<-paste0(metadata2$site,"soil")
metadata$method<-paste0(metadata$method,"DNA")
metadata2$method<-paste0(metadata2$method,"soil")
metadata<-rbind(metadata,metadata2)
metadata.normal2<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
metadata<-metadata[,c(1,5)]
metadata<-rbind(metadata,metadata.normal2)

normal.fre2<-merge(normal.fre,metadata,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$COI_Otu1),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method!="Meastre",]
normal.fre2<-aggregate(normal.fre2[,2:(ncol(normal.fre2)-2)],by=list(method=normal.fre2$method,site2=normal.fre2$site2),sum)
###for GSMc40
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-normal[,grepl(GSMc40A,names(normal))]
names(A40)<-paste0(names(A40),"40A")

B40<-normal[,grepl(GSMc40B,names(normal))]
names(B40)<-paste0(names(B40),"40B")

metadata.normal<-rbind(
  data.frame(method="GSMc_40A",site=names(A40)),
  data.frame(method="GSMc_40B",site=names(B40)),
  metadata)

try<-function(x,y){
  merge(x,y,by="OTU",all=T)
}
A40$OTU<-row.names(A40)
B40$OTU<-row.names(B40)

a<-list(A40,B40)
animal.LZ<-Reduce(try,a)
row.names(animal.LZ)<-animal.LZ$OTU
animal.LZ<-animal.LZ[,-1]
animal.LZ<-as.data.frame(t(animal.LZ))
animal.LZ$site<-row.names(animal.LZ)

animal.LZ2<-merge(animal.LZ,metadata.normal,by="site")
animal.LZ2$site2<-substr(animal.LZ2$site,1,2)
animal.LZ2<-aggregate(animal.LZ2[,2:(ncol(animal.LZ2)-2)],by=list(method=animal.LZ2$method,site2=animal.LZ2$site2),sum)

animal.LZ2<-as.data.frame(t(animal.LZ2))
normal.fre2<-as.data.frame(t(normal.fre2))
normal.fre2$me<-row.names(normal.fre2)
animal.LZ2$me<-row.names(animal.LZ2)
a<-merge(normal.fre2,animal.LZ2,by="me",all=T)
row.names(a)<-a$me
a<-a[,-1]
a<-as.data.frame(t(a))
normal.fre2<-a
###
deep_SUCC<-normal.fre2[grepl("Deep|SUCC",normal.fre2$method),]
deep_SUCC<-deep_SUCC[!grepl("DNA|soil",deep_SUCC$method),]
deep_SUCC<-deep_SUCC[!grepl("deep_SUCC",deep_SUCC$method),]
deep_SUCC[,1:(ncol(deep_SUCC)-2)]<-apply(deep_SUCC[,1:(ncol(deep_SUCC)-2)],2,as.numeric)
deep_SUCC<-aggregate(deep_SUCC[,1:(ncol(deep_SUCC)-2)],by=list(site2=deep_SUCC$site2),sum)
deep_SUCC$method<-"deep_SUCC"
normal.fre2<-rbind(normal.fre2,deep_SUCC)

##LZ
LZ <- normal.fre2[grepl("LZ",normal.fre2$site2),]
row.names(LZ)<-LZ$method
LZ <-LZ[,1:(ncol(LZ)-2)]
## cumulative abundance
LZ <- as.data.frame(t(LZ))
LZ<-as.data.frame(apply(LZ,2,as.numeric))
for (i in 1:length(colnames(LZ))){
  LZ[,i] <- LZ[,i][order(LZ[,i], decreasing = T)]
  
}

method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LZ <- LZ[rowMeans(LZ)!= 0,]
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40A"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40A"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40B"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40B"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="Maestre"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"Maestre"
LZ.Maestre$index<-row.names(LZ.Maestre)


LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDiv"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDiv"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCC"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCC"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCAS"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCAS"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="Deep"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deep"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCC"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCC"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="Zobel"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"Zobel"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
                  )
all_data$index<-as.numeric(all_data$index)
library(ggplot2)
cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
library(dplyr)
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

###for the pooled
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc_62soil"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc_62soil"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40Asoil"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40Asoil"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40Bsoil"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40Bsoil"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="Maestresoil"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"Maestresoil"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDivsoil"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDivsoil"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCCsoil"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCCsoil"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCASsoil"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCASsoil"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="deepsoil"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deepsoil"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCCsoil"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCCsoil"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.MDB15<-as.data.frame(LZ[,names(LZ)=="MDB15soil"])
LZ.MDB15<- as.data.frame(LZ.MDB15[rowMeans(LZ.MDB15)!= 0,])
LZ.MDB15$mean<-round(rowMeans(LZ.MDB15))
LZ.MDB15$cumulate<-LZ.MDB15$mean
for (i in 2:length(LZ.MDB15$cumulate))
  LZ.MDB15$cumulate[i] <- LZ.MDB15$mean[i]+LZ.MDB15$cumulate[i-1]
LZ.MDB15$method<-"MDB15soil"
LZ.MDB15$index<-row.names(LZ.MDB15)

LZ.MDB5<-as.data.frame(LZ[,names(LZ)=="MDB5soil"])
LZ.MDB5<- as.data.frame(LZ.MDB5[rowMeans(LZ.MDB5)!= 0,])
LZ.MDB5$mean<-round(rowMeans(LZ.MDB5))
LZ.MDB5$cumulate<-LZ.MDB5$mean
for (i in 2:length(LZ.MDB5$cumulate))
  LZ.MDB5$cumulate[i] <- LZ.MDB5$mean[i]+LZ.MDB5$cumulate[i-1]
LZ.MDB5$method<-"MDB5soil"
LZ.MDB5$index<-row.names(LZ.MDB5)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="Zobelsoil"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"Zobelsoil"
LZ.Zobel$index<-row.names(LZ.Zobel)
#

all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)
cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc_62DNA"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc_62DNA"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40ADNA"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40ADNA"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40BDNA"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40BDNA"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="MaestreDNA"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"MaestreDNA"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDivDNA"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDivDNA"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCCDNA"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCCDNA"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCASDNA"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCASDNA"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="deepDNA"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deepDNA"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCCDNA"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCCDNA"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.MDB15<-as.data.frame(LZ[,names(LZ)=="MDB15DNA"])
LZ.MDB15<- as.data.frame(LZ.MDB15[rowMeans(LZ.MDB15)!= 0,])
LZ.MDB15$mean<-round(rowMeans(LZ.MDB15))
LZ.MDB15$cumulate<-LZ.MDB15$mean
for (i in 2:length(LZ.MDB15$cumulate))
  LZ.MDB15$cumulate[i] <- LZ.MDB15$mean[i]+LZ.MDB15$cumulate[i-1]
LZ.MDB15$method<-"MDB15DNA"
LZ.MDB15$index<-row.names(LZ.MDB15)

LZ.MDB5<-as.data.frame(LZ[,names(LZ)=="MDB5DNA"])
LZ.MDB5<- as.data.frame(LZ.MDB5[rowMeans(LZ.MDB5)!= 0,])
LZ.MDB5$mean<-round(rowMeans(LZ.MDB5))
LZ.MDB5$cumulate<-LZ.MDB5$mean
for (i in 2:length(LZ.MDB5$cumulate))
  LZ.MDB5$cumulate[i] <- LZ.MDB5$mean[i]+LZ.MDB5$cumulate[i-1]
LZ.MDB5$method<-"MDB5DNA"
LZ.MDB5$index<-row.names(LZ.MDB5)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="ZobelDNA"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"ZobelDNA"
LZ.Zobel$index<-row.names(LZ.Zobel)
#

all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)
cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##LV
LV <- normal.fre2[grepl("LV",normal.fre2$site2),]
row.names(LV)<-LV$method
LV <-LV[,1:(ncol(LV)-2)]
## cumulative abundance
LV <- as.data.frame(t(LV))
LV<-as.data.frame(apply(LV,2,as.numeric))
for (i in 1:length(colnames(LV))){
  LV[,i] <- LV[,i][order(LV[,i], decreasing = T)]
  
}

method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LV <- LV[rowMeans(LV)!= 0,]
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40A"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40A"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40B"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40B"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="Maestre"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"Maestre"
LV.Maestre$index<-row.names(LV.Maestre)

LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDiv"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDiv"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCC"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCC"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCAS"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCAS"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="Deep"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deep"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCC"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCC"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.Zobel<-as.data.frame(LV[,names(LV)=="Zobel"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"Zobel"
LV.Zobel$index<-row.names(LV.Zobel)
#
all_data<-rbind(  LV.Zobel[,3:5],
                  LV.LUCAS[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

###for the pooled
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc_62soil"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc_62soil"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40Asoil"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40Asoil"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40Bsoil"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40Bsoil"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="Maestresoil"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"Maestresoil"
LV.Maestre$index<-row.names(LV.Maestre)

LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDivsoil"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDivsoil"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCCsoil"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCCsoil"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCASsoil"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCASsoil"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="deepsoil"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deepsoil"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCCsoil"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCCsoil"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.MDB15<-as.data.frame(LV[,names(LV)=="MDB15soil"])
LV.MDB15<- as.data.frame(LV.MDB15[rowMeans(LV.MDB15)!= 0,])
LV.MDB15$mean<-round(rowMeans(LV.MDB15))
LV.MDB15$cumulate<-LV.MDB15$mean
for (i in 2:length(LV.MDB15$cumulate))
  LV.MDB15$cumulate[i] <- LV.MDB15$mean[i]+LV.MDB15$cumulate[i-1]
LV.MDB15$method<-"MDB15soil"
LV.MDB15$index<-row.names(LV.MDB15)

LV.MDB5<-as.data.frame(LV[,names(LV)=="MDB5soil"])
LV.MDB5<- as.data.frame(LV.MDB5[rowMeans(LV.MDB5)!= 0,])
LV.MDB5$mean<-round(rowMeans(LV.MDB5))
LV.MDB5$cumulate<-LV.MDB5$mean
for (i in 2:length(LV.MDB5$cumulate))
  LV.MDB5$cumulate[i] <- LV.MDB5$mean[i]+LV.MDB5$cumulate[i-1]
LV.MDB5$method<-"MDB5soil"
LV.MDB5$index<-row.names(LV.MDB5)

LV.Zobel<-as.data.frame(LV[,names(LV)=="Zobelsoil"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"Zobelsoil"
LV.Zobel$index<-row.names(LV.Zobel)
#

all_data<-rbind(  LV.Zobel[,3:5],
                  LV.LUCAS[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc_62DNA"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc_62DNA"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40ADNA"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40ADNA"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40BDNA"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40BDNA"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="MaestreDNA"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"MaestreDNA"
LV.Maestre$index<-row.names(LV.Maestre)

LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDivDNA"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDivDNA"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCCDNA"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCCDNA"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCASDNA"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCASDNA"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="deepDNA"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deepDNA"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCCDNA"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCCDNA"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.MDB15<-as.data.frame(LV[,names(LV)=="MDB15DNA"])
LV.MDB15<- as.data.frame(LV.MDB15[rowMeans(LV.MDB15)!= 0,])
LV.MDB15$mean<-round(rowMeans(LV.MDB15))
LV.MDB15$cumulate<-LV.MDB15$mean
for (i in 2:length(LV.MDB15$cumulate))
  LV.MDB15$cumulate[i] <- LV.MDB15$mean[i]+LV.MDB15$cumulate[i-1]
LV.MDB15$method<-"MDB15DNA"
LV.MDB15$index<-row.names(LV.MDB15)

LV.MDB5<-as.data.frame(LV[,names(LV)=="MDB5DNA"])
LV.MDB5<- as.data.frame(LV.MDB5[rowMeans(LV.MDB5)!= 0,])
LV.MDB5$mean<-round(rowMeans(LV.MDB5))
LV.MDB5$cumulate<-LV.MDB5$mean
for (i in 2:length(LV.MDB5$cumulate))
  LV.MDB5$cumulate[i] <- LV.MDB5$mean[i]+LV.MDB5$cumulate[i-1]
LV.MDB5$method<-"MDB5DNA"
LV.MDB5$index<-row.names(LV.MDB5)

LV.Zobel<-as.data.frame(LV[,names(LV)=="ZobelDNA"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"ZobelDNA"
LV.Zobel$index<-row.names(LV.Zobel)
#

all_data<-rbind(  LV.Zobel[,3:5],
                  LV.LUCAS[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##LW
LW <- normal.fre2[grepl("LW",normal.fre2$site2),]
row.names(LW)<-LW$method
LW <-LW[,1:(ncol(LW)-2)]
## cumulative abundance
LW <- as.data.frame(t(LW))
LW<-as.data.frame(apply(LW,2,as.numeric))
for (i in 1:length(colnames(LW))){
  LW[,i] <- LW[,i][order(LW[,i], decreasing = T)]
  
}

method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LW <- LW[rowMeans(LW)!= 0,]
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40A"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40A"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40B"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40B"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDiv"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDiv"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCC"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCC"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCAS"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCAS"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.Zobel<-as.data.frame(LW[,names(LW)=="Zobel"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"Zobel"
LW.Zobel$index<-row.names(LW.Zobel)
#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)



cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

###for the pooled
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc_62soil"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc_62soil"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40Asoil"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40Asoil"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40Bsoil"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40Bsoil"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.Maestre<-as.data.frame(LW[,names(LW)=="Maestresoil"])
LW.Maestre<- as.data.frame(LW.Maestre[rowMeans(LW.Maestre)!= 0,])
LW.Maestre$mean<-round(rowMeans(LW.Maestre))
LW.Maestre$cumulate<-LW.Maestre$mean
for (i in 2:length(LW.Maestre$cumulate))
  LW.Maestre$cumulate[i] <- LW.Maestre$mean[i]+LW.Maestre$cumulate[i-1]
LW.Maestre$method<-"Maestresoil"
LW.Maestre$index<-row.names(LW.Maestre)


LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDivsoil"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDivsoil"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCCsoil"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCCsoil"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCASsoil"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCASsoil"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.deep<-as.data.frame(LW[,names(LW)=="deepsoil"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])

LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deepsoil"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCCsoil"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCCsoil"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.MDB15<-as.data.frame(LW[,names(LW)=="MDB15soil"])
LW.MDB15<- as.data.frame(LW.MDB15[rowMeans(LW.MDB15)!= 0,])
LW.MDB15$mean<-round(rowMeans(LW.MDB15))
LW.MDB15$cumulate<-LW.MDB15$mean
for (i in 2:length(LW.MDB15$cumulate))
  LW.MDB15$cumulate[i] <- LW.MDB15$mean[i]+LW.MDB15$cumulate[i-1]
LW.MDB15$method<-"MDB15soil"
LW.MDB15$index<-row.names(LW.MDB15)

LW.MDB5<-as.data.frame(LW[,names(LW)=="MDB5soil"])
LW.MDB5<- as.data.frame(LW.MDB5[rowMeans(LW.MDB5)!= 0,])
LW.MDB5$mean<-round(rowMeans(LW.MDB5))
LW.MDB5$cumulate<-LW.MDB5$mean
for (i in 2:length(LW.MDB5$cumulate))
  LW.MDB5$cumulate[i] <- LW.MDB5$mean[i]+LW.MDB5$cumulate[i-1]
LW.MDB5$method<-"MDB5soil"
LW.MDB5$index<-row.names(LW.MDB5)

LW.Zobel<-as.data.frame(LW[,names(LW)=="Zobelsoil"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"Zobelsoil"
LW.Zobel$index<-row.names(LW.Zobel)
#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc_62DNA"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc_62DNA"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40ADNA"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40ADNA"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40BDNA"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40BDNA"
LW.GSMc40B$index<-row.names(LW.GSMc40B)



LW.Maestre<-as.data.frame(LW[,names(LW)=="MaestreDNA"])
LW.Maestre<- as.data.frame(LW.Maestre[rowMeans(LW.Maestre)!= 0,])
LW.Maestre$mean<-round(rowMeans(LW.Maestre))
LW.Maestre$cumulate<-LW.Maestre$mean
for (i in 2:length(LW.Maestre$cumulate))
  LW.Maestre$cumulate[i] <- LW.Maestre$mean[i]+LW.Maestre$cumulate[i-1]
LW.Maestre$method<-"MaestreDNA"
LW.Maestre$index<-row.names(LW.Maestre)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDivDNA"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDivDNA"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCCDNA"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCCDNA"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCASDNA"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCASDNA"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.deep<-as.data.frame(LW[,names(LW)=="deepDNA"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])

LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deepDNA"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCCDNA"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCCDNA"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.MDB15<-as.data.frame(LW[,names(LW)=="MDB15DNA"])
LW.MDB15<- as.data.frame(LW.MDB15[rowMeans(LW.MDB15)!= 0,])
LW.MDB15$mean<-round(rowMeans(LW.MDB15))
LW.MDB15$cumulate<-LW.MDB15$mean
for (i in 2:length(LW.MDB15$cumulate))
  LW.MDB15$cumulate[i] <- LW.MDB15$mean[i]+LW.MDB15$cumulate[i-1]
LW.MDB15$method<-"MDB15DNA"
LW.MDB15$index<-row.names(LW.MDB15)

LW.MDB5<-as.data.frame(LW[,names(LW)=="MDB5DNA"])
LW.MDB5<- as.data.frame(LW.MDB5[rowMeans(LW.MDB5)!= 0,])
LW.MDB5$mean<-round(rowMeans(LW.MDB5))
LW.MDB5$cumulate<-LW.MDB5$mean
for (i in 2:length(LW.MDB5$cumulate))
  LW.MDB5$cumulate[i] <- LW.MDB5$mean[i]+LW.MDB5$cumulate[i-1]
LW.MDB5$method<-"MDB5DNA"
LW.MDB5$index<-row.names(LW.MDB5)

LW.Zobel<-as.data.frame(LW[,names(LW)=="ZobelDNA"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"ZobelDNA"
LW.Zobel$index<-row.names(LW.Zobel)
#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)
library(ggplot2)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

#fungi
normal<-read.table("fungi.rarefy.table2.csv",header = T,row.names = 1,sep = ",")

normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata2<-metadata
metadata$site<-paste0(metadata$site,"DNA")
metadata2$site<-paste0(metadata2$site,"soil")
metadata$method<-paste0(metadata$method,"DNA")
metadata2$method<-paste0(metadata2$method,"soil")
metadata<-rbind(metadata,metadata2)

metadata.normal2<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
metadata<-metadata[,c(1,5)]

metadata<-rbind(metadata,metadata.normal2)

normal.fre2<-merge(normal.fre,metadata,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$'0003b6ab41bde80a3cc58f45791695956d51cb52'),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]

normal.fre2<-normal.fre2[normal.fre2$method!="Meastre",]

normal.fre2<-aggregate(normal.fre2[,2:(ncol(normal.fre2)-2)],by=list(method=normal.fre2$method,site2=normal.fre2$site2),sum)
###for GSMc40
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-normal[,grepl(GSMc40A,names(normal))]
names(A40)<-paste0(names(A40),"40A")

B40<-normal[,grepl(GSMc40B,names(normal))]
names(B40)<-paste0(names(B40),"40B")

metadata.normal<-rbind(
  data.frame(method="GSMc_40A",site=names(A40)),
  data.frame(method="GSMc_40B",site=names(B40)),
  metadata)

try<-function(x,y){
  merge(x,y,by="OTU",all=T)
}
A40$OTU<-row.names(A40)
B40$OTU<-row.names(B40)

a<-list(A40,B40)

fungi.LZ<-Reduce(try,a)
row.names(fungi.LZ)<-fungi.LZ$OTU
fungi.LZ<-fungi.LZ[,-1]
fungi.LZ<-as.data.frame(t(fungi.LZ))
fungi.LZ$site<-row.names(fungi.LZ)

fungi.LZ2<-merge(fungi.LZ,metadata.normal,by="site")
fungi.LZ2$site2<-substr(fungi.LZ2$site,1,2)
fungi.LZ2<-aggregate(fungi.LZ2[,2:(ncol(fungi.LZ2)-2)],by=list(method=fungi.LZ2$method,site2=fungi.LZ2$site2),sum)

fungi.LZ2<-as.data.frame(t(fungi.LZ2))
normal.fre2<-as.data.frame(t(normal.fre2))
normal.fre2$me<-row.names(normal.fre2)
fungi.LZ2$me<-row.names(fungi.LZ2)
a<-merge(normal.fre2,fungi.LZ2,by="me",all=T)
row.names(a)<-a$me
a<-a[,-1]
a<-as.data.frame(t(a))
normal.fre2<-a
###
deep_SUCC<-normal.fre2[grepl("Deep|SUCC",normal.fre2$method),]
deep_SUCC<-deep_SUCC[!grepl("DNA|soil",deep_SUCC$method),]
deep_SUCC<-deep_SUCC[!grepl("deep_SUCC",deep_SUCC$method),]
deep_SUCC[,1:(ncol(deep_SUCC)-2)]<-apply(deep_SUCC[,1:(ncol(deep_SUCC)-2)],2,as.numeric)
deep_SUCC<-aggregate(deep_SUCC[,1:(ncol(deep_SUCC)-2)],by=list(site2=deep_SUCC$site2),sum)
deep_SUCC$method<-"deep_SUCC"
normal.fre2<-rbind(normal.fre2,deep_SUCC)
##LZ
LZ <- normal.fre2[grepl("LZ",normal.fre2$site2),]
row.names(LZ)<-LZ$method
LZ <-LZ[,1:(ncol(LZ)-2)]
## cumulative abundance
LZ <- as.data.frame(t(LZ))
LZ<-as.data.frame(apply(LZ,2,as.numeric))
for (i in 1:length(colnames(LZ))){
  LZ[,i] <- LZ[,i][order(LZ[,i], decreasing = T)]
  
}
method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LZ <- LZ[rowMeans(LZ)!= 0,]
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40A"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40A"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40B"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40B"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="Maestre"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"Maestre"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDiv"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDiv"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCC"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCC"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCAS"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCAS"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="Deep"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])
LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deep"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCC"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCC"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="Zobel"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"Zobel"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
###for the pooled
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc_62soil"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc_62soil"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40Asoil"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40Asoil"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40Bsoil"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40Bsoil"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="Maestresoil"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"Maestresoil"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDivsoil"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDivsoil"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCCsoil"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCCsoil"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCASsoil"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCASsoil"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="deepsoil"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deepsoil"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCCsoil"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCCsoil"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.MDB15<-as.data.frame(LZ[,names(LZ)=="MDB15soil"])
LZ.MDB15<- as.data.frame(LZ.MDB15[rowMeans(LZ.MDB15)!= 0,])
LZ.MDB15$mean<-round(rowMeans(LZ.MDB15))
LZ.MDB15$cumulate<-LZ.MDB15$mean
for (i in 2:length(LZ.MDB15$cumulate))
  LZ.MDB15$cumulate[i] <- LZ.MDB15$mean[i]+LZ.MDB15$cumulate[i-1]
LZ.MDB15$method<-"MDB15soil"
LZ.MDB15$index<-row.names(LZ.MDB15)

LZ.MDB5<-as.data.frame(LZ[,names(LZ)=="MDB5soil"])
LZ.MDB5<- as.data.frame(LZ.MDB5[rowMeans(LZ.MDB5)!= 0,])
LZ.MDB5$mean<-round(rowMeans(LZ.MDB5))
LZ.MDB5$cumulate<-LZ.MDB5$mean
for (i in 2:length(LZ.MDB5$cumulate))
  LZ.MDB5$cumulate[i] <- LZ.MDB5$mean[i]+LZ.MDB5$cumulate[i-1]
LZ.MDB5$method<-"MDB5soil"
LZ.MDB5$index<-row.names(LZ.MDB5)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="Zobelsoil"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"Zobelsoil"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc_62DNA"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc_62DNA"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40ADNA"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40ADNA"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40BDNA"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40BDNA"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="MaestreDNA"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"MaestreDNA"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDivDNA"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDivDNA"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCCDNA"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCCDNA"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCASDNA"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCASDNA"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="deepDNA"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deepDNA"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCCDNA"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCCDNA"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.MDB15<-as.data.frame(LZ[,names(LZ)=="MDB15DNA"])
LZ.MDB15<- as.data.frame(LZ.MDB15[rowMeans(LZ.MDB15)!= 0,])
LZ.MDB15$mean<-round(rowMeans(LZ.MDB15))
LZ.MDB15$cumulate<-LZ.MDB15$mean
for (i in 2:length(LZ.MDB15$cumulate))
  LZ.MDB15$cumulate[i] <- LZ.MDB15$mean[i]+LZ.MDB15$cumulate[i-1]
LZ.MDB15$method<-"MDB15DNA"
LZ.MDB15$index<-row.names(LZ.MDB15)

LZ.MDB5<-as.data.frame(LZ[,names(LZ)=="MDB5DNA"])
LZ.MDB5<- as.data.frame(LZ.MDB5[rowMeans(LZ.MDB5)!= 0,])
LZ.MDB5$mean<-round(rowMeans(LZ.MDB5))
LZ.MDB5$cumulate<-LZ.MDB5$mean
for (i in 2:length(LZ.MDB5$cumulate))
  LZ.MDB5$cumulate[i] <- LZ.MDB5$mean[i]+LZ.MDB5$cumulate[i-1]
LZ.MDB5$method<-"MDB5DNA"
LZ.MDB5$index<-row.names(LZ.MDB5)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="ZobelDNA"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"ZobelDNA"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)
cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##LV
LV <- normal.fre2[grepl("LV",normal.fre2$site2),]
row.names(LV)<-LV$method
LV <-LV[,1:(ncol(LV)-2)]
## cumulative abundance
LV <- as.data.frame(t(LV))
LV<-as.data.frame(apply(LV,2,as.numeric))
for (i in 1:length(colnames(LV))){
  LV[,i] <- LV[,i][order(LV[,i], decreasing = T)]
  
}

method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LV <- LV[rowMeans(LV)!= 0,]
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40A"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40A"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40B"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40B"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="Maestre"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"Maestre"
LV.Maestre$index<-row.names(LV.Maestre)


LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDiv"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDiv"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCC"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCC"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCAS"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCAS"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="Deep"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deep"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCC"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCC"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.Zobel<-as.data.frame(LV[,names(LV)=="Zobel"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"Zobel"
LV.Zobel$index<-row.names(LV.Zobel)
#
all_data<-rbind(  LV.Zobel[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

###for the pooled
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc_62soil"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc_62soil"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40Asoil"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40Asoil"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40Bsoil"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40Bsoil"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="Maestresoil"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"Maestresoil"
LV.Maestre$index<-row.names(LV.Maestre)

LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDivsoil"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDivsoil"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCCsoil"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCCsoil"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCASsoil"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCASsoil"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="deepsoil"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deepsoil"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCCsoil"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCCsoil"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.MDB15<-as.data.frame(LV[,names(LV)=="MDB15soil"])
LV.MDB15<- as.data.frame(LV.MDB15[rowMeans(LV.MDB15)!= 0,])
LV.MDB15$mean<-round(rowMeans(LV.MDB15))
LV.MDB15$cumulate<-LV.MDB15$mean
for (i in 2:length(LV.MDB15$cumulate))
  LV.MDB15$cumulate[i] <- LV.MDB15$mean[i]+LV.MDB15$cumulate[i-1]
LV.MDB15$method<-"MDB15soil"
LV.MDB15$index<-row.names(LV.MDB15)

LV.MDB5<-as.data.frame(LV[,names(LV)=="MDB5soil"])
LV.MDB5<- as.data.frame(LV.MDB5[rowMeans(LV.MDB5)!= 0,])
LV.MDB5$mean<-round(rowMeans(LV.MDB5))
LV.MDB5$cumulate<-LV.MDB5$mean
for (i in 2:length(LV.MDB5$cumulate))
  LV.MDB5$cumulate[i] <- LV.MDB5$mean[i]+LV.MDB5$cumulate[i-1]
LV.MDB5$method<-"MDB5soil"
LV.MDB5$index<-row.names(LV.MDB5)

LV.Zobel<-as.data.frame(LV[,names(LV)=="Zobelsoil"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"Zobelsoil"
LV.Zobel$index<-row.names(LV.Zobel)
#
all_data<-rbind( 
  LV.SUCC[,3:5],
  LV.DarkDiv[,3:5],
  LV.GSMc[,3:5],
  LV.deep[,3:5],
  LV.deep_SUCC[,3:5],
  LV.GSMc40A[,3:5],
  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc_62DNA"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc_62DNA"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40ADNA"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40ADNA"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40BDNA"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40BDNA"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="MaestreDNA"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"MaestreDNA"
LV.Maestre$index<-row.names(LV.Maestre)


LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDivDNA"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDivDNA"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCCDNA"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCCDNA"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCASDNA"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCASDNA"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="deepDNA"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deepDNA"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCCDNA"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCCDNA"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.MDB15<-as.data.frame(LV[,names(LV)=="MDB15DNA"])
LV.MDB15<- as.data.frame(LV.MDB15[rowMeans(LV.MDB15)!= 0,])
LV.MDB15$mean<-round(rowMeans(LV.MDB15))
LV.MDB15$cumulate<-LV.MDB15$mean
for (i in 2:length(LV.MDB15$cumulate))
  LV.MDB15$cumulate[i] <- LV.MDB15$mean[i]+LV.MDB15$cumulate[i-1]
LV.MDB15$method<-"MDB15DNA"
LV.MDB15$index<-row.names(LV.MDB15)

LV.MDB5<-as.data.frame(LV[,names(LV)=="MDB5DNA"])
LV.MDB5<- as.data.frame(LV.MDB5[rowMeans(LV.MDB5)!= 0,])
LV.MDB5$mean<-round(rowMeans(LV.MDB5))
LV.MDB5$cumulate<-LV.MDB5$mean
for (i in 2:length(LV.MDB5$cumulate))
  LV.MDB5$cumulate[i] <- LV.MDB5$mean[i]+LV.MDB5$cumulate[i-1]
LV.MDB5$method<-"MDB5DNA"
LV.MDB5$index<-row.names(LV.MDB5)

LV.Zobel<-as.data.frame(LV[,names(LV)=="ZobelDNA"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"ZobelDNA"
LV.Zobel$index<-row.names(LV.Zobel)
#
all_data<-rbind(  LV.Zobel[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)
library(ggplot2)
out<-ggplot(all_data, aes(x = index, y = cumulate, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)


out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##LW
LW <- normal.fre2[grepl("LW",normal.fre2$site2),]
row.names(LW)<-LW$method
LW <-LW[,1:(ncol(LW)-2)]
## cumulative abundance
LW <- as.data.frame(t(LW))
LW<-as.data.frame(apply(LW,2,as.numeric))
for (i in 1:length(colnames(LW))){
  LW[,i] <- LW[,i][order(LW[,i], decreasing = T)]
  
}

method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LW <- LW[rowMeans(LW)!= 0,]
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40A"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40A"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40B"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40B"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDiv"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDiv"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCC"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCC"
LW.SUCC$index<-row.names(LW.SUCC)

LW.deep<-as.data.frame(LW[,names(LW)=="Deep"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])
LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deep"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCC"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCC"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCAS"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCAS"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.Zobel<-as.data.frame(LW[,names(LW)=="Zobel"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"Zobel"
LW.Zobel$index<-row.names(LW.Zobel)

#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

###for the pooled
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc_62soil"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc_62soil"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40Asoil"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40Asoil"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40Bsoil"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40Bsoil"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.Maestre<-as.data.frame(LW[,names(LW)=="Maestresoil"])
LW.Maestre<- as.data.frame(LW.Maestre[rowMeans(LW.Maestre)!= 0,])
LW.Maestre$mean<-round(rowMeans(LW.Maestre))
LW.Maestre$cumulate<-LW.Maestre$mean
for (i in 2:length(LW.Maestre$cumulate))
  LW.Maestre$cumulate[i] <- LW.Maestre$mean[i]+LW.Maestre$cumulate[i-1]
LW.Maestre$method<-"Maestresoil"
LW.Maestre$index<-row.names(LW.Maestre)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDivsoil"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDivsoil"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCCsoil"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCCsoil"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCASsoil"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCASsoil"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.deep<-as.data.frame(LW[,names(LW)=="deepsoil"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])

LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deepsoil"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCCsoil"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCCsoil"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.MDB15<-as.data.frame(LW[,names(LW)=="MDB15soil"])
LW.MDB15<- as.data.frame(LW.MDB15[rowMeans(LW.MDB15)!= 0,])
LW.MDB15$mean<-round(rowMeans(LW.MDB15))
LW.MDB15$cumulate<-LW.MDB15$mean
for (i in 2:length(LW.MDB15$cumulate))
  LW.MDB15$cumulate[i] <- LW.MDB15$mean[i]+LW.MDB15$cumulate[i-1]
LW.MDB15$method<-"MDB15soil"
LW.MDB15$index<-row.names(LW.MDB15)

LW.MDB5<-as.data.frame(LW[,names(LW)=="MDB5soil"])
LW.MDB5<- as.data.frame(LW.MDB5[rowMeans(LW.MDB5)!= 0,])
LW.MDB5$mean<-round(rowMeans(LW.MDB5))
LW.MDB5$cumulate<-LW.MDB5$mean
for (i in 2:length(LW.MDB5$cumulate))
  LW.MDB5$cumulate[i] <- LW.MDB5$mean[i]+LW.MDB5$cumulate[i-1]
LW.MDB5$method<-"MDB5soil"
LW.MDB5$index<-row.names(LW.MDB5)

LW.Zobel<-as.data.frame(LW[,names(LW)=="Zobelsoil"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"Zobelsoil"
LW.Zobel$index<-row.names(LW.Zobel)
#

all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))
methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc_62DNA"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc_62DNA"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40ADNA"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40ADNA"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40BDNA"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40BDNA"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.Maestre<-as.data.frame(LW[,names(LW)=="MaestreDNA"])
LW.Maestre<- as.data.frame(LW.Maestre[rowMeans(LW.Maestre)!= 0,])
LW.Maestre$mean<-round(rowMeans(LW.Maestre))
LW.Maestre$cumulate<-LW.Maestre$mean
for (i in 2:length(LW.Maestre$cumulate))
  LW.Maestre$cumulate[i] <- LW.Maestre$mean[i]+LW.Maestre$cumulate[i-1]
LW.Maestre$method<-"MaestreDNA"
LW.Maestre$index<-row.names(LW.Maestre)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDivDNA"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDivDNA"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCCDNA"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCCDNA"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCASDNA"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCASDNA"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.deep<-as.data.frame(LW[,names(LW)=="deepDNA"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])

LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deepDNA"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCCDNA"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCCDNA"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.MDB15<-as.data.frame(LW[,names(LW)=="MDB15DNA"])
LW.MDB15<- as.data.frame(LW.MDB15[rowMeans(LW.MDB15)!= 0,])
LW.MDB15$mean<-round(rowMeans(LW.MDB15))
LW.MDB15$cumulate<-LW.MDB15$mean
for (i in 2:length(LW.MDB15$cumulate))
  LW.MDB15$cumulate[i] <- LW.MDB15$mean[i]+LW.MDB15$cumulate[i-1]
LW.MDB15$method<-"MDB15DNA"
LW.MDB15$index<-row.names(LW.MDB15)

LW.MDB5<-as.data.frame(LW[,names(LW)=="MDB5DNA"])
LW.MDB5<- as.data.frame(LW.MDB5[rowMeans(LW.MDB5)!= 0,])
LW.MDB5$mean<-round(rowMeans(LW.MDB5))
LW.MDB5$cumulate<-LW.MDB5$mean
for (i in 2:length(LW.MDB5$cumulate))
  LW.MDB5$cumulate[i] <- LW.MDB5$mean[i]+LW.MDB5$cumulate[i-1]
LW.MDB5$method<-"MDB5DNA"
LW.MDB5$index<-row.names(LW.MDB5)

LW.Zobel<-as.data.frame(LW[,names(LW)=="ZobelDNA"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"ZobelDNA"
LW.Zobel$index<-row.names(LW.Zobel)
#

all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
#bacteria
normal<-read.table("bacteria.rarefy.table.csv",header = T,row.names = 1,sep = ",")

normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata2<-metadata
metadata$site<-paste0(metadata$site,"DNA")
metadata2$site<-paste0(metadata2$site,"soil")
metadata$method<-paste0(metadata$method,"DNA")
metadata2$method<-paste0(metadata2$method,"soil")
metadata<-rbind(metadata,metadata2)

metadata.normal2<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
metadata<-metadata[,c(1,5)]
metadata<-rbind(metadata,metadata.normal2)

normal.fre2<-merge(normal.fre,metadata,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$'0000294ab66fb96003e4afd4237026e0'),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]

normal.fre2<-normal.fre2[normal.fre2$method!="Meastre",]

normal.fre2<-aggregate(normal.fre2[,2:(ncol(normal.fre2)-2)],by=list(method=normal.fre2$method,site2=normal.fre2$site2),sum)
###for GSMc40
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-normal[,grepl(GSMc40A,names(normal))]
names(A40)<-paste0(names(A40),"40A")

B40<-normal[,grepl(GSMc40B,names(normal))]
names(B40)<-paste0(names(B40),"40B")

metadata.normal<-rbind(
  data.frame(method="GSMc_40A",site=names(A40)),
  data.frame(method="GSMc_40B",site=names(B40)),
  metadata)

try<-function(x,y){
  merge(x,y,by="OTU",all=T)
}
A40$OTU<-row.names(A40)
B40$OTU<-row.names(B40)

a<-list(A40,B40)

bacteria.LZ<-Reduce(try,a)
row.names(bacteria.LZ)<-bacteria.LZ$OTU
bacteria.LZ<-bacteria.LZ[,-1]
bacteria.LZ<-as.data.frame(t(bacteria.LZ))
bacteria.LZ$site<-row.names(bacteria.LZ)

bacteria.LZ2<-merge(bacteria.LZ,metadata.normal,by="site")
bacteria.LZ2$site2<-substr(bacteria.LZ2$site,1,2)
bacteria.LZ2<-aggregate(bacteria.LZ2[,2:(ncol(bacteria.LZ2)-2)],by=list(method=bacteria.LZ2$method,site2=bacteria.LZ2$site2),sum)

bacteria.LZ2<-as.data.frame(t(bacteria.LZ2))
normal.fre2<-as.data.frame(t(normal.fre2))
normal.fre2$me<-row.names(normal.fre2)
bacteria.LZ2$me<-row.names(bacteria.LZ2)
a<-merge(normal.fre2,bacteria.LZ2,by="me",all=T)
row.names(a)<-a$me
a<-a[,-1]
a<-as.data.frame(t(a))
normal.fre2<-a
###
deep_SUCC<-normal.fre2[grepl("Deep|SUCC",normal.fre2$method),]
deep_SUCC<-deep_SUCC[!grepl("DNA|soil",deep_SUCC$method),]
deep_SUCC<-deep_SUCC[!grepl("deep_SUCC",deep_SUCC$method),]
deep_SUCC[,1:(ncol(deep_SUCC)-2)]<-apply(deep_SUCC[,1:(ncol(deep_SUCC)-2)],2,as.numeric)
deep_SUCC<-aggregate(deep_SUCC[,1:(ncol(deep_SUCC)-2)],by=list(site2=deep_SUCC$site2),sum)
deep_SUCC$method<-"deep_SUCC"
normal.fre2<-rbind(normal.fre2,deep_SUCC)

##LZ
LZ <- normal.fre2[grepl("LZ",normal.fre2$site2),]
row.names(LZ)<-LZ$method
LZ <-LZ[,1:(ncol(LZ)-2)]
## cumulative abundance
LZ <- as.data.frame(t(LZ))
LZ<-as.data.frame(apply(LZ,2,as.numeric))
for (i in 1:length(colnames(LZ))){
  LZ[,i] <- LZ[,i][order(LZ[,i], decreasing = T)]
  
}
method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LZ <- LZ[rowMeans(LZ)!= 0,]
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40A"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40A"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40B"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40B"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="Maestre"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"Maestre"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDiv"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDiv"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCC"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCC"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCAS"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCAS"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="Deep"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])
LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deep"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCC"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCC"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="Zobel"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"Zobel"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
###for the pooled
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc_62soil"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc_62soil"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40Asoil"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40Asoil"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40Bsoil"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40Bsoil"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="Maestresoil"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"Maestresoil"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDivsoil"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDivsoil"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCCsoil"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCCsoil"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCASsoil"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCASsoil"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="deepsoil"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deepsoil"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCCsoil"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCCsoil"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.MDB15<-as.data.frame(LZ[,names(LZ)=="MDB15soil"])
LZ.MDB15<- as.data.frame(LZ.MDB15[rowMeans(LZ.MDB15)!= 0,])
LZ.MDB15$mean<-round(rowMeans(LZ.MDB15))
LZ.MDB15$cumulate<-LZ.MDB15$mean
for (i in 2:length(LZ.MDB15$cumulate))
  LZ.MDB15$cumulate[i] <- LZ.MDB15$mean[i]+LZ.MDB15$cumulate[i-1]
LZ.MDB15$method<-"MDB15soil"
LZ.MDB15$index<-row.names(LZ.MDB15)

LZ.MDB5<-as.data.frame(LZ[,names(LZ)=="MDB5soil"])
LZ.MDB5<- as.data.frame(LZ.MDB5[rowMeans(LZ.MDB5)!= 0,])
LZ.MDB5$mean<-round(rowMeans(LZ.MDB5))
LZ.MDB5$cumulate<-LZ.MDB5$mean
for (i in 2:length(LZ.MDB5$cumulate))
  LZ.MDB5$cumulate[i] <- LZ.MDB5$mean[i]+LZ.MDB5$cumulate[i-1]
LZ.MDB5$method<-"MDB5soil"
LZ.MDB5$index<-row.names(LZ.MDB5)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="Zobelsoil"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"Zobelsoil"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LZ.GSMc<-as.data.frame(LZ[,names(LZ)=="GSMc_62DNA"])
LZ.GSMc<- as.data.frame(LZ.GSMc[rowMeans(LZ.GSMc)!= 0,])
LZ.GSMc$mean<-round(rowMeans(LZ.GSMc))
LZ.GSMc$cumulate<-LZ.GSMc$mean
for (i in 2:length(LZ.GSMc$cumulate))
  LZ.GSMc$cumulate[i] <- LZ.GSMc$mean[i]+LZ.GSMc$cumulate[i-1]
LZ.GSMc$method<-"GSMc_62DNA"
LZ.GSMc$index<-row.names(LZ.GSMc)

LZ.GSMc40A<-as.data.frame(LZ[,names(LZ)=="GSMc_40ADNA"])
LZ.GSMc40A<- as.data.frame(LZ.GSMc40A[rowMeans(LZ.GSMc40A)!= 0,])
LZ.GSMc40A$mean<-round(rowMeans(LZ.GSMc40A))
LZ.GSMc40A$cumulate<-LZ.GSMc40A$mean
for (i in 2:length(LZ.GSMc40A$cumulate))
  LZ.GSMc40A$cumulate[i] <- LZ.GSMc40A$mean[i]+LZ.GSMc40A$cumulate[i-1]
LZ.GSMc40A$method<-"GSMc_40ADNA"
LZ.GSMc40A$index<-row.names(LZ.GSMc40A)

LZ.GSMc40B<-as.data.frame(LZ[,names(LZ)=="GSMc_40BDNA"])
LZ.GSMc40B<- as.data.frame(LZ.GSMc40B[rowMeans(LZ.GSMc40B)!= 0,])
LZ.GSMc40B$mean<-round(rowMeans(LZ.GSMc40B))
LZ.GSMc40B$cumulate<-LZ.GSMc40B$mean
for (i in 2:length(LZ.GSMc40B$cumulate))
  LZ.GSMc40B$cumulate[i] <- LZ.GSMc40B$mean[i]+LZ.GSMc40B$cumulate[i-1]
LZ.GSMc40B$method<-"GSMc_40BDNA"
LZ.GSMc40B$index<-row.names(LZ.GSMc40B)

LZ.Maestre<-as.data.frame(LZ[,names(LZ)=="MaestreDNA"])
LZ.Maestre<- as.data.frame(LZ.Maestre[rowMeans(LZ.Maestre)!= 0,])
LZ.Maestre$mean<-round(rowMeans(LZ.Maestre))
LZ.Maestre$cumulate<-LZ.Maestre$mean
for (i in 2:length(LZ.Maestre$cumulate))
  LZ.Maestre$cumulate[i] <- LZ.Maestre$mean[i]+LZ.Maestre$cumulate[i-1]
LZ.Maestre$method<-"MaestreDNA"
LZ.Maestre$index<-row.names(LZ.Maestre)

LZ.DarkDiv<-as.data.frame(LZ[,names(LZ)=="DarkDivDNA"])
LZ.DarkDiv<- as.data.frame(LZ.DarkDiv[rowMeans(LZ.DarkDiv)!= 0,])
LZ.DarkDiv$mean<-round(rowMeans(LZ.DarkDiv))
LZ.DarkDiv$cumulate<-LZ.DarkDiv$mean
for (i in 2:length(LZ.DarkDiv$cumulate))
  LZ.DarkDiv$cumulate[i] <- LZ.DarkDiv$mean[i]+LZ.DarkDiv$cumulate[i-1]
LZ.DarkDiv$method<-"DarkDivDNA"
LZ.DarkDiv$index<-row.names(LZ.DarkDiv)

LZ.SUCC<-as.data.frame(LZ[,names(LZ)=="SUCCDNA"])
LZ.SUCC<- as.data.frame(LZ.SUCC[rowMeans(LZ.SUCC)!= 0,])
LZ.SUCC$mean<-round(rowMeans(LZ.SUCC))
LZ.SUCC$cumulate<-LZ.SUCC$mean
for (i in 2:length(LZ.SUCC$cumulate))
  LZ.SUCC$cumulate[i] <- LZ.SUCC$mean[i]+LZ.SUCC$cumulate[i-1]
LZ.SUCC$method<-"SUCCDNA"
LZ.SUCC$index<-row.names(LZ.SUCC)

LZ.LUCAS<-as.data.frame(LZ[,names(LZ)=="LUCASDNA"])
LZ.LUCAS<- as.data.frame(LZ.LUCAS[rowMeans(LZ.LUCAS)!= 0,])
LZ.LUCAS$mean<-round(rowMeans(LZ.LUCAS))
LZ.LUCAS$cumulate<-LZ.LUCAS$mean
for (i in 2:length(LZ.LUCAS$cumulate))
  LZ.LUCAS$cumulate[i] <- LZ.LUCAS$mean[i]+LZ.LUCAS$cumulate[i-1]
LZ.LUCAS$method<-"LUCASDNA"
LZ.LUCAS$index<-row.names(LZ.LUCAS)

LZ.deep<-as.data.frame(LZ[,names(LZ)=="deepDNA"])
LZ.deep<- as.data.frame(LZ.deep[rowMeans(LZ.deep)!= 0,])

LZ.deep$mean<-round(rowMeans(LZ.deep))
LZ.deep$cumulate<-LZ.deep$mean
for (i in 2:length(LZ.deep$cumulate))
  LZ.deep$cumulate[i] <- LZ.deep$mean[i]+LZ.deep$cumulate[i-1]
LZ.deep$method<-"deepDNA"
LZ.deep$index<-row.names(LZ.deep)

LZ.deep_SUCC<-as.data.frame(LZ[,names(LZ)=="deep_SUCCDNA"])
LZ.deep_SUCC<- as.data.frame(LZ.deep_SUCC[rowMeans(LZ.deep_SUCC)!= 0,])
LZ.deep_SUCC$mean<-round(rowMeans(LZ.deep_SUCC))
LZ.deep_SUCC$cumulate<-LZ.deep_SUCC$mean
for (i in 2:length(LZ.deep_SUCC$cumulate))
  LZ.deep_SUCC$cumulate[i] <- LZ.deep_SUCC$mean[i]+LZ.deep_SUCC$cumulate[i-1]
LZ.deep_SUCC$method<-"deep_SUCCDNA"
LZ.deep_SUCC$index<-row.names(LZ.deep_SUCC)

LZ.MDB15<-as.data.frame(LZ[,names(LZ)=="MDB15DNA"])
LZ.MDB15<- as.data.frame(LZ.MDB15[rowMeans(LZ.MDB15)!= 0,])
LZ.MDB15$mean<-round(rowMeans(LZ.MDB15))
LZ.MDB15$cumulate<-LZ.MDB15$mean
for (i in 2:length(LZ.MDB15$cumulate))
  LZ.MDB15$cumulate[i] <- LZ.MDB15$mean[i]+LZ.MDB15$cumulate[i-1]
LZ.MDB15$method<-"MDB15DNA"
LZ.MDB15$index<-row.names(LZ.MDB15)

LZ.MDB5<-as.data.frame(LZ[,names(LZ)=="MDB5DNA"])
LZ.MDB5<- as.data.frame(LZ.MDB5[rowMeans(LZ.MDB5)!= 0,])
LZ.MDB5$mean<-round(rowMeans(LZ.MDB5))
LZ.MDB5$cumulate<-LZ.MDB5$mean
for (i in 2:length(LZ.MDB5$cumulate))
  LZ.MDB5$cumulate[i] <- LZ.MDB5$mean[i]+LZ.MDB5$cumulate[i-1]
LZ.MDB5$method<-"MDB5DNA"
LZ.MDB5$index<-row.names(LZ.MDB5)

LZ.Zobel<-as.data.frame(LZ[,names(LZ)=="ZobelDNA"])
LZ.Zobel<- as.data.frame(LZ.Zobel[rowMeans(LZ.Zobel)!= 0,])
LZ.Zobel$mean<-round(rowMeans(LZ.Zobel))
LZ.Zobel$cumulate<-LZ.Zobel$mean
for (i in 2:length(LZ.Zobel$cumulate))
  LZ.Zobel$cumulate[i] <- LZ.Zobel$mean[i]+LZ.Zobel$cumulate[i-1]
LZ.Zobel$method<-"ZobelDNA"
LZ.Zobel$index<-row.names(LZ.Zobel)
#
all_data<-rbind(  LZ.Zobel[,3:5],
                  LZ.LUCAS[,3:5],
                  LZ.SUCC[,3:5],
                  LZ.DarkDiv[,3:5],
                  LZ.GSMc[,3:5],
                  LZ.deep[,3:5],
                  LZ.deep_SUCC[,3:5],
                  LZ.GSMc40A[,3:5],
                  LZ.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )
out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##LV
LV <- normal.fre2[grepl("LV",normal.fre2$site2),]
row.names(LV)<-LV$method
LV <-LV[,1:(ncol(LV)-2)]
## cumulative abundance
LV <- as.data.frame(t(LV))
LV<-as.data.frame(apply(LV,2,as.numeric))
for (i in 1:length(colnames(LV))){
  LV[,i] <- LV[,i][order(LV[,i], decreasing = T)]
  
}
method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LV <- LV[rowMeans(LV)!= 0,]
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40A"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40A"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40B"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40B"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="Maestre"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"Maestre"
LV.Maestre$index<-row.names(LV.Maestre)

LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDiv"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDiv"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCC"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCC"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCAS"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCAS"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="Deep"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deep"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCC"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCC"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.Zobel<-as.data.frame(LV[,names(LV)=="Zobel"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"Zobel"
LV.Zobel$index<-row.names(LV.Zobel)
#
all_data<-rbind(  LV.Zobel[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.LUCAS[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)
###get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
###for the pooled
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc_62soil"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc_62soil"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40Asoil"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40Asoil"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40Bsoil"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40Bsoil"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="Maestresoil"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"Maestresoil"
LV.Maestre$index<-row.names(LV.Maestre)


LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDivsoil"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDivsoil"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCCsoil"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCCsoil"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCASsoil"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCASsoil"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="deepsoil"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deepsoil"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCCsoil"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCCsoil"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.MDB15<-as.data.frame(LV[,names(LV)=="MDB15soil"])
LV.MDB15<- as.data.frame(LV.MDB15[rowMeans(LV.MDB15)!= 0,])
LV.MDB15$mean<-round(rowMeans(LV.MDB15))
LV.MDB15$cumulate<-LV.MDB15$mean
for (i in 2:length(LV.MDB15$cumulate))
  LV.MDB15$cumulate[i] <- LV.MDB15$mean[i]+LV.MDB15$cumulate[i-1]
LV.MDB15$method<-"MDB15soil"
LV.MDB15$index<-row.names(LV.MDB15)

LV.MDB5<-as.data.frame(LV[,names(LV)=="MDB5soil"])
LV.MDB5<- as.data.frame(LV.MDB5[rowMeans(LV.MDB5)!= 0,])
LV.MDB5$mean<-round(rowMeans(LV.MDB5))
LV.MDB5$cumulate<-LV.MDB5$mean
for (i in 2:length(LV.MDB5$cumulate))
  LV.MDB5$cumulate[i] <- LV.MDB5$mean[i]+LV.MDB5$cumulate[i-1]
LV.MDB5$method<-"MDB5soil"
LV.MDB5$index<-row.names(LV.MDB5)

LV.Zobel<-as.data.frame(LV[,names(LV)=="Zobelsoil"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"Zobelsoil"
LV.Zobel$index<-row.names(LV.Zobel)
#

all_data<-rbind(  LV.Zobel[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.LUCAS[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LV.GSMc<-as.data.frame(LV[,names(LV)=="GSMc_62DNA"])
LV.GSMc<- as.data.frame(LV.GSMc[rowMeans(LV.GSMc)!= 0,])
LV.GSMc$mean<-round(rowMeans(LV.GSMc))
LV.GSMc$cumulate<-LV.GSMc$mean
for (i in 2:length(LV.GSMc$cumulate))
  LV.GSMc$cumulate[i] <- LV.GSMc$mean[i]+LV.GSMc$cumulate[i-1]
LV.GSMc$method<-"GSMc_62DNA"
LV.GSMc$index<-row.names(LV.GSMc)

LV.GSMc40A<-as.data.frame(LV[,names(LV)=="GSMc_40ADNA"])
LV.GSMc40A<- as.data.frame(LV.GSMc40A[rowMeans(LV.GSMc40A)!= 0,])
LV.GSMc40A$mean<-round(rowMeans(LV.GSMc40A))
LV.GSMc40A$cumulate<-LV.GSMc40A$mean
for (i in 2:length(LV.GSMc40A$cumulate))
  LV.GSMc40A$cumulate[i] <- LV.GSMc40A$mean[i]+LV.GSMc40A$cumulate[i-1]
LV.GSMc40A$method<-"GSMc_40ADNA"
LV.GSMc40A$index<-row.names(LV.GSMc40A)

LV.GSMc40B<-as.data.frame(LV[,names(LV)=="GSMc_40BDNA"])
LV.GSMc40B<- as.data.frame(LV.GSMc40B[rowMeans(LV.GSMc40B)!= 0,])
LV.GSMc40B$mean<-round(rowMeans(LV.GSMc40B))
LV.GSMc40B$cumulate<-LV.GSMc40B$mean
for (i in 2:length(LV.GSMc40B$cumulate))
  LV.GSMc40B$cumulate[i] <- LV.GSMc40B$mean[i]+LV.GSMc40B$cumulate[i-1]
LV.GSMc40B$method<-"GSMc_40BDNA"
LV.GSMc40B$index<-row.names(LV.GSMc40B)

LV.Maestre<-as.data.frame(LV[,names(LV)=="MaestreDNA"])
LV.Maestre<- as.data.frame(LV.Maestre[rowMeans(LV.Maestre)!= 0,])
LV.Maestre$mean<-round(rowMeans(LV.Maestre))
LV.Maestre$cumulate<-LV.Maestre$mean
for (i in 2:length(LV.Maestre$cumulate))
  LV.Maestre$cumulate[i] <- LV.Maestre$mean[i]+LV.Maestre$cumulate[i-1]
LV.Maestre$method<-"MaestreDNA"
LV.Maestre$index<-row.names(LV.Maestre)

LV.DarkDiv<-as.data.frame(LV[,names(LV)=="DarkDivDNA"])
LV.DarkDiv<- as.data.frame(LV.DarkDiv[rowMeans(LV.DarkDiv)!= 0,])
LV.DarkDiv$mean<-round(rowMeans(LV.DarkDiv))
LV.DarkDiv$cumulate<-LV.DarkDiv$mean
for (i in 2:length(LV.DarkDiv$cumulate))
  LV.DarkDiv$cumulate[i] <- LV.DarkDiv$mean[i]+LV.DarkDiv$cumulate[i-1]
LV.DarkDiv$method<-"DarkDivDNA"
LV.DarkDiv$index<-row.names(LV.DarkDiv)

LV.SUCC<-as.data.frame(LV[,names(LV)=="SUCCDNA"])
LV.SUCC<- as.data.frame(LV.SUCC[rowMeans(LV.SUCC)!= 0,])
LV.SUCC$mean<-round(rowMeans(LV.SUCC))
LV.SUCC$cumulate<-LV.SUCC$mean
for (i in 2:length(LV.SUCC$cumulate))
  LV.SUCC$cumulate[i] <- LV.SUCC$mean[i]+LV.SUCC$cumulate[i-1]
LV.SUCC$method<-"SUCCDNA"
LV.SUCC$index<-row.names(LV.SUCC)

LV.LUCAS<-as.data.frame(LV[,names(LV)=="LUCASDNA"])
LV.LUCAS<- as.data.frame(LV.LUCAS[rowMeans(LV.LUCAS)!= 0,])
LV.LUCAS$mean<-round(rowMeans(LV.LUCAS))
LV.LUCAS$cumulate<-LV.LUCAS$mean
for (i in 2:length(LV.LUCAS$cumulate))
  LV.LUCAS$cumulate[i] <- LV.LUCAS$mean[i]+LV.LUCAS$cumulate[i-1]
LV.LUCAS$method<-"LUCASDNA"
LV.LUCAS$index<-row.names(LV.LUCAS)

LV.deep<-as.data.frame(LV[,names(LV)=="deepDNA"])
LV.deep<- as.data.frame(LV.deep[rowMeans(LV.deep)!= 0,])

LV.deep$mean<-round(rowMeans(LV.deep))
LV.deep$cumulate<-LV.deep$mean
for (i in 2:length(LV.deep$cumulate))
  LV.deep$cumulate[i] <- LV.deep$mean[i]+LV.deep$cumulate[i-1]
LV.deep$method<-"deepDNA"
LV.deep$index<-row.names(LV.deep)

LV.deep_SUCC<-as.data.frame(LV[,names(LV)=="deep_SUCCDNA"])
LV.deep_SUCC<- as.data.frame(LV.deep_SUCC[rowMeans(LV.deep_SUCC)!= 0,])
LV.deep_SUCC$mean<-round(rowMeans(LV.deep_SUCC))
LV.deep_SUCC$cumulate<-LV.deep_SUCC$mean
for (i in 2:length(LV.deep_SUCC$cumulate))
  LV.deep_SUCC$cumulate[i] <- LV.deep_SUCC$mean[i]+LV.deep_SUCC$cumulate[i-1]
LV.deep_SUCC$method<-"deep_SUCCDNA"
LV.deep_SUCC$index<-row.names(LV.deep_SUCC)

LV.MDB15<-as.data.frame(LV[,names(LV)=="MDB15DNA"])
LV.MDB15<- as.data.frame(LV.MDB15[rowMeans(LV.MDB15)!= 0,])
LV.MDB15$mean<-round(rowMeans(LV.MDB15))
LV.MDB15$cumulate<-LV.MDB15$mean
for (i in 2:length(LV.MDB15$cumulate))
  LV.MDB15$cumulate[i] <- LV.MDB15$mean[i]+LV.MDB15$cumulate[i-1]
LV.MDB15$method<-"MDB15DNA"
LV.MDB15$index<-row.names(LV.MDB15)

LV.MDB5<-as.data.frame(LV[,names(LV)=="MDB5DNA"])
LV.MDB5<- as.data.frame(LV.MDB5[rowMeans(LV.MDB5)!= 0,])
LV.MDB5$mean<-round(rowMeans(LV.MDB5))
LV.MDB5$cumulate<-LV.MDB5$mean
for (i in 2:length(LV.MDB5$cumulate))
  LV.MDB5$cumulate[i] <- LV.MDB5$mean[i]+LV.MDB5$cumulate[i-1]
LV.MDB5$method<-"MDB5DNA"
LV.MDB5$index<-row.names(LV.MDB5)

LV.Zobel<-as.data.frame(LV[,names(LV)=="ZobelDNA"])
LV.Zobel<- as.data.frame(LV.Zobel[rowMeans(LV.Zobel)!= 0,])
LV.Zobel$mean<-round(rowMeans(LV.Zobel))
LV.Zobel$cumulate<-LV.Zobel$mean
for (i in 2:length(LV.Zobel$cumulate))
  LV.Zobel$cumulate[i] <- LV.Zobel$mean[i]+LV.Zobel$cumulate[i-1]
LV.Zobel$method<-"ZobelDNA"
LV.Zobel$index<-row.names(LV.Zobel)
#
all_data<-rbind(  LV.Zobel[,3:5],
                  LV.SUCC[,3:5],
                  LV.DarkDiv[,3:5],
                  LV.LUCAS[,3:5],
                  LV.GSMc[,3:5],
                  LV.deep[,3:5],
                  LV.deep_SUCC[,3:5],
                  LV.GSMc40A[,3:5],
                  LV.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)
cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)

####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##LW
LW <- normal.fre2[grepl("LW",normal.fre2$site2),]
row.names(LW)<-LW$method
LW <-LW[,1:(ncol(LW)-2)]
## cumulative abundance
LW <- as.data.frame(t(LW))
LW<-as.data.frame(apply(LW,2,as.numeric))
for (i in 1:length(colnames(LW))){
  LW[,i] <- LW[,i][order(LW[,i], decreasing = T)]
  
}

method<-unique(metadata.normal$method)
method<-c(method,"deep_SUCC")
## cumulative abundance
LW <- LW[rowMeans(LW)!= 0,]
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40A"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40A"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40B"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40B"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDiv"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDiv"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCC"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCC"
LW.SUCC$index<-row.names(LW.SUCC)

LW.deep<-as.data.frame(LW[,names(LW)=="Deep"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])
LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deep"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCC"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCC"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCAS"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCAS"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.Zobel<-as.data.frame(LW[,names(LW)=="Zobel"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"Zobel"
LW.Zobel$index<-row.names(LW.Zobel)
#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))
methods<-c("DarkDiv","GSMc","GSMc_40A","GSMc_40B","LUCAS","SUCC","Zobel","deep","deep_SUCC")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)

###for the pooled
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc_62soil"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc_62soil"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40Asoil"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40Asoil"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40Bsoil"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40Bsoil"
LW.GSMc40B$index<-row.names(LW.GSMc40B)



LW.Maestre<-as.data.frame(LW[,names(LW)=="Maestresoil"])
LW.Maestre<- as.data.frame(LW.Maestre[rowMeans(LW.Maestre)!= 0,])
LW.Maestre$mean<-round(rowMeans(LW.Maestre))
LW.Maestre$cumulate<-LW.Maestre$mean
for (i in 2:length(LW.Maestre$cumulate))
  LW.Maestre$cumulate[i] <- LW.Maestre$mean[i]+LW.Maestre$cumulate[i-1]
LW.Maestre$method<-"Maestresoil"
LW.Maestre$index<-row.names(LW.Maestre)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDivsoil"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDivsoil"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCCsoil"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCCsoil"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCASsoil"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCASsoil"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.deep<-as.data.frame(LW[,names(LW)=="deepsoil"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])

LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deepsoil"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCCsoil"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCCsoil"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.MDB15<-as.data.frame(LW[,names(LW)=="MDB15soil"])
LW.MDB15<- as.data.frame(LW.MDB15[rowMeans(LW.MDB15)!= 0,])
LW.MDB15$mean<-round(rowMeans(LW.MDB15))
LW.MDB15$cumulate<-LW.MDB15$mean
for (i in 2:length(LW.MDB15$cumulate))
  LW.MDB15$cumulate[i] <- LW.MDB15$mean[i]+LW.MDB15$cumulate[i-1]
LW.MDB15$method<-"MDB15soil"
LW.MDB15$index<-row.names(LW.MDB15)

LW.MDB5<-as.data.frame(LW[,names(LW)=="MDB5soil"])
LW.MDB5<- as.data.frame(LW.MDB5[rowMeans(LW.MDB5)!= 0,])
LW.MDB5$mean<-round(rowMeans(LW.MDB5))
LW.MDB5$cumulate<-LW.MDB5$mean
for (i in 2:length(LW.MDB5$cumulate))
  LW.MDB5$cumulate[i] <- LW.MDB5$mean[i]+LW.MDB5$cumulate[i-1]
LW.MDB5$method<-"MDB5soil"
LW.MDB5$index<-row.names(LW.MDB5)

LW.Zobel<-as.data.frame(LW[,names(LW)=="Zobelsoil"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"Zobelsoil"
LW.Zobel$index<-row.names(LW.Zobel)
#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))

methods<-c("DarkDivsoil","GSMc_62soil","GSMc_40Asoil","GSMc_40Bsoil","LUCASsoil","SUCCsoil","Zobelsoil","deepsoil","deep_SUCCsoil")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
##DNA
LW.GSMc<-as.data.frame(LW[,names(LW)=="GSMc_62DNA"])
LW.GSMc<- as.data.frame(LW.GSMc[rowMeans(LW.GSMc)!= 0,])
LW.GSMc$mean<-round(rowMeans(LW.GSMc))
LW.GSMc$cumulate<-LW.GSMc$mean
for (i in 2:length(LW.GSMc$cumulate))
  LW.GSMc$cumulate[i] <- LW.GSMc$mean[i]+LW.GSMc$cumulate[i-1]
LW.GSMc$method<-"GSMc_62DNA"
LW.GSMc$index<-row.names(LW.GSMc)

LW.GSMc40A<-as.data.frame(LW[,names(LW)=="GSMc_40ADNA"])
LW.GSMc40A<- as.data.frame(LW.GSMc40A[rowMeans(LW.GSMc40A)!= 0,])
LW.GSMc40A$mean<-round(rowMeans(LW.GSMc40A))
LW.GSMc40A$cumulate<-LW.GSMc40A$mean
for (i in 2:length(LW.GSMc40A$cumulate))
  LW.GSMc40A$cumulate[i] <- LW.GSMc40A$mean[i]+LW.GSMc40A$cumulate[i-1]
LW.GSMc40A$method<-"GSMc_40ADNA"
LW.GSMc40A$index<-row.names(LW.GSMc40A)

LW.GSMc40B<-as.data.frame(LW[,names(LW)=="GSMc_40BDNA"])
LW.GSMc40B<- as.data.frame(LW.GSMc40B[rowMeans(LW.GSMc40B)!= 0,])
LW.GSMc40B$mean<-round(rowMeans(LW.GSMc40B))
LW.GSMc40B$cumulate<-LW.GSMc40B$mean
for (i in 2:length(LW.GSMc40B$cumulate))
  LW.GSMc40B$cumulate[i] <- LW.GSMc40B$mean[i]+LW.GSMc40B$cumulate[i-1]
LW.GSMc40B$method<-"GSMc_40BDNA"
LW.GSMc40B$index<-row.names(LW.GSMc40B)

LW.Maestre<-as.data.frame(LW[,names(LW)=="MaestreDNA"])
LW.Maestre<- as.data.frame(LW.Maestre[rowMeans(LW.Maestre)!= 0,])
LW.Maestre$mean<-round(rowMeans(LW.Maestre))
LW.Maestre$cumulate<-LW.Maestre$mean
for (i in 2:length(LW.Maestre$cumulate))
  LW.Maestre$cumulate[i] <- LW.Maestre$mean[i]+LW.Maestre$cumulate[i-1]
LW.Maestre$method<-"MaestreDNA"
LW.Maestre$index<-row.names(LW.Maestre)

LW.DarkDiv<-as.data.frame(LW[,names(LW)=="DarkDivDNA"])
LW.DarkDiv<- as.data.frame(LW.DarkDiv[rowMeans(LW.DarkDiv)!= 0,])
LW.DarkDiv$mean<-round(rowMeans(LW.DarkDiv))
LW.DarkDiv$cumulate<-LW.DarkDiv$mean
for (i in 2:length(LW.DarkDiv$cumulate))
  LW.DarkDiv$cumulate[i] <- LW.DarkDiv$mean[i]+LW.DarkDiv$cumulate[i-1]
LW.DarkDiv$method<-"DarkDivDNA"
LW.DarkDiv$index<-row.names(LW.DarkDiv)

LW.SUCC<-as.data.frame(LW[,names(LW)=="SUCCDNA"])
LW.SUCC<- as.data.frame(LW.SUCC[rowMeans(LW.SUCC)!= 0,])
LW.SUCC$mean<-round(rowMeans(LW.SUCC))
LW.SUCC$cumulate<-LW.SUCC$mean
for (i in 2:length(LW.SUCC$cumulate))
  LW.SUCC$cumulate[i] <- LW.SUCC$mean[i]+LW.SUCC$cumulate[i-1]
LW.SUCC$method<-"SUCCDNA"
LW.SUCC$index<-row.names(LW.SUCC)

LW.LUCAS<-as.data.frame(LW[,names(LW)=="LUCASDNA"])
LW.LUCAS<- as.data.frame(LW.LUCAS[rowMeans(LW.LUCAS)!= 0,])
LW.LUCAS$mean<-round(rowMeans(LW.LUCAS))
LW.LUCAS$cumulate<-LW.LUCAS$mean
for (i in 2:length(LW.LUCAS$cumulate))
  LW.LUCAS$cumulate[i] <- LW.LUCAS$mean[i]+LW.LUCAS$cumulate[i-1]
LW.LUCAS$method<-"LUCASDNA"
LW.LUCAS$index<-row.names(LW.LUCAS)

LW.deep<-as.data.frame(LW[,names(LW)=="deepDNA"])
LW.deep<- as.data.frame(LW.deep[rowMeans(LW.deep)!= 0,])

LW.deep$mean<-round(rowMeans(LW.deep))
LW.deep$cumulate<-LW.deep$mean
for (i in 2:length(LW.deep$cumulate))
  LW.deep$cumulate[i] <- LW.deep$mean[i]+LW.deep$cumulate[i-1]
LW.deep$method<-"deepDNA"
LW.deep$index<-row.names(LW.deep)

LW.deep_SUCC<-as.data.frame(LW[,names(LW)=="deep_SUCCDNA"])
LW.deep_SUCC<- as.data.frame(LW.deep_SUCC[rowMeans(LW.deep_SUCC)!= 0,])
LW.deep_SUCC$mean<-round(rowMeans(LW.deep_SUCC))
LW.deep_SUCC$cumulate<-LW.deep_SUCC$mean
for (i in 2:length(LW.deep_SUCC$cumulate))
  LW.deep_SUCC$cumulate[i] <- LW.deep_SUCC$mean[i]+LW.deep_SUCC$cumulate[i-1]
LW.deep_SUCC$method<-"deep_SUCCDNA"
LW.deep_SUCC$index<-row.names(LW.deep_SUCC)

LW.MDB15<-as.data.frame(LW[,names(LW)=="MDB15DNA"])
LW.MDB15<- as.data.frame(LW.MDB15[rowMeans(LW.MDB15)!= 0,])
LW.MDB15$mean<-round(rowMeans(LW.MDB15))
LW.MDB15$cumulate<-LW.MDB15$mean
for (i in 2:length(LW.MDB15$cumulate))
  LW.MDB15$cumulate[i] <- LW.MDB15$mean[i]+LW.MDB15$cumulate[i-1]
LW.MDB15$method<-"MDB15DNA"
LW.MDB15$index<-row.names(LW.MDB15)

LW.MDB5<-as.data.frame(LW[,names(LW)=="MDB5DNA"])
LW.MDB5<- as.data.frame(LW.MDB5[rowMeans(LW.MDB5)!= 0,])
LW.MDB5$mean<-round(rowMeans(LW.MDB5))
LW.MDB5$cumulate<-LW.MDB5$mean
for (i in 2:length(LW.MDB5$cumulate))
  LW.MDB5$cumulate[i] <- LW.MDB5$mean[i]+LW.MDB5$cumulate[i-1]
LW.MDB5$method<-"MDB5DNA"
LW.MDB5$index<-row.names(LW.MDB5)

LW.Zobel<-as.data.frame(LW[,names(LW)=="ZobelDNA"])
LW.Zobel<- as.data.frame(LW.Zobel[rowMeans(LW.Zobel)!= 0,])
LW.Zobel$mean<-round(rowMeans(LW.Zobel))
LW.Zobel$cumulate<-LW.Zobel$mean
for (i in 2:length(LW.Zobel$cumulate))
  LW.Zobel$cumulate[i] <- LW.Zobel$mean[i]+LW.Zobel$cumulate[i-1]
LW.Zobel$method<-"ZobelDNA"
LW.Zobel$index<-row.names(LW.Zobel)
#
all_data<-rbind(  LW.Zobel[,3:5],
                  LW.LUCAS[,3:5],
                  LW.SUCC[,3:5],
                  LW.DarkDiv[,3:5],
                  LW.GSMc[,3:5],
                  LW.deep[,3:5],
                  LW.deep_SUCC[,3:5],
                  LW.GSMc40A[,3:5],
                  LW.GSMc40B[,3:5]
)
all_data$index<-as.numeric(all_data$index)

cbPalette <- rev(c("#9E9E88","#E69F99","#009E73","#1188D3","#E69F00",
                   "#9Fdd99","#A66EFF" ,"#991166","#E1E199","#000099","#FF88D3"))


methods<-c("DarkDivDNA","GSMc_62DNA","GSMc_40ADNA","GSMc_40BDNA","LUCASDNA","SUCCDNA","ZobelDNA","deepDNA","deep_SUCCDNA")
cbPalette <- setNames(cbPalette, methods)
####get the proportion
table <- all_data%>%
  group_by(method) %>%
  summarise(max=max(cumulate))

a<-merge(table,all_data,by="method",all=T)
a$ratio<-a$cumulate/a$max

names(a)[1]<-"method"
out<-ggplot(a, aes(x = index, y = ratio, color = method)) +
  geom_point(size = 1)+
  scale_color_discrete(name = "method") +
  labs(x="Ranking OTU",y="Accumulative abundance ratio")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )
out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)
###rank of the site-richness for each designs
##animal
animal.normal.m<-read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
library(dplyr)
mean.animal <- animal.normal.m %>%
  group_by(method) %>%
  summarise(mean=mean(mean))
all<-merge(animal.normal.m,mean.animal,by="method")

all$precentage<-(all$mean.x-all$mean.y)/all$mean.y
animal<-all
animal$community<-"animal"
##fungi
fungi.normal.m<-read.csv("fungi.totalrichness.csv",header = TRUE,row.names = 1,sep = ",")
mean.fungi <- fungi.normal.m %>%
  group_by(method) %>%
  summarise(mean=mean(mean))

all<-merge(fungi.normal.m,mean.fungi,by="method")

all$precentage<-(all$mean.x-all$mean.y)/all$mean.y
fungi<-all
fungi$community<-"fungi"
##bacteria
bacteria.normal.m<-read.csv("bacteria.totalrichness.csv",header = TRUE,row.names = 1,sep = ",")
mean.bacteria <- bacteria.normal.m %>%
  group_by(method) %>%
  summarise(mean=mean(mean))
all<-merge(bacteria.normal.m,mean.bacteria,by="method")

all$precentage<-(all$mean.x-all$mean.y)/all$mean.y
bacteria<-all
bacteria$community<-"bacteria"
all<-rbind(animal,bacteria,fungi)
all$method<-gsub("GSMc40A","N40D0-5A",all$method)
all$method<-gsub("GSMc40B","N40D0-5B",all$method)
all$method<-gsub("GSMc","N62D0-5",all$method)
all$method<-gsub("deep_SUCC","N9D0-40",all$method)
all$method<-gsub("deep","N9D30-40",all$method)
all$method<-gsub("SUCC","N9D0-10",all$method)
all$method<-gsub("Zobel","N9D0-1",all$method)
all$method<-gsub("LUCAS","N5D0-20",all$method)
all$method<-gsub("DarkDiv","N8D0-10",all$method)
###
all2 <- all %>%
  group_by(method,community) %>%
  mutate(rank = rank(precentage, ties.method = "first")) %>%
  arrange(method,community, rank) 
all2$rank<-as.numeric(all2$rank)
bacteria<-all2[all2$community %in% "bacteria",]
fungi<-all2[all2$community %in% "fungi",]
animal<-all2[all2$community %in% "animal",]
#
library(ggplot2)
  ggplot(bacteria[,c(1,2,9)], aes( site,method)) + 
  geom_tile(aes(fill = rank))+
  scale_fill_gradient(low = "#EEFFFF", high = "#00EEFF")+
  theme_light()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

  ggplot(animal[,c(1,2,9)], aes( site,method)) + 
    geom_tile(aes(fill = rank))+
    scale_fill_gradient(low = "#EEFFFF", high = "#00EEFF")+
    theme_light()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggplot(fungi[,c(1,2,9)], aes( site,method)) + 
    geom_tile(aes(fill = rank))+
    scale_fill_gradient(low = "#EEFFFF", high = "#00EEFF")+
    theme_light()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

####pooled richness and shannon diversity index comparison
##richness
##bacteria
bacteria.normal.m<-read.csv("bacteria.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m$site<-row.names(bacteria.normal.m)
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
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
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
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
  labs(x="method",y="Richness")+
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
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
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
  labs(x="method",y="Richness")+
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
set.seed(888)
mod01<-lmer(log(estimate.x)~  type+(1|method.x) + (1|site2),data = bacteria2)

b<-avg_comparisons(mod01, variables = list(type = "pairwise")) 
plot_predictions(mod01,condition = "type", transform = exp )+
  labs(x="method",y="Richness")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

##fungi
fungi.normal.m<-read.csv("fungi.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m$site<-row.names(fungi.normal.m)
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
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
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
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
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(limits=c(0,1500))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###soil
soil<-fungi.normal.m[grepl("soil",fungi.normal.m$site),]
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
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
  labs(x="method",y="Richness")+
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
set.seed(888)
mod01<-lmer(log(estimate.x)~  type+(1|method.x) + (1|site2),data = fungi2)

b<-avg_comparisons(mod01, variables = list(type = "pairwise")) 
plot_predictions(mod01,condition = "type", transform = exp )+
  labs(x="method",y="Richness")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
##animal
animal.normal.m<-read.csv("animal.pooled.richness.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m$site<-row.names(animal.normal.m)
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
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

names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
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
  labs(x="method",y="Richness")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","Maestre","DarkDiv","MDB15","MDB5","LUCAS","SUCC","deep_SUCC")))+ 
  scale_y_continuous(limits = c(0,1500))+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###soil
soil<-animal.normal.m[grepl("soil",animal.normal.m$site),]
soil<-soil[soil$method!="deep",]
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
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
  labs(x="method",y="Richness")+
  theme_light() +
  scale_y_continuous(limits = c(0,1500))+
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
set.seed(888)
mod01<-lmer(log(estimate.x)~  type+(1|method.x) + (1|site2),data = animal2)

b<-avg_comparisons(mod01, variables = list(type = "pairwise")) 
plot_predictions(mod01,condition = "type", transform = exp )+
  labs(x="method",y="Richness")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
###shannon diversity index
##bacteria
bacteria.normal.m<-read.csv("bacteria.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.normal.m$site<-row.names(bacteria.normal.m)
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
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

names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
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

names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
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
fungi.normal.m<-read.csv("fungi.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.normal.m$site<-row.names(fungi.normal.m)
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
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
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
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
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
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
##animal
animal.normal.m<-read.csv("animal.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.normal.m$site<-row.names(animal.normal.m)
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
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
names(DNA)[names(DNA)=="method"]<-"method.x"
names(DNA)[names(DNA)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = DNA)
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
names(soil)[names(soil)=="method"]<-"method.x"
names(soil)[names(soil)=="estimate"]<-"estimate.x"
mod01<-lm(estimate.x~  method.x + site2,data = soil)
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

###pooling effect calculation at the summary sequencing depths
##effective number
##animal
#DNA
all<-read.csv("proportion.shannon.all.csv",header = TRUE,row.names = 1,sep = ",")
all<-all[grepl("100%",all$proportion),]
all<-all[,c(2:5,9)]
bacteria<-all[grepl("bacteria",all$community),]
fungi<-all[grepl("fungi",all$community),]
animal<-all[grepl("animal",all$community),]
DNA<-animal[animal$pooling=="DNA",]
DNA<-DNA[DNA$method!="deep",]
DNA<-DNA[DNA$method!="deep_SUCC",]
mod01<-lmer(ratio ~ method + (1|site),data = DNA)
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
mod01<-lmer(ratio ~ method + (1|site),data = soil)
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
  filter(!(site == "LV" & method == "LUCAS"))
DNA <- DNA %>%
  filter(!(site == "LV" & method == "deep"))
DNA <- DNA %>%
  filter(!(site == "LV" & method == "deep_SUCC"))
mod01<-lmer(ratio ~ method + (1|site),data = DNA)
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
mod01<-lmer(ratio ~ method + (1|site),data = soil)
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
mod01<-lmer(ratio ~ method + (1|site),data = DNA)
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
soil<-all[all$pooling=="soil",]
mod01<-lmer(ratio ~ method + (1|site),data = soil)
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

##richness
##animal
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
#DNA
all<-read.csv("proportion.all.csv",header = TRUE,row.names = 1,sep = ",")
all<-all[grepl("100%",all$proportion),]
all<-all[,c(2:5,9)]
bacteria<-all[grepl("bacteria",all$community),]
fungi<-all[grepl("fungi",all$community),]
animal<-all[grepl("animal",all$community),]
DNA<-animal[animal$pooling=="DNA",]
DNA<-DNA[DNA$method!="deep",]
DNA<-DNA[DNA$method!="deep_SUCC",]
mod01<-lmer(ratio ~ method + (1|site),data = DNA)
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
mod01<-lmer(ratio ~ method + (1|site),data = soil)
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
  filter(!(site == "LV" & method == "LUCAS"))
DNA <- DNA %>%
  filter(!(site == "LV" & method == "deep"))
DNA <- DNA %>%
  filter(!(site == "LV" & method == "deep_SUCC"))
mod01<-lmer(ratio ~ method + (1|site),data = DNA)
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
mod01<-lmer(ratio ~ method + (1|site),data = soil)
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
mod01<-lmer(ratio ~ method + (1|site),data = DNA)
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
soil<-all[all$pooling=="soil",]
mod01<-lmer(ratio ~ method + (1|site),data = soil)
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
###pooling effect at different sequencing depths
##richness
#bacteria
library(scales)
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
library(ggplot2)
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
library(scales)
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
library(dplyr)
soil <- soil %>%
  filter(!(method == "deep"))
soil <- soil %>%
  filter(!(method == "deep_SUCC"))

mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = soil)
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
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = DNA)
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
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
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
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = soil)
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
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = DNA)
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
bacteria[bacteria$number==0.25,]$try<-"25%"
bacteria[bacteria$number>=0.50&bacteria$number<0.75,]$try<-"=50%~<75%"
bacteria[bacteria$number==0.75,]$try<-"75%"
bacteria[bacteria$number<0.95&bacteria$number>0.75,]$try<-"<75%~95%"
bacteria[bacteria$number==1,]$try<-"100%"
bacteria[bacteria$number>1,]$try<-">100%"
soil<-bacteria[grepl("soil",bacteria$pooling),]
soil<-soil[!is.na(soil$ratio),]

mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = soil)
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
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
  scale_x_discrete(limits=as.character(c( "<25%","25%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
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
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = DNA)
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
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
  scale_x_discrete(limits=as.character(c( "<25%","25%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
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
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
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

###effective number
#Shannon
##bacteria soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.soil.m<-read.csv("bacteria.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.soil.m<-bacteria.soil.m[grepl("soil",row.names(bacteria.soil.m)),]
row.names(bacteria.soil.m)<-gsub("soil","",row.names(bacteria.soil.m))
bacteria.soil.m$site<-row.names(bacteria.soil.m)
bacteria.soil.m<-merge(bacteria.soil.m,soil.metadata,by="site",all=T)
bacteria.soil.m$site2<-substr(bacteria.soil.m$site,1,2)
bacteria.soil.m$estimate<-exp(bacteria.soil.m$estimate)
bacteria.normal.m <-read.csv("unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
##visualize
bacteria.soilLV<-bacteria.soil.m[grepl("LV",bacteria.soil.m$site2),]
bacteria.soilLV<-bacteria.soilLV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="shannon"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,6,7)]
bacteria.soil<-merge(bacteria.soilLV,bacteria.normalLV,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LV"
bacteria.soil$pooling.effect<-((bacteria.soil$estimate))/(bacteria.soil$normalestimate)
bacteria.soilLV<-bacteria.soil
##LW
bacteria.soilLW<-bacteria.soil.m[grepl("LW",bacteria.soil.m$site2),]
bacteria.soilLW<-bacteria.soilLW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="shannon"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,6,7)]
bacteria.soil<-merge(bacteria.soilLW,bacteria.normalLW,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LW"
bacteria.soil$pooling.effect<-((bacteria.soil$estimate))/(bacteria.soil$normalestimate)
bacteria.soilLW<-bacteria.soil
##LZ
bacteria.soilLZ<-bacteria.soil.m[grepl("LZ",bacteria.soil.m$site2),]
bacteria.soilLZ<-bacteria.soilLZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="shannon"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,6,7)]
bacteria.soil<-merge(bacteria.soilLZ,bacteria.normalLZ,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LZ"
bacteria.soil$pooling.effect<-((bacteria.soil$estimate))/(bacteria.soil$normalestimate)
bacteria.soilLZ<-bacteria.soil
bacteria.all<-rbind(bacteria.soilLZ,bacteria.soilLV,bacteria.soilLW)
##bacteria DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.DNA.m<-read.csv("bacteria.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.DNA.m<-bacteria.DNA.m[grepl("DNA",row.names(bacteria.DNA.m)),]
row.names(bacteria.DNA.m)<-gsub("DNA","",row.names(bacteria.DNA.m))
bacteria.DNA.m$site<-row.names(bacteria.DNA.m)
bacteria.DNA.m<-merge(bacteria.DNA.m,DNA.metadata,by="site",all=T)
bacteria.DNA.m$site2<-substr(bacteria.DNA.m$site,1,2)
bacteria.DNA.m$estimate<-exp(bacteria.DNA.m$estimate)
bacteria.normal.m <-read.csv("unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
##visualize
bacteria.DNALV<-bacteria.DNA.m[grepl("LV",bacteria.DNA.m$site2),]
bacteria.DNALV<-bacteria.DNALV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="shannon"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,6,7)]
bacteria.DNA<-merge(bacteria.DNALV,bacteria.normalLV,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LV"
bacteria.DNA$pooling.effect<-((bacteria.DNA$estimate))/((bacteria.DNA$normalestimate))
bacteria.DNALV<-bacteria.DNA
##LW
bacteria.DNALW<-bacteria.DNA.m[grepl("LW",bacteria.DNA.m$site2),]
bacteria.DNALW<-bacteria.DNALW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="shannon"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,6,7)]
bacteria.DNA<-merge(bacteria.DNALW,bacteria.normalLW,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LW"
bacteria.DNA$pooling.effect<-((bacteria.DNA$estimate))/((bacteria.DNA$normalestimate))
bacteria.DNALW<-bacteria.DNA
##LZ
bacteria.DNALZ<-bacteria.DNA.m[grepl("LZ",bacteria.DNA.m$site2),]
bacteria.DNALZ<-bacteria.DNALZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="shannon"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,6,7)]
bacteria.DNA<-merge(bacteria.DNALZ,bacteria.normalLZ,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LZ"
bacteria.DNA$pooling.effect<-((bacteria.DNA$estimate))/((bacteria.DNA$normalestimate))
bacteria.DNALZ<-bacteria.DNA
bacteria.all2<-rbind(bacteria.DNALZ,bacteria.DNALV,bacteria.DNALW)
##separate the site
bacteria.all$pooling<-"soil"
bacteria.all2$pooling<-"DNA"
all<-rbind(bacteria.all,bacteria.all2)
all2<-aggregate(all[,6],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
library(ggplot2)
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
bacteria.normal.m <-read.csv("unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
##bacteria soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.soil.m<-read.csv("bacteria.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.soil.m<-bacteria.soil.m[grepl("soil",row.names(bacteria.soil.m)),]
row.names(bacteria.soil.m)<-gsub("soil","",row.names(bacteria.soil.m))
bacteria.soil.m$site<-row.names(bacteria.soil.m)

bacteria.soil.m<-merge(bacteria.soil.m,soil.metadata,by="site",all=T)
bacteria.soil.m$site2<-substr(bacteria.soil.m$site,1,2)
bacteria.soil.m$estimate<-exp(bacteria.soil.m$estimate)
##visualize
bacteria.soilLV<-bacteria.soil.m[grepl("LV",bacteria.soil.m$site2),]
bacteria.soilLV<-bacteria.soilLV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="shannon"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,6,7)]
bacteria.soil<-merge(bacteria.soilLV,bacteria.normalLV,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LV"
bacteria.soil$pooling.effect<-((bacteria.soil$estimate))/((bacteria.soil$normalestimate))
bacteria.soilLV<-bacteria.soil
##LW
bacteria.soilLW<-bacteria.soil.m[grepl("LW",bacteria.soil.m$site2),]
bacteria.soilLW<-bacteria.soilLW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="shannon"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,6,7)]
bacteria.soil<-merge(bacteria.soilLW,bacteria.normalLW,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LW"
bacteria.soil$pooling.effect<-((bacteria.soil$estimate))/((bacteria.soil$normalestimate))
bacteria.soilLW<-bacteria.soil
##LZ
bacteria.soilLZ<-bacteria.soil.m[grepl("LZ",bacteria.soil.m$site2),]
bacteria.soilLZ<-bacteria.soilLZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="shannon"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,6,7)]
bacteria.soil<-merge(bacteria.soilLZ,bacteria.normalLZ,by="method",all=T)
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$normalestimate),]
bacteria.soil<-bacteria.soil[!is.na(bacteria.soil$site),]
bacteria.soil$site2<-"LZ"
bacteria.soil$pooling.effect<-(bacteria.soil$estimate)/((bacteria.soil$normalestimate))
bacteria.soilLZ<-bacteria.soil
bacteria.all<-rbind(bacteria.soilLZ,bacteria.soilLV,bacteria.soilLW)
##bacteria DNA
bacteria.normal.m <-read.csv("unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv",header = TRUE,row.names = 1,sep = ",")
##bacteria DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
bacteria.DNA.m<-read.csv("bacteria.shannon.csv",header = TRUE,row.names = 1,sep = ",")
bacteria.DNA.m<-bacteria.DNA.m[grepl("DNA",row.names(bacteria.DNA.m)),]
row.names(bacteria.DNA.m)<-gsub("DNA","",row.names(bacteria.DNA.m))
bacteria.DNA.m$site<-row.names(bacteria.DNA.m)

bacteria.DNA.m<-merge(bacteria.DNA.m,DNA.metadata,by="site",all=T)
bacteria.DNA.m$site2<-substr(bacteria.DNA.m$site,1,2)
bacteria.DNA.m$estimate<-exp(bacteria.DNA.m$estimate)
##visualize
bacteria.DNALV<-bacteria.DNA.m[grepl("LV",bacteria.DNA.m$site2),]
bacteria.DNALV<-bacteria.DNALV[,c(1,2,8)]
bacteria.normalLV<-bacteria.normal.m[grepl("LV",bacteria.normal.m$site),]
names(bacteria.normalLV)[names(bacteria.normalLV)=="shannon"]<-"normalestimate"
names(bacteria.normalLV)[names(bacteria.normalLV)=="site"]<-"site2"
bacteria.normalLV<-bacteria.normalLV[,c(1,6,7)]
bacteria.DNA<-merge(bacteria.DNALV,bacteria.normalLV,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LV"
bacteria.DNA$pooling.effect<-(bacteria.DNA$estimate)/((bacteria.DNA$normalestimate))
bacteria.DNALV<-bacteria.DNA
##LW
bacteria.DNALW<-bacteria.DNA.m[grepl("LW",bacteria.DNA.m$site2),]
bacteria.DNALW<-bacteria.DNALW[,c(1,2,8)]
bacteria.normalLW<-bacteria.normal.m[grepl("LW",bacteria.normal.m$site),]
names(bacteria.normalLW)[names(bacteria.normalLW)=="shannon"]<-"normalestimate"
names(bacteria.normalLW)[names(bacteria.normalLW)=="site"]<-"site2"
bacteria.normalLW<-bacteria.normalLW[,c(1,6,7)]
bacteria.DNA<-merge(bacteria.DNALW,bacteria.normalLW,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LW"
bacteria.DNA$pooling.effect<-(bacteria.DNA$estimate)/((bacteria.DNA$normalestimate))
bacteria.DNALW<-bacteria.DNA
##LZ
bacteria.DNALZ<-bacteria.DNA.m[grepl("LZ",bacteria.DNA.m$site2),]
bacteria.DNALZ<-bacteria.DNALZ[,c(1,2,8)]
bacteria.normalLZ<-bacteria.normal.m[grepl("LZ",bacteria.normal.m$site),]
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="shannon"]<-"normalestimate"
names(bacteria.normalLZ)[names(bacteria.normalLZ)=="site"]<-"site2"
bacteria.normalLZ<-bacteria.normalLZ[,c(1,6,7)]
bacteria.DNA<-merge(bacteria.DNALZ,bacteria.normalLZ,by="method",all=T)
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$normalestimate),]
bacteria.DNA<-bacteria.DNA[!is.na(bacteria.DNA$site),]
bacteria.DNA$site2<-"LZ"
bacteria.DNA$pooling.effect<-(bacteria.DNA$estimate)/((bacteria.DNA$normalestimate))
bacteria.DNALZ<-bacteria.DNA

bacteria.all2<-rbind(bacteria.DNALZ,bacteria.DNALV,bacteria.DNALW)
##separate the site
bacteria.all$pooling<-"soil"
bacteria.all2$pooling<-"DNA"
all<-rbind(bacteria.all,bacteria.all2)
all2<-aggregate(all[,6],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
bacteria.normal.ratio<-all2
bacteria.normal.ratio$type<-"low"
bacteria.normal.ratio$community<-"bacteria"
bacteria.pool.ratio$type<-"high"
#summary sequencing
bacteria.summary.m<-read.csv("proportion.bacteria.shannon.atsummary sequencing depth nounrarefied sample.csv",header = TRUE,row.names = 1,sep = ",")
infor<-bacteria.summary.m[,c(1,8,13)]
infor$mer<-paste0(infor$site,infor$method)
infor<-infor[!duplicated(infor),]
bacteria.normal.ratio$mer<-paste0(bacteria.normal.ratio$site2,bacteria.normal.ratio$method)
bacteria.normal.ratio<-merge(infor[,c(2,4)],bacteria.normal.ratio,by="mer")
bacteria.normal.ratio$pooled.seqs<-1067
library(scales)
bacteria.normal.ratio$proportion<-percent(bacteria.normal.ratio$pooled.seqs/bacteria.normal.ratio$unpooledsequencing.depth)
bacteria.pool.ratio$mer<-paste0(bacteria.pool.ratio$site2,bacteria.pool.ratio$method)
bacteria.pool.ratio<-merge(infor[,c(2,4)],bacteria.pool.ratio,by="mer")
bacteria.pool.ratio$pooled.seqs<-39290
bacteria.pool.ratio$proportion<-percent(bacteria.pool.ratio$pooled.seqs/bacteria.pool.ratio$unpooledsequencing.depth)
##
bacteria.summary.m<-bacteria.summary.m[,c(1,3,5,8,11:13)]
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
animal.normal.m <- read.csv("unpooled.animal.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##animal soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.soil.m<-read.csv("animal.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.soil.m<-animal.soil.m[grepl("soil",row.names(animal.soil.m)),]
row.names(animal.soil.m)<-gsub("soil","",row.names(animal.soil.m))
animal.soil.m$site<-row.names(animal.soil.m)

animal.soil.m<-merge(animal.soil.m,soil.metadata,by="site",all=T)
animal.soil.m$site2<-substr(animal.soil.m$site,1,2)
animal.soil.m$estimate<-exp(animal.soil.m$estimate)
##visualize
animal.soilLV<-animal.soil.m[grepl("LV",animal.soil.m$site2),]
animal.soilLV<-animal.soilLV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="shannon"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,6,7)]
animal.soil<-merge(animal.soilLV,animal.normalLV,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LV"
animal.soil$pooling.effect<-(animal.soil$estimate)/((animal.soil$normalestimate))
animal.soilLV<-animal.soil
##LW
animal.soilLW<-animal.soil.m[grepl("LW",animal.soil.m$site2),]
animal.soilLW<-animal.soilLW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="shannon"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,6,7)]
animal.soil<-merge(animal.soilLW,animal.normalLW,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LW"
animal.soil$pooling.effect<-(animal.soil$estimate)/((animal.soil$normalestimate))
animal.soilLW<-animal.soil
##LZ
animal.soilLZ<-animal.soil.m[grepl("LZ",animal.soil.m$site2),]
animal.soilLZ<-animal.soilLZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="shannon"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,6,7)]
animal.soil<-merge(animal.soilLZ,animal.normalLZ,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LZ"
animal.soil$pooling.effect<-(animal.soil$estimate)/((animal.soil$normalestimate))
animal.soilLZ<-animal.soil
animal.all<-rbind(animal.soilLZ,animal.soilLV,animal.soilLW)
##animal DNA
animal.normal.m <- read.csv("unpooled.animal.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##animal DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.DNA.m<-read.csv("animal.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.DNA.m<-animal.DNA.m[grepl("DNA",row.names(animal.DNA.m)),]
row.names(animal.DNA.m)<-gsub("DNA","",row.names(animal.DNA.m))
animal.DNA.m$site<-row.names(animal.DNA.m)

animal.DNA.m<-merge(animal.DNA.m,DNA.metadata,by="site",all=T)
animal.DNA.m$site2<-substr(animal.DNA.m$site,1,2)
animal.DNA.m$estimate<-exp(animal.DNA.m$estimate)
##visualize
animal.DNALV<-animal.DNA.m[grepl("LV",animal.DNA.m$site2),]
animal.DNALV<-animal.DNALV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="shannon"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,6,7)]
animal.DNA<-merge(animal.DNALV,animal.normalLV,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LV"
animal.DNA$pooling.effect<-(animal.DNA$estimate)/((animal.DNA$normalestimate))
animal.DNALV<-animal.DNA
##LW
animal.DNALW<-animal.DNA.m[grepl("LW",animal.DNA.m$site2),]
animal.DNALW<-animal.DNALW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="shannon"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,6,7)]
animal.DNA<-merge(animal.DNALW,animal.normalLW,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LW"
animal.DNA$pooling.effect<-(animal.DNA$estimate)/((animal.DNA$normalestimate))
animal.DNALW<-animal.DNA
##LZ
animal.DNALZ<-animal.DNA.m[grepl("LZ",animal.DNA.m$site2),]
animal.DNALZ<-animal.DNALZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="shannon"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,6,7)]
animal.DNA<-merge(animal.DNALZ,animal.normalLZ,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LZ"
animal.DNA$pooling.effect<-(animal.DNA$estimate)/((animal.DNA$normalestimate))
animal.DNALZ<-animal.DNA
animal.all2<-rbind(animal.DNALZ,animal.DNALV,animal.DNALW)
##separate the site
animal.all$pooling<-"soil"
animal.all2$pooling<-"DNA"
all<-rbind(animal.all,animal.all2)
all2<-aggregate(all[,6],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
library(ggplot2)
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
animal.pool.ratio<-all2
animal.pool.ratio$community<-"animal"
###
#animal
#divide each method
animal.normal.m <- read.csv("unpooled.animal.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##animal soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.soil.m<-read.csv("animal.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.soil.m<-animal.soil.m[grepl("soil",row.names(animal.soil.m)),]
row.names(animal.soil.m)<-gsub("soil","",row.names(animal.soil.m))
animal.soil.m$site<-row.names(animal.soil.m)

animal.soil.m<-merge(animal.soil.m,soil.metadata,by="site",all=T)
animal.soil.m$site2<-substr(animal.soil.m$site,1,2)
animal.soil.m$estimate<-exp(animal.soil.m$estimate)
##visualize
animal.soilLV<-animal.soil.m[grepl("LV",animal.soil.m$site2),]
animal.soilLV<-animal.soilLV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="shannon"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,6,7)]
animal.soil<-merge(animal.soilLV,animal.normalLV,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LV"
animal.soil$pooling.effect<-(animal.soil$estimate)/(animal.soil$normalestimate)
animal.soilLV<-animal.soil
##LW
animal.soilLW<-animal.soil.m[grepl("LW",animal.soil.m$site2),]
animal.soilLW<-animal.soilLW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="shannon"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,6,7)]
animal.soil<-merge(animal.soilLW,animal.normalLW,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LW"
animal.soil$pooling.effect<-(animal.soil$estimate)/((animal.soil$normalestimate))
animal.soilLW<-animal.soil
##LZ
animal.soilLZ<-animal.soil.m[grepl("LZ",animal.soil.m$site2),]
animal.soilLZ<-animal.soilLZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="shannon"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,6,7)]
animal.soil<-merge(animal.soilLZ,animal.normalLZ,by="method",all=T)
animal.soil<-animal.soil[!is.na(animal.soil$normalestimate),]
animal.soil<-animal.soil[!is.na(animal.soil$site),]
animal.soil$site2<-"LZ"
animal.soil$pooling.effect<-(animal.soil$estimate)/((animal.soil$normalestimate))
animal.soilLZ<-animal.soil
animal.all<-rbind(animal.soilLZ,animal.soilLV,animal.soilLW)
##animal DNA
animal.normal.m <- read.csv("unpooled.animal.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
animal.DNA.m<-read.csv("animal.shannon.csv",header = TRUE,row.names = 1,sep = ",")
animal.DNA.m<-animal.DNA.m[grepl("DNA",row.names(animal.DNA.m)),]
row.names(animal.DNA.m)<-gsub("DNA","",row.names(animal.DNA.m))
animal.DNA.m$site<-row.names(animal.DNA.m)

animal.DNA.m<-merge(animal.DNA.m,DNA.metadata,by="site",all=T)
animal.DNA.m$site2<-substr(animal.DNA.m$site,1,2)
animal.DNA.m$estimate<-exp(animal.DNA.m$estimate)
##visualize
animal.DNALV<-animal.DNA.m[grepl("LV",animal.DNA.m$site2),]
animal.DNALV<-animal.DNALV[,c(1,2,8)]
animal.normalLV<-animal.normal.m[grepl("LV",animal.normal.m$site),]
names(animal.normalLV)[names(animal.normalLV)=="shannon"]<-"normalestimate"
names(animal.normalLV)[names(animal.normalLV)=="site"]<-"site2"
animal.normalLV<-animal.normalLV[,c(1,6,7)]
animal.DNA<-merge(animal.DNALV,animal.normalLV,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LV"
animal.DNA$pooling.effect<-(animal.DNA$estimate)/((animal.DNA$normalestimate))
animal.DNALV<-animal.DNA
##LW
animal.DNALW<-animal.DNA.m[grepl("LW",animal.DNA.m$site2),]
animal.DNALW<-animal.DNALW[,c(1,2,8)]
animal.normalLW<-animal.normal.m[grepl("LW",animal.normal.m$site),]
names(animal.normalLW)[names(animal.normalLW)=="shannon"]<-"normalestimate"
names(animal.normalLW)[names(animal.normalLW)=="site"]<-"site2"
animal.normalLW<-animal.normalLW[,c(1,6,7)]
animal.DNA<-merge(animal.DNALW,animal.normalLW,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LW"
animal.DNA$pooling.effect<-(animal.DNA$estimate)/((animal.DNA$normalestimate))
animal.DNALW<-animal.DNA
##LZ
animal.DNALZ<-animal.DNA.m[grepl("LZ",animal.DNA.m$site2),]
animal.DNALZ<-animal.DNALZ[,c(1,2,8)]
animal.normalLZ<-animal.normal.m[grepl("LZ",animal.normal.m$site),]
names(animal.normalLZ)[names(animal.normalLZ)=="shannon"]<-"normalestimate"
names(animal.normalLZ)[names(animal.normalLZ)=="site"]<-"site2"
animal.normalLZ<-animal.normalLZ[,c(1,6,7)]
animal.DNA<-merge(animal.DNALZ,animal.normalLZ,by="method",all=T)
animal.DNA<-animal.DNA[!is.na(animal.DNA$normalestimate),]
animal.DNA<-animal.DNA[!is.na(animal.DNA$site),]
animal.DNA$site2<-"LZ"
animal.DNA$pooling.effect<-(animal.DNA$estimate)/((animal.DNA$normalestimate))
animal.DNALZ<-animal.DNA

animal.all2<-rbind(animal.DNALZ,animal.DNALV,animal.DNALW)
##separate the site
animal.all$pooling<-"soil"
animal.all2$pooling<-"DNA"
all<-rbind(animal.all,animal.all2)
all2<-aggregate(all[,6],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
animal.normal.ratio<-all2

animal.normal.ratio$type<-"low"
animal.normal.ratio$community<-"animal"
animal.pool.ratio$type<-"high"
#summary sequencing
animal.summary.m<-read.csv("proportio.animal.shannon.atsummary sequencing depth nounrarefied sample.csv",header = TRUE,row.names = 1,sep = ",")
infor<-animal.summary.m[,c(1,8,13)]
infor$mer<-paste0(infor$site,infor$method)
infor<-infor[!duplicated(infor),]
animal.normal.ratio$mer<-paste0(animal.normal.ratio$site2,animal.normal.ratio$method)
animal.normal.ratio<-merge(infor[,c(2,4)],animal.normal.ratio,by="mer")
animal.normal.ratio$pooled.seqs<-880
library(scales)
animal.normal.ratio$proportion<-percent(animal.normal.ratio$pooled.seqs/animal.normal.ratio$unpooledsequencing.depth)
animal.pool.ratio$mer<-paste0(animal.pool.ratio$site2,animal.pool.ratio$method)
animal.pool.ratio<-merge(infor[,c(2,4)],animal.pool.ratio,by="mer")
animal.pool.ratio$pooled.seqs<-15962
animal.pool.ratio$proportion<-percent(animal.pool.ratio$pooled.seqs/animal.pool.ratio$unpooledsequencing.depth)
##
animal.summary.m<-animal.summary.m[,c(1,3,5,8,11:13)]
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
fungi.normal.m <- read.csv("unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##fungi soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.soil.m<-read.csv("fungi.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.soil.m<-fungi.soil.m[grepl("soil",row.names(fungi.soil.m)),]
row.names(fungi.soil.m)<-gsub("soil","",row.names(fungi.soil.m))
fungi.soil.m$site<-row.names(fungi.soil.m)

fungi.soil.m<-merge(fungi.soil.m,soil.metadata,by="site",all=T)
fungi.soil.m$site2<-substr(fungi.soil.m$site,1,2)
fungi.soil.m$estimate<-exp(fungi.soil.m$estimate)
##
fungi.soilLV<-fungi.soil.m[grepl("LV",fungi.soil.m$site2),]
fungi.soilLV<-fungi.soilLV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="shannon"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,6,7)]
fungi.soil<-merge(fungi.soilLV,fungi.normalLV,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LV"
fungi.soil$pooling.effect<-(fungi.soil$estimate)/((fungi.soil$normalestimate))
fungi.soilLV<-fungi.soil
##LW
fungi.soilLW<-fungi.soil.m[grepl("LW",fungi.soil.m$site2),]
fungi.soilLW<-fungi.soilLW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="shannon"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,6,7)]
fungi.soil<-merge(fungi.soilLW,fungi.normalLW,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LW"
fungi.soil$pooling.effect<-(fungi.soil$estimate)/((fungi.soil$normalestimate))
fungi.soilLW<-fungi.soil
##LZ
fungi.soilLZ<-fungi.soil.m[grepl("LZ",fungi.soil.m$site2),]
fungi.soilLZ<-fungi.soilLZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="shannon"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,6,7)]
fungi.soil<-merge(fungi.soilLZ,fungi.normalLZ,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LZ"
fungi.soil$pooling.effect<-(fungi.soil$estimate)/((fungi.soil$normalestimate))
fungi.soilLZ<-fungi.soil
fungi.all<-rbind(fungi.soilLZ,fungi.soilLV,fungi.soilLW)
##fungi DNA
fungi.normal.m <- read.csv("unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##fungi DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.DNA.m<-read.csv("fungi.pooled.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.DNA.m<-fungi.DNA.m[grepl("DNA",row.names(fungi.DNA.m)),]
row.names(fungi.DNA.m)<-gsub("DNA","",row.names(fungi.DNA.m))
fungi.DNA.m$site<-row.names(fungi.DNA.m)

fungi.DNA.m<-merge(fungi.DNA.m,DNA.metadata,by="site",all=T)
fungi.DNA.m$site2<-substr(fungi.DNA.m$site,1,2)
fungi.DNA.m$estimate<-exp(fungi.DNA.m$estimate)
##visualize
fungi.DNALV<-fungi.DNA.m[grepl("LV",fungi.DNA.m$site2),]
fungi.DNALV<-fungi.DNALV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="shannon"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,6,7)]
fungi.DNA<-merge(fungi.DNALV,fungi.normalLV,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LV"
fungi.DNA$pooling.effect<-(fungi.DNA$estimate)/((fungi.DNA$normalestimate))
fungi.DNALV<-fungi.DNA
##LW
fungi.DNALW<-fungi.DNA.m[grepl("LW",fungi.DNA.m$site2),]
fungi.DNALW<-fungi.DNALW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="shannon"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,6,7)]
fungi.DNA<-merge(fungi.DNALW,fungi.normalLW,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LW"
fungi.DNA$pooling.effect<-(fungi.DNA$estimate)/((fungi.DNA$normalestimate))
fungi.DNALW<-fungi.DNA
##LZ
fungi.DNALZ<-fungi.DNA.m[grepl("LZ",fungi.DNA.m$site2),]
fungi.DNALZ<-fungi.DNALZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="shannon"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,6,7)]
fungi.DNA<-merge(fungi.DNALZ,fungi.normalLZ,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LZ"
fungi.DNA$pooling.effect<-(fungi.DNA$estimate)/((fungi.DNA$normalestimate))
fungi.DNALZ<-fungi.DNA
fungi.all2<-rbind(fungi.DNALZ,fungi.DNALV,fungi.DNALW)
##separate the site
fungi.all$pooling<-"soil"
fungi.all2$pooling<-"DNA"
all<-rbind(fungi.all,fungi.all2)
all2<-aggregate(all[,6],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
library(ggplot2)
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
fungi.pool.ratio<-all2
fungi.pool.ratio$community<-"fungi"
###
#fungi
#divide each method
fungi.normal.m <- read.csv("unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##fungi soil
soil.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.soil.m<-read.csv("fungi.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.soil.m<-fungi.soil.m[grepl("soil",row.names(fungi.soil.m)),]
row.names(fungi.soil.m)<-gsub("soil","",row.names(fungi.soil.m))
fungi.soil.m$site<-row.names(fungi.soil.m)

fungi.soil.m<-merge(fungi.soil.m,soil.metadata,by="site",all=T)
fungi.soil.m$site2<-substr(fungi.soil.m$site,1,2)
fungi.soil.m$estimate<-exp(fungi.soil.m$estimate)
##visualize
fungi.soilLV<-fungi.soil.m[grepl("LV",fungi.soil.m$site2),]
fungi.soilLV<-fungi.soilLV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="shannon"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,6,7)]
fungi.soil<-merge(fungi.soilLV,fungi.normalLV,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LV"
fungi.soil$pooling.effect<-fungi.soil$estimate/((fungi.soil$normalestimate))
fungi.soilLV<-fungi.soil
##LW
fungi.soilLW<-fungi.soil.m[grepl("LW",fungi.soil.m$site2),]
fungi.soilLW<-fungi.soilLW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="shannon"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,6,7)]
fungi.soil<-merge(fungi.soilLW,fungi.normalLW,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LW"
fungi.soil$pooling.effect<-(fungi.soil$estimate)/((fungi.soil$normalestimate))
fungi.soilLW<-fungi.soil
##LZ
fungi.soilLZ<-fungi.soil.m[grepl("LZ",fungi.soil.m$site2),]
fungi.soilLZ<-fungi.soilLZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="shannon"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,6,7)]
fungi.soil<-merge(fungi.soilLZ,fungi.normalLZ,by="method",all=T)
fungi.soil<-fungi.soil[!is.na(fungi.soil$normalestimate),]
fungi.soil<-fungi.soil[!is.na(fungi.soil$site),]
fungi.soil$site2<-"LZ"
fungi.soil$pooling.effect<-(fungi.soil$estimate)/((fungi.soil$normalestimate))
fungi.soilLZ<-fungi.soil
fungi.all<-rbind(fungi.soilLZ,fungi.soilLV,fungi.soilLW)
##fungi DNA
fungi.normal.m <- read.csv("unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv",row.names = 1)
##fungi DNA
DNA.metadata <- read.csv("soil.metadata.csv",header = TRUE,sep = ",")
fungi.DNA.m<-read.csv("fungi.shannon.csv",header = TRUE,row.names = 1,sep = ",")
fungi.DNA.m<-fungi.DNA.m[grepl("DNA",row.names(fungi.DNA.m)),]
row.names(fungi.DNA.m)<-gsub("DNA","",row.names(fungi.DNA.m))
fungi.DNA.m$site<-row.names(fungi.DNA.m)

fungi.DNA.m<-merge(fungi.DNA.m,DNA.metadata,by="site",all=T)
fungi.DNA.m$site2<-substr(fungi.DNA.m$site,1,2)
fungi.DNA.m$estimate<-exp(fungi.DNA.m$estimate)
##visualize
fungi.DNALV<-fungi.DNA.m[grepl("LV",fungi.DNA.m$site2),]
fungi.DNALV<-fungi.DNALV[,c(1,2,8)]
fungi.normalLV<-fungi.normal.m[grepl("LV",fungi.normal.m$site),]
names(fungi.normalLV)[names(fungi.normalLV)=="shannon"]<-"normalestimate"
names(fungi.normalLV)[names(fungi.normalLV)=="site"]<-"site2"
fungi.normalLV<-fungi.normalLV[,c(1,6,7)]
fungi.DNA<-merge(fungi.DNALV,fungi.normalLV,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LV"
fungi.DNA$pooling.effect<-(fungi.DNA$estimate)/((fungi.DNA$normalestimate))
fungi.DNALV<-fungi.DNA
##LW
fungi.DNALW<-fungi.DNA.m[grepl("LW",fungi.DNA.m$site2),]
fungi.DNALW<-fungi.DNALW[,c(1,2,8)]
fungi.normalLW<-fungi.normal.m[grepl("LW",fungi.normal.m$site),]
names(fungi.normalLW)[names(fungi.normalLW)=="shannon"]<-"normalestimate"
names(fungi.normalLW)[names(fungi.normalLW)=="site"]<-"site2"
fungi.normalLW<-fungi.normalLW[,c(1,6,7)]
fungi.DNA<-merge(fungi.DNALW,fungi.normalLW,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LW"
fungi.DNA$pooling.effect<-(fungi.DNA$estimate)/((fungi.DNA$normalestimate))
fungi.DNALW<-fungi.DNA
##LZ
fungi.DNALZ<-fungi.DNA.m[grepl("LZ",fungi.DNA.m$site2),]
fungi.DNALZ<-fungi.DNALZ[,c(1,2,8)]
fungi.normalLZ<-fungi.normal.m[grepl("LZ",fungi.normal.m$site),]
names(fungi.normalLZ)[names(fungi.normalLZ)=="shannon"]<-"normalestimate"
names(fungi.normalLZ)[names(fungi.normalLZ)=="site"]<-"site2"
fungi.normalLZ<-fungi.normalLZ[,c(1,6,7)]
fungi.DNA<-merge(fungi.DNALZ,fungi.normalLZ,by="method",all=T)
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$normalestimate),]
fungi.DNA<-fungi.DNA[!is.na(fungi.DNA$site),]
fungi.DNA$site2<-"LZ"
fungi.DNA$pooling.effect<-(fungi.DNA$estimate)/((fungi.DNA$normalestimate))
fungi.DNALZ<-fungi.DNA

fungi.all2<-rbind(fungi.DNALZ,fungi.DNALV,fungi.DNALW)
##separate the site
fungi.all$pooling<-"soil"
fungi.all2$pooling<-"DNA"
all<-rbind(fungi.all,fungi.all2)
all2<-aggregate(all[,6],by=list(method=all$method,site2=all$site2,pooling=all$pooling),mean)
fungi.normal.ratio<-all2

fungi.normal.ratio$type<-"low"
fungi.normal.ratio$community<-"fungi"
fungi.pool.ratio$type<-"high"
#summary sequencing
fungi.summary.m<-read.csv("proportion.fungi.shannon.atsummary sequencing depth nounrarefied sample.csv",header = TRUE,row.names = 1,sep = ",")
infor<-fungi.summary.m[,c(1,8,13)]
infor$mer<-paste0(infor$site,infor$method)
infor<-infor[!duplicated(infor),]
fungi.normal.ratio$mer<-paste0(fungi.normal.ratio$site2,fungi.normal.ratio$method)
fungi.normal.ratio<-merge(infor[,c(2,4)],fungi.normal.ratio,by="mer")
fungi.normal.ratio$pooled.seqs<-814

fungi.normal.ratio$proportion<-percent(fungi.normal.ratio$pooled.seqs/fungi.normal.ratio$unpooledsequencing.depth)
fungi.pool.ratio$mer<-paste0(fungi.pool.ratio$site2,fungi.pool.ratio$method)
fungi.pool.ratio<-merge(infor[,c(2,4)],fungi.pool.ratio,by="mer")
fungi.pool.ratio$pooled.seqs<-14409
fungi.pool.ratio$proportion<-percent(fungi.pool.ratio$pooled.seqs/fungi.pool.ratio$unpooledsequencing.depth)
##
fungi.summary.m<-fungi.summary.m[,c(1,3,5,8,11:13)]
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
write.csv(all,"proportion.shannon.all.csv")
##animal
library(dplyr)
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

soil <- soil %>%
  filter(!(type == "high"&number==0.25))
soil <- soil %>%
  filter(!(type == "high"&number==0.5))
soil <- soil %>%
  filter(!(type == "high"&number==0.75))
soil <- soil %>%
  filter(!(type == "high"&number==1))
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = soil)
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
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = DNA)
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
soil <- soil %>%
  filter(!(type == "high"&number==0.25))
soil <- soil %>%
  filter(!(type == "high"&number==0.5))
soil <- soil %>%
  filter(!(type == "high"&number==0.75))
soil <- soil %>%
  filter(!(type == "high"&number==1))
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = soil)
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
  scale_y_continuous(breaks = c(0.75,1,1.25,1.5,1.75,2))+
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
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = DNA)
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
  scale_y_continuous(breaks = c(0.75,1,1.25,1.5,1.75,2))+
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
bacteria[bacteria$number==0.25,]$try<-"25%"
bacteria[bacteria$number>=0.50&bacteria$number<0.75,]$try<-"=50%~<75%"
bacteria[bacteria$number==0.75,]$try<-"75%"
bacteria[bacteria$number<0.95&bacteria$number>0.75,]$try<-"<75%~95%"
bacteria[bacteria$number==1,]$try<-"100%"
bacteria[bacteria$number>1,]$try<-">100%"
soil<-bacteria[grepl("soil",bacteria$pooling),]
soil<-soil[!is.na(soil$ratio),]
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = soil)
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
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
  scale_x_discrete(limits=as.character(c( "<25%","25%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  scale_y_continuous(breaks = c(0.8,1,1.5,2,2.5))+
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
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
mod01<-lmer(ratio~ try+(1|method )+ (1|site),data = DNA)
b<-avg_comparisons(mod01, variables = list(try = "pairwise")) 
letter<-cldList(p.value~contrast,threshold = 0.05,data=b)
letter$Group[1]<-"100%"
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
  scale_x_discrete(limits=as.character(c( "<25%","25%","=50%~<75%","75%","<75%~95%","100%",">100%")))+ 
  geom_hline(yintercept = 1,linetype="dashed",color="grey")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))

###beta diversity comparison across the designs
##bacteria
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal<-bacteria.normal[,!grepl("DNA|soil",names(bacteria.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
bacteria.normal <-bacteria.normal[rowSums(bacteria.normal)!=0,]
bacteria.normal <-bacteria.normal[,!grepl("LV7|LZ7|LW7",names(bacteria.normal))]
bacteria.normal <-bacteria.normal[,!grepl("LVB|LZB|LWB",names(bacteria.normal))]
###
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc_40A<-bacteria.normal[,grepl(GSMc40A,names(bacteria.normal))]

GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
GSMc_40B<-bacteria.normal[,grepl(GSMc40B,names(bacteria.normal))]
names(GSMc_40A)<-paste0(names(GSMc_40A),"GSMc40A")
names(GSMc_40B)<-paste0(names(GSMc_40B),"GSMc40B")
GSMc_40B$OTU<-row.names(GSMc_40B)
GSMc_40A$OTU<-row.names(GSMc_40A)

bacteria.normal$OTU<-row.names(bacteria.normal)
a<-list(bacteria.normal,GSMc_40A,GSMc_40B)
mer<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all<-Reduce(mer,a)
bacteria.normal<-all
##
GSMc40A<-data.frame(site=paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),"GSMc40A"),
          method="GSMc40A")
GSMc40B<-data.frame(site=paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),"GSMc40B"),
          method="GSMc40B")
library(dplyr)
metadata.normal<-select(metadata.normal, c("site","method"))
metadata2<-rbind(GSMc40A,GSMc40B,metadata.normal)

library(vegan)
row.names(bacteria.normal)<-bacteria.normal$OTU
bacteria.normal<-bacteria.normal[,-1]
bacteria.normal<-as.data.frame(t(bacteria.normal))
bacteria.normal3<-decostand(bacteria.normal,method = "total")
dist_normal<- vegdist(bacteria.normal3, method = "bray")
normal.diversity2<- metadata2[match(row.names(bacteria.normal3),metadata2$site), ]
normal.diversity2$site2<-substr(normal.diversity2$site,1,2)
##
set.seed(9999)
adonis_result<-adonis2(dist_normal~normal.diversity2$method,strata=normal.diversity2$site2,data=normal.diversity2,permutations = 999)
#using site as random factor
pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
{
  
  library(vegan)
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
  
}

b<-pairwise.adonis(bacteria.normal3, normal.diversity2$method, sim.method="bray", p.adjust.m= "bonferroni")
##fungi
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal<-fungi.normal[,!grepl("DNA|soil",names(fungi.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
fungi.normal <-fungi.normal[rowSums(fungi.normal)!=0,]
fungi.normal <-fungi.normal[,!grepl("LV7|LZ7|LW7",names(fungi.normal))]
fungi.normal <-fungi.normal[,!grepl("LVB|LZB|LWB",names(fungi.normal))]
###
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc_40A<-fungi.normal[,grepl(GSMc40A,names(fungi.normal))]

GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
GSMc_40B<-fungi.normal[,grepl(GSMc40B,names(fungi.normal))]
names(GSMc_40A)<-paste0(names(GSMc_40A),"GSMc40A")
names(GSMc_40B)<-paste0(names(GSMc_40B),"GSMc40B")
GSMc_40B$OTU<-row.names(GSMc_40B)
GSMc_40A$OTU<-row.names(GSMc_40A)

fungi.normal$OTU<-row.names(fungi.normal)
a<-list(fungi.normal,GSMc_40A,GSMc_40B)
mer<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all<-Reduce(mer,a)
fungi.normal<-all
##
GSMc40A<-data.frame(site=paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),"GSMc40A"),
                    method="GSMc40A")
GSMc40B<-data.frame(site=paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),"GSMc40B"),
                    method="GSMc40B")
library(dplyr)
metadata.normal<-select(metadata.normal, c("site","method"))
metadata2<-rbind(GSMc40A,GSMc40B,metadata.normal)

library(vegan)
row.names(fungi.normal)<-fungi.normal$OTU
fungi.normal<-fungi.normal[,-1]
fungi.normal<-as.data.frame(t(fungi.normal))
fungi.normal3<-decostand(fungi.normal,method = "total")
dist_normal<- vegdist(fungi.normal3, method = "bray")
normal.diversity2<- metadata2[match(row.names(fungi.normal3),metadata2$site), ]
normal.diversity2$site2<-substr(normal.diversity2$site,1,2)
##
set.seed(9999)
adonis_result<-adonis2(dist_normal~normal.diversity2$method,strata=normal.diversity2$site2,permutations = 999)
adonis_result

b<-pairwise.adonis(fungi.normal3, normal.diversity2$method, sim.method="bray", p.adjust.m= "bonferroni")
##animal
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
animal.normal<-animal.normal[,!grepl("DNA|soil",names(animal.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
animal.normal <-animal.normal[rowSums(animal.normal)!=0,]
animal.normal <-animal.normal[,!grepl("LV7|LZ7|LW7",names(animal.normal))]
animal.normal <-animal.normal[,!grepl("LVB|LZB|LWB",names(animal.normal))]
###
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc_40A<-animal.normal[,grepl(GSMc40A,names(animal.normal))]

GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
GSMc_40B<-animal.normal[,grepl(GSMc40B,names(animal.normal))]
names(GSMc_40A)<-paste0(names(GSMc_40A),"GSMc40A")
names(GSMc_40B)<-paste0(names(GSMc_40B),"GSMc40B")
GSMc_40B$OTU<-row.names(GSMc_40B)
GSMc_40A$OTU<-row.names(GSMc_40A)

animal.normal$OTU<-row.names(animal.normal)
a<-list(animal.normal,GSMc_40A,GSMc_40B)
mer<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all<-Reduce(mer,a)
animal.normal<-all
##
GSMc40A<-data.frame(site=paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),"GSMc40A"),
                    method="GSMc40A")
GSMc40B<-data.frame(site=paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),"GSMc40B"),
                    method="GSMc40B")
metadata.normal<-select(metadata.normal, c("site","method"))
metadata2<-rbind(GSMc40A,GSMc40B,metadata.normal)

row.names(animal.normal)<-animal.normal$OTU
animal.normal<-animal.normal[,-1]
animal.normal<-as.data.frame(t(animal.normal))
animal.normal3<-decostand(animal.normal,method = "total")
dist_normal<- vegdist(animal.normal3, method = "bray")
normal.diversity2<- metadata2[match(row.names(animal.normal3),metadata2$site), ]
normal.diversity2$site2<-substr(normal.diversity2$site,1,2)
##
set.seed(9999)
adonis_result<-adonis2(dist_normal~normal.diversity2$method,strata=normal.diversity2$site2,permutations = 999)
adonis_result

b<-pairwise.adonis(animal.normal3, normal.diversity2$method, sim.method="bray", p.adjust.m= "bonferroni")
b
###calculate community dispersion
##bacteria
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal<-bacteria.normal[,!grepl("DNA|soil",names(bacteria.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
bacteria.normal <-bacteria.normal[rowSums(bacteria.normal)!=0,]
bacteria.normal <-bacteria.normal[,!grepl("LV7|LZ7|LW7",names(bacteria.normal))]
bacteria.normal <-bacteria.normal[,!grepl("LVB|LZB|LWB",names(bacteria.normal))]
###
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc_40A<-bacteria.normal[,grepl(GSMc40A,names(bacteria.normal))]

GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
GSMc_40B<-bacteria.normal[,grepl(GSMc40B,names(bacteria.normal))]
names(GSMc_40A)<-paste0(names(GSMc_40A),"GSMc40A")
names(GSMc_40B)<-paste0(names(GSMc_40B),"GSMc40B")
GSMc_40B$OTU<-row.names(GSMc_40B)
GSMc_40A$OTU<-row.names(GSMc_40A)

bacteria.normal$OTU<-row.names(bacteria.normal)
a<-list(bacteria.normal,GSMc_40A,GSMc_40B)
mer<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all<-Reduce(mer,a)
bacteria.normal<-all
##
GSMc40A<-data.frame(site=paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),"GSMc40A"),
                    method="GSMc40A")
GSMc40B<-data.frame(site=paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),"GSMc40B"),
                    method="GSMc40B")
library(dplyr)
metadata.normal<-select(metadata.normal, c("site","method"))
metadata2<-rbind(GSMc40A,GSMc40B,metadata.normal)

library(vegan)
row.names(bacteria.normal)<-bacteria.normal$OTU
bacteria.normal<-bacteria.normal[,-1]
bacteria.normal<-as.data.frame(t(bacteria.normal))
bacteria.normal3<-decostand(bacteria.normal,method = "total")
dist_normal<- vegdist(bacteria.normal3, method = "bray")
normal.diversity2<- metadata2[match(row.names(bacteria.normal3),metadata2$site), ]
normal.diversity2$site2<-substr(normal.diversity2$site,1,2)
##LZ
a<-as.data.frame(as.matrix(dist_normal))
LZ<-a[,grepl("LZ",names(a))]
LZ<-LZ[grepl("LZ",row.names(LZ)),]
normal.diversity.LZ2<-normal.diversity2[match(row.names(LZ),normal.diversity2$site), ]
LZ<-as.dist(LZ)
set.seed(9999)
dispersion_method <- betadisper(LZ, group=normal.diversity.LZ2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

###LV
LV<-a[,grepl("LV",names(a))]
LV<-LV[grepl("LV",row.names(LV)),]
normal.diversity.LV2<-normal.diversity2[match(row.names(LV),normal.diversity2$site), ]
LV<-as.dist(LV)
set.seed(9999)
dispersion_method <- betadisper(LV, group=normal.diversity.LV2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

#LW

LW<-a[,grepl("LW",names(a))]
LW<-LW[grepl("LW",row.names(LW)),]
normal.diversity.LW2<-normal.diversity2[match(row.names(LW),normal.diversity2$site), ]
LW<-as.dist(LW)
set.seed(9999)
dispersion_method <- betadisper(LW, group=normal.diversity.LW2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

##fungi
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal<-fungi.normal[,!grepl("DNA|soil",names(fungi.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
fungi.normal <-fungi.normal[rowSums(fungi.normal)!=0,]
fungi.normal <-fungi.normal[,!grepl("LV7|LZ7|LW7",names(fungi.normal))]
fungi.normal <-fungi.normal[,!grepl("LVB|LZB|LWB",names(fungi.normal))]
###
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc_40A<-fungi.normal[,grepl(GSMc40A,names(fungi.normal))]

GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
GSMc_40B<-fungi.normal[,grepl(GSMc40B,names(fungi.normal))]
names(GSMc_40A)<-paste0(names(GSMc_40A),"GSMc40A")
names(GSMc_40B)<-paste0(names(GSMc_40B),"GSMc40B")
GSMc_40B$OTU<-row.names(GSMc_40B)
GSMc_40A$OTU<-row.names(GSMc_40A)

fungi.normal$OTU<-row.names(fungi.normal)
a<-list(fungi.normal,GSMc_40A,GSMc_40B)
mer<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all<-Reduce(mer,a)
fungi.normal<-all
##
GSMc40A<-data.frame(site=paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),"GSMc40A"),
                    method="GSMc40A")
GSMc40B<-data.frame(site=paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),"GSMc40B"),
                    method="GSMc40B")

metadata.normal<-select(metadata.normal, c("site","method"))
metadata2<-rbind(GSMc40A,GSMc40B,metadata.normal)

row.names(fungi.normal)<-fungi.normal$OTU
fungi.normal<-fungi.normal[,-1]
fungi.normal<-as.data.frame(t(fungi.normal))
fungi.normal3<-decostand(fungi.normal,method = "total")
dist_normal<- vegdist(fungi.normal3, method = "bray")
normal.diversity2<- metadata2[match(row.names(fungi.normal3),metadata2$site), ]
normal.diversity2$site2<-substr(normal.diversity2$site,1,2)
##
##LZ
a<-as.data.frame(as.matrix(dist_normal))
LZ<-a[,grepl("LZ",names(a))]
LZ<-LZ[grepl("LZ",row.names(LZ)),]
normal.diversity.LZ2<-normal.diversity2[match(row.names(LZ),normal.diversity2$site), ]
LZ<-as.dist(LZ)
set.seed(9999)
dispersion_method <- betadisper(LZ, group=normal.diversity.LZ2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

###LV
LV<-a[,grepl("LV",names(a))]
LV<-LV[grepl("LV",row.names(LV)),]
normal.diversity.LV2<-normal.diversity2[match(row.names(LV),normal.diversity2$site), ]
LV<-as.dist(LV)
set.seed(9999)
dispersion_method <- betadisper(LV, group=normal.diversity.LV2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

#LW
LW<-a[,grepl("LW",names(a))]
LW<-LW[grepl("LW",row.names(LW)),]
normal.diversity.LW2<-normal.diversity2[match(row.names(LW),normal.diversity2$site), ]
LW<-as.dist(LW)
set.seed(9999)
dispersion_method <- betadisper(LW, group=normal.diversity.LW2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

##animal
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
animal.normal<-animal.normal[,!grepl("DNA|soil",names(animal.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
animal.normal <-animal.normal[rowSums(animal.normal)!=0,]
animal.normal <-animal.normal[,!grepl("LV7|LZ7|LW7",names(animal.normal))]
animal.normal <-animal.normal[,!grepl("LVB|LZB|LWB",names(animal.normal))]
###
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc_40A<-animal.normal[,grepl(GSMc40A,names(animal.normal))]

GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
GSMc_40B<-animal.normal[,grepl(GSMc40B,names(animal.normal))]
names(GSMc_40A)<-paste0(names(GSMc_40A),"GSMc40A")
names(GSMc_40B)<-paste0(names(GSMc_40B),"GSMc40B")
GSMc_40B$OTU<-row.names(GSMc_40B)
GSMc_40A$OTU<-row.names(GSMc_40A)

animal.normal$OTU<-row.names(animal.normal)
a<-list(animal.normal,GSMc_40A,GSMc_40B)
mer<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all<-Reduce(mer,a)
animal.normal<-all
##
GSMc40A<-data.frame(site=paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X),
                                  gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X)),"GSMc40A"),
                    method="GSMc40A")
GSMc40B<-data.frame(site=paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X),
                                  gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X)),"GSMc40B"),
                    method="GSMc40B")
library(dplyr)
metadata.normal<-select(metadata.normal, c("site","method"))
metadata2<-rbind(GSMc40A,GSMc40B,metadata.normal)

library(vegan)
row.names(animal.normal)<-animal.normal$OTU
animal.normal<-animal.normal[,-1]
animal.normal<-as.data.frame(t(animal.normal))
animal.normal3<-decostand(animal.normal,method = "total")
dist_normal<- vegdist(animal.normal3, method = "bray")
normal.diversity2<- metadata2[match(row.names(animal.normal3),metadata2$site), ]
normal.diversity2$site2<-substr(normal.diversity2$site,1,2)
##LZ
a<-as.data.frame(as.matrix(dist_normal))
LZ<-a[,grepl("LZ",names(a))]
LZ<-LZ[grepl("LZ",row.names(LZ)),]
normal.diversity.LZ2<-normal.diversity2[match(row.names(LZ),normal.diversity2$site), ]
LZ<-as.dist(LZ)
set.seed(9999)
dispersion_method <- betadisper(LZ, group=normal.diversity.LZ2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)

###LV
LV<-a[,grepl("LV",names(a))]
LV<-LV[grepl("LV",row.names(LV)),]
normal.diversity.LV2<-normal.diversity2[match(row.names(LV),normal.diversity2$site), ]
LV<-as.dist(LV)
set.seed(9999)
dispersion_method <- betadisper(LV, group=normal.diversity.LV2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)
#LW

LW<-a[,grepl("LW",names(a))]
LW<-LW[grepl("LW",row.names(LW)),]
normal.diversity.LW2<-normal.diversity2[match(row.names(LW),normal.diversity2$site), ]
LW<-as.dist(LW)
set.seed(9999)
dispersion_method <- betadisper(LW, group=normal.diversity.LW2$method,
                                type = c("centroid"), bias.adjust = FALSE,
                                sqrt.dist = FALSE, add = FALSE)
###significance
set.seed(888)
anova(dispersion_method)
####comparison among pooling types in community compositional similarity between sampling designs
##bacteria
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal<-bacteria.normal[,!grepl("DNA|soil",names(bacteria.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
bacteria.normal <-bacteria.normal[rowSums(bacteria.normal)!=0,]
bacteria.normal <-bacteria.normal[,!grepl("LV7|LZ7|LW7",names(bacteria.normal))]
bacteria.normal <-bacteria.normal[,!grepl("LVB|LZB|LWB",names(bacteria.normal))]
metadata2<-metadata.normal

bacteria.normal<-as.data.frame(t(bacteria.normal))
bacteria.normal$site<-row.names(bacteria.normal)
bacteria.normal<-merge(bacteria.normal,metadata2[,c(1,2)],by="site")
bacteria.normal$site<-substr(bacteria.normal$site,1,2)
bacteria.normal<-aggregate(bacteria.normal[,2:(ncol(bacteria.normal)-1)],by=list(method=bacteria.normal$method,site=bacteria.normal$site),sum)
row.names(bacteria.normal)<-paste0(bacteria.normal$site,bacteria.normal$method)
bacteria.normal<-bacteria.normal[,-c(1:2)]
bacteria.normal3<-decostand(bacteria.normal,method = "total")
dist_normal<- vegdist(bacteria.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <-dist_normal2
dist_normal3$site1<-substr(dist_normal3$row,1,2)
dist_normal3$site2<-substr(dist_normal3$col,1,2)
between<- dist_normal3[dist_normal3$row != dist_normal3$col,]
between<- between[substr(dist_normal3$row,1,2) == substr(dist_normal3$col,1,2),]
bacteria.between<-between
bacteria.between$col<-gsub("LZ|LW|LV","",bacteria.between$col)
bacteria.between$row<-gsub("LZ|LW|LV","",bacteria.between$row)
bacteria.between<-data.frame(similarity=bacteria.between$value,site=bacteria.between$site1,comparison=paste0(bacteria.between$col,"-",bacteria.between$row),pooling="unpooled")
##animal
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
animal.normal<-animal.normal[,!grepl("DNA|soil",names(animal.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
animal.normal <-animal.normal[rowSums(animal.normal)!=0,]
animal.normal <-animal.normal[,!grepl("LV7|LZ7|LW7",names(animal.normal))]
animal.normal <-animal.normal[,!grepl("LVB|LZB|LWB",names(animal.normal))]
metadata2<-metadata.normal
#
animal.normal<-as.data.frame(t(animal.normal))
animal.normal$site<-row.names(animal.normal)
animal.normal<-merge(animal.normal,metadata2[,c(1,2)],by="site")
animal.normal$site<-substr(animal.normal$site,1,2)
animal.normal<-aggregate(animal.normal[,2:(ncol(animal.normal)-1)],by=list(method=animal.normal$method,site=animal.normal$site),sum)
row.names(animal.normal)<-paste0(animal.normal$site,animal.normal$method)
animal.normal<-animal.normal[,-c(1:2)]
animal.normal3<-decostand(animal.normal,method = "total")
dist_normal<- vegdist(animal.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <-dist_normal2
dist_normal3$site1<-substr(dist_normal3$row,1,2)
dist_normal3$site2<-substr(dist_normal3$col,1,2)
between<- dist_normal3[dist_normal3$row != dist_normal3$col,]
between<- between[substr(dist_normal3$row,1,2) == substr(dist_normal3$col,1,2),]
animal.between<-between
animal.between$col<-gsub("LZ|LW|LV","",animal.between$col)
animal.between$row<-gsub("LZ|LW|LV","",animal.between$row)
animal.between<-data.frame(similarity=animal.between$value,site=animal.between$site1,comparison=paste0(animal.between$col,"-",animal.between$row),pooling="unpooled")

##fungi
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal<-fungi.normal[,!grepl("DNA|soil",names(fungi.normal))]
metadata.normal<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata.normal$type<-"unpooled"
fungi.normal <-fungi.normal[rowSums(fungi.normal)!=0,]
fungi.normal <-fungi.normal[,!grepl("LV7|LZ7|LW7",names(fungi.normal))]
fungi.normal <-fungi.normal[,!grepl("LVB|LZB|LWB",names(fungi.normal))]
metadata2<-metadata.normal
###
fungi.normal<-as.data.frame(t(fungi.normal))
fungi.normal$site<-row.names(fungi.normal)
fungi.normal<-merge(fungi.normal,metadata2[,c(1,2)],by="site")
fungi.normal$site<-substr(fungi.normal$site,1,2)
fungi.normal<-aggregate(fungi.normal[,2:(ncol(fungi.normal)-1)],by=list(method=fungi.normal$method,site=fungi.normal$site),sum)
row.names(fungi.normal)<-paste0(fungi.normal$site,fungi.normal$method)
fungi.normal<-fungi.normal[,-c(1:2)]
fungi.normal3<-decostand(fungi.normal,method = "total")
dist_normal<- vegdist(fungi.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <-dist_normal2
dist_normal3$site1<-substr(dist_normal3$row,1,2)
dist_normal3$site2<-substr(dist_normal3$col,1,2)
between<- dist_normal3[dist_normal3$row != dist_normal3$col,]
between<- between[substr(dist_normal3$row,1,2) == substr(dist_normal3$col,1,2),]
fungi.between<-between
fungi.between$col<-gsub("LZ|LW|LV","",fungi.between$col)
fungi.between$row<-gsub("LZ|LW|LV","",fungi.between$row)
fungi.between<-data.frame(similarity=fungi.between$value,site=fungi.between$site1,comparison=paste0(fungi.between$col,"-",fungi.between$row),pooling="unpooled")
###for pool
##bacteria
library(dplyr)
#DNA
DNA<-read.table("bacteria.poolat100.csv",header = T,row.names = 1,sep = ",")
DNA<-DNA[,grepl("DNA",names(DNA))]
names(DNA)<-gsub("DNA","",names(DNA))
metadata.DNA<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata.DNA[,c(1,5)]
bacteria.normal<-as.data.frame(t(DNA))
bacteria.normal3<-decostand(bacteria.normal,method = "total")
dist_normal<- vegdist(bacteria.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <- merge(dist_normal2,metadata.DNA,by.x="col",by.y="site",all=T)
names(dist_normal3)[4] <- "method.col"
dist_normal3 <- merge(dist_normal3,metadata.DNA,by.x="row",by.y="site",all=T)
names(dist_normal3)[5] <- "method.row"
dist_normal3 <- dist_normal3[!is.na(dist_normal3$value),]
between<- dist_normal3[dist_normal3$method.row != dist_normal3$method.col,]
between$site1<-substr(between$row,1,2)
between$site2<-substr(between$col,1,2)
between<- between[between$site2 == between$site1,]
between2 <- data.frame(constract=paste0(between$method.col,"-",between$method.row),value=between$value,site=between$site2)
DNA.bacteria.between<-between2
#soil
soil<-read.table("bacteria.poolat100.csv",header = T,row.names = 1,sep = ",")
soil<-soil[,grepl("soil",names(soil))]
names(soil)<-gsub("soil","",names(soil))
metadata.soil<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.soil<-metadata.soil[,c(1,5)]
bacteria.normal<-as.data.frame(t(soil))
bacteria.normal3<-decostand(bacteria.normal,method = "total")
dist_normal<- vegdist(bacteria.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <- merge(dist_normal2,metadata.soil,by.x="col",by.y="site",all=T)
names(dist_normal3)[4] <- "method.col"
dist_normal3 <- merge(dist_normal3,metadata.soil,by.x="row",by.y="site",all=T)
names(dist_normal3)[5] <- "method.row"
dist_normal3 <- dist_normal3[!is.na(dist_normal3$value),]
between<- dist_normal3[dist_normal3$method.row != dist_normal3$method.col,]

between$site1<-substr(between$row,1,2)
between$site2<-substr(between$col,1,2)
between<- between[between$site2 == between$site1,]
between2 <- data.frame(constract=paste0(between$method.col,"-",between$method.row),value=between$value,site=between$site2)
soil.bacteria.between<-between2
soil.bacteria.between$pooling<-"soil"
DNA.bacteria.between$pooling<-"DNA"
pool.bacteria.between<-rbind(soil.bacteria.between,DNA.bacteria.between)
names(pool.bacteria.between)<-c("comparison","similarity", "site","pooling")
##animal
#DNA
DNA<-read.table("animal.poolat100.csv",header = T,row.names = 1,sep = ",")
DNA<-DNA[,grepl("DNA",names(DNA))]
names(DNA)<-gsub("DNA","",names(DNA))
metadata.DNA<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata.DNA[,c(1,5)]
animal.normal<-as.data.frame(t(DNA))
animal.normal3<-decostand(animal.normal,method = "total")
dist_normal<- vegdist(animal.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <- merge(dist_normal2,metadata.DNA,by.x="col",by.y="site",all=T)
names(dist_normal3)[4] <- "method.col"
dist_normal3 <- merge(dist_normal3,metadata.DNA,by.x="row",by.y="site",all=T)
names(dist_normal3)[5] <- "method.row"
dist_normal3 <- dist_normal3[!is.na(dist_normal3$value),]
within<- dist_normal3[dist_normal3$method.row == dist_normal3$method.col,]
between<- dist_normal3[dist_normal3$method.row != dist_normal3$method.col,]

between$site1<-substr(between$row,1,2)
between$site2<-substr(between$col,1,2)
between<- between[between$site2 == between$site1,]
between2 <- data.frame(constract=paste0(between$method.col,"-",between$method.row),value=between$value,site=between$site2)
DNA.animal.between<-between2
#soil
soil<-read.table("animal.poolat100.csv",header = T,row.names = 1,sep = ",")
soil<-soil[,grepl("soil",names(soil))]
names(soil)<-gsub("soil","",names(soil))
metadata.soil<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.soil<-metadata.soil[,c(1,5)]
animal.normal<-as.data.frame(t(soil))
animal.normal3<-decostand(animal.normal,method = "total")
dist_normal<- vegdist(animal.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <- merge(dist_normal2,metadata.soil,by.x="col",by.y="site",all=T)
names(dist_normal3)[4] <- "method.col"
dist_normal3 <- merge(dist_normal3,metadata.soil,by.x="row",by.y="site",all=T)
names(dist_normal3)[5] <- "method.row"
dist_normal3 <- dist_normal3[!is.na(dist_normal3$value),]

between<- dist_normal3[dist_normal3$method.row != dist_normal3$method.col,]
between$site1<-substr(between$row,1,2)
between$site2<-substr(between$col,1,2)
between<- between[between$site2 == between$site1,]

between2 <- data.frame(constract=paste0(between$method.col,"-",between$method.row),value=between$value,site=between$site2)
soil.animal.between<-between2
soil.animal.between$pooling<-"soil"
DNA.animal.between$pooling<-"DNA"
pool.animal.between<-rbind(soil.animal.between,DNA.animal.between)
names(pool.animal.between)<-c("comparison","similarity", "site","pooling")
##fungi
#DNA
DNA<-read.table("fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
DNA<-DNA[,grepl("DNA",names(DNA))]
names(DNA)<-gsub("DNA","",names(DNA))
metadata.DNA<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata.DNA[,c(1,5)]
fungi.normal<-as.data.frame(t(DNA))
fungi.normal3<-decostand(fungi.normal,method = "total")
dist_normal<- vegdist(fungi.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <- merge(dist_normal2,metadata.DNA,by.x="col",by.y="site",all=T)
names(dist_normal3)[4] <- "method.col"
dist_normal3 <- merge(dist_normal3,metadata.DNA,by.x="row",by.y="site",all=T)
names(dist_normal3)[5] <- "method.row"
dist_normal3 <- dist_normal3[!is.na(dist_normal3$value),]
between<- dist_normal3[dist_normal3$method.row != dist_normal3$method.col,]

between$site1<-substr(between$row,1,2)
between$site2<-substr(between$col,1,2)
between<- between[between$site2 == between$site1,]

between2 <- data.frame(constract=paste0(between$method.col,"-",between$method.row),value=between$value,site=between$site2)
DNA.fungi.between<-between2

#soil
soil<-read.table("fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
soil<-soil[,grepl("soil",names(soil))]
names(soil)<-gsub("soil","",names(soil))
metadata.soil<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.soil<-metadata.soil[,c(1,5)]
fungi.normal<-as.data.frame(t(soil))
fungi.normal3<-decostand(fungi.normal,method = "total")
dist_normal<- vegdist(fungi.normal3, method = "bray")
##
dist_normal2 <- dist2list(dist_normal)
dist_normal2$value <- 1-dist_normal2$value
##
dist_normal3 <- merge(dist_normal2,metadata.soil,by.x="col",by.y="site",all=T)
names(dist_normal3)[4] <- "method.col"
dist_normal3 <- merge(dist_normal3,metadata.soil,by.x="row",by.y="site",all=T)
names(dist_normal3)[5] <- "method.row"
dist_normal3 <- dist_normal3[!is.na(dist_normal3$value),]
between<- dist_normal3[dist_normal3$method.row != dist_normal3$method.col,]

between$site1<-substr(between$row,1,2)
between$site2<-substr(between$col,1,2)
between<- between[between$site2 == between$site1,]
between2 <- data.frame(constract=paste0(between$method.col,"-",between$method.row),value=between$value,site=between$site2)
soil.fungi.between<-between2

soil.fungi.between$pooling<-"soil"
DNA.fungi.between$pooling<-"DNA"
pool.fungi.between<-rbind(soil.fungi.between,DNA.fungi.between)
names(pool.fungi.between)<-c("comparison","similarity", "site","pooling")
###animal
animal<-rbind(pool.animal.between,animal.between)
mod01<-lmer(similarity ~  pooling + (1|site),data = animal)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###bacteria
bacteria<-rbind(pool.bacteria.between,bacteria.between)
mod01<-lmer(similarity ~  pooling + (1|site),data = bacteria)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 

###fungi
fungi<-rbind(pool.fungi.between,fungi.between)
mod01<-lmer(similarity ~  pooling + (1|site),data = fungi)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###mantel correlogram
animal<- read.table("animal.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
animal2<-animal[,!grepl("DNA|soil",names(animal))]
animal2<-as.data.frame(t(animal2))
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
animal2<-animal2[grepl(paste0(metadata[grepl("GSMc|SUCC|DarkDiv",metadata$method),]$site,collapse = "|"),row.names(animal2)),]
LV<-animal2[grepl("LV",row.names(animal2)),]
LZ<-animal2[grepl("LZ",row.names(animal2)),]
LW<-animal2[grepl("LW",row.names(animal2)),]
LZ<-LZ[rowSums(LZ)>0,]
coordinate<-read.table('coordinate.csv',sep=",",header = T)
coordinate.LZ<-coordinate[match(row.names(LZ),coordinate$sample),]
coordinate$sample<-gsub("LZ","LV",coordinate$sample)
coordinate.LV<-coordinate[match(row.names(LV),coordinate$sample),]
coordinate$sample<-gsub("LV","LW",coordinate$sample)
coordinate.LW<-coordinate[match(row.names(LW),coordinate$sample),]
#LV
set.seed(888)
LV.D <- dist(LV)
LV.correlog <- mantel.correlog(LV.D, XY=coordinate.LV[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LV.correlog)
summary(LV.correlog)
LV.correlog 
LV.correlog[["break.pts"]]
#LW
LW.D <- dist(LW)
set.seed(888)
LW.correlog <- mantel.correlog(LW.D, XY=coordinate.LW[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LW.correlog)
summary(LW.correlog)
LW.correlog  
LW.correlog[["break.pts"]]
#LZ
LZ.D <- dist(LZ)
set.seed(888)
LZ.correlog <- mantel.correlog(LZ.D, XY=coordinate.LZ[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LZ.correlog)
summary(LZ.correlog)
LZ.correlog  
LZ.correlog[["break.pts"]]
##bacteria
bacteria<- read.table("bacteria.rarefy.table.csv",header = T,row.names = 1,sep = ",")
bacteria2<-bacteria[,!grepl("DNA|soil",names(bacteria))]
bacteria2<-as.data.frame(t(bacteria2))
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
bacteria2<-bacteria2[grepl(paste0(metadata[grepl("GSMc|SUCC|DarkDiv",metadata$method),]$site,collapse = "|"),row.names(bacteria2)),]
LV<-bacteria2[grepl("LV",row.names(bacteria2)),]
LZ<-bacteria2[grepl("LZ",row.names(bacteria2)),]
LW<-bacteria2[grepl("LW",row.names(bacteria2)),]
LZ<-LZ[rowSums(LZ)>0,]

coordinate<-read.table('coordinate.csv',sep=",",header = T)
coordinate.LZ<-coordinate[match(row.names(LZ),coordinate$sample),]
coordinate$sample<-gsub("LZ","LV",coordinate$sample)
coordinate.LV<-coordinate[match(row.names(LV),coordinate$sample),]
coordinate$sample<-gsub("LV","LW",coordinate$sample)
coordinate.LW<-coordinate[match(row.names(LW),coordinate$sample),]
#LV
LV.D <- dist(LV)
set.seed(888)
LV.correlog <- mantel.correlog(LV.D, XY=coordinate.LV[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LV.correlog)
summary(LV.correlog)
LV.correlog  
LV.correlog[["break.pts"]]
#LW
library(vegan)
LW.D <- dist(LW)
set.seed(888)
LW.correlog <- mantel.correlog(LW.D, XY=coordinate.LW[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LW.correlog)
summary(LW.correlog)
LW.correlog 
LW.correlog[["break.pts"]]
#LZ
library(vegan)
LZ.D <- dist(LZ)
set.seed(888)
LZ.correlog <- mantel.correlog(LZ.D, XY=coordinate.LZ[,2:3],n.class=9,cutoff = F,nperm=999)
plot(LZ.correlog)
summary(LZ.correlog)
LZ.correlog  
LZ.correlog[["break.pts"]]
#fungi
fungi<- read.table("fungi.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
fungi2<-fungi[,!grepl("DNA|soil",names(fungi))]
fungi2<-as.data.frame(t(fungi2))
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
fungi2<-fungi2[grepl(paste0(metadata[grepl("GSMc|SUCC|DarkDiv",metadata$method),]$site,collapse = "|"),row.names(fungi2)),]
LV<-fungi2[grepl("LV",row.names(fungi2)),]
LZ<-fungi2[grepl("LZ",row.names(fungi2)),]
LW<-fungi2[grepl("LW",row.names(fungi2)),]
LZ<-LZ[rowSums(LZ)>0,]

coordinate<-read.table('coordinate.csv',sep=",",header = T)
coordinate.LZ<-coordinate[match(row.names(LZ),coordinate$sample),]
coordinate$sample<-gsub("LZ","LV",coordinate$sample)
coordinate.LV<-coordinate[match(row.names(LV),coordinate$sample),]
coordinate$sample<-gsub("LV","LW",coordinate$sample)
coordinate.LW<-coordinate[match(row.names(LW),coordinate$sample),]
#LV
LV.D <- dist(LV)
set.seed(888)
LV.correlog <- mantel.correlog(LV.D, XY=coordinate.LV[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LV.correlog)
summary(LV.correlog)
LV.correlog  
LV.correlog[["break.pts"]]
#LW
LW.D <- dist(LW)
set.seed(888)
LW.correlog <- mantel.correlog(LW.D, XY=coordinate.LW[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LW.correlog)
summary(LW.correlog)
LW.correlog  
LV.correlog[["break.pts"]]
#LZ
LZ.D <- dist(LZ)
set.seed(888)
LZ.correlog <- mantel.correlog(LZ.D, XY=coordinate.LZ[,2:3],n.class=9,cutoff = F, nperm=999)
plot(LZ.correlog)
summary(LZ.correlog)
LZ.correlog  
LZ.correlog[["break.pts"]]
###analysis about artefact
###artefact from three families
a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
b<-read.csv("otu.select.txt",header=F)
a<-a[!duplicated(a$V1),]
a1<-a[a$V1 %in% b$V1,]
a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
###
normal<-read.csv("fungi.rarefy.table2.csv")
n.r<-normal[row.names(normal) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[,!grepl("DNA|soil",names(n.r))]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[rowSums(n.r)>0,]
n.r<-n.r[,!grepl("B",names(n.r))]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(350:353)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-1)],by=list(design=nr$method,site=nr$site),sum)
nr <- nr[, c(rep(TRUE,2),(colSums(nr[,-c(1:2)]) != 0))]
arichness<-data.frame(richness=rowSums(nr[,-c(1:2)]>0),design=paste0(nr$design,nr$site))
###GSMc40AB
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-n.r[grepl(GSMc40A,row.names(n.r)),]
A40$site<-substr(row.names(A40),1,2)
A40<-aggregate(A40[,1:(ncol(A40)-1)],by=list(site=A40$site),sum)
arichness1<-data.frame(richness=rowSums(A40[,-1]>0),design=paste0("GSMC_40A",A40$site))

B40<-n.r[grepl(GSMc40B,row.names(n.r)),]
B40$site<-substr(row.names(B40),1,2)
B40<-aggregate(B40[,1:(ncol(B40)-1)],by=list(site=B40$site),sum)
arichness2<-data.frame(richness=rowSums(B40[,-1]>0),design=paste0("GSMC_40B",B40$site))
arichness<-rbind(arichness,arichness1,arichness2)
###get the artifect ratio
ar<-read.table("artefact.txt")
ar1<-read.table("thelephoraceae.artfect.txt")
ar2<-read.table("sebacinaceae.artfect.txt")
ar<-rbind(ar,ar1,ar2)
ar<-data.frame(V1=gsub("_.*","",ar$V1))
ar<-ar[!grepl(";k",ar$V1),]
nr.no<-nr[,names(nr) %in% c(ar,names(nr)[1:2])]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:2)]>0),design=paste0(nr.no$design,nr.no$site))
###GSMc40AB
A40<-n.r[grepl(GSMc40A,row.names(n.r)),]
A40$site<-substr(row.names(A40),1,2)
A40.no<-A40[,names(A40) %in% c(ar,names(A40)[ncol(A40)])]
A40<-aggregate(A40.no[,1:(ncol(A40.no)-1)],by=list(site=A40.no$site),sum)
nrichness1<-data.frame(richness=rowSums(A40[,-1]>0),design=paste0("GSMC_40A",A40$site))

B40<-n.r[grepl(GSMc40B,row.names(n.r)),]
B40$site<-substr(row.names(B40),1,2)
B40.no<-B40[,names(B40) %in% c(ar,names(B40)[ncol(B40)])]
B40<-aggregate(B40.no[,1:(ncol(B40.no)-1)],by=list(site=B40.no$site),sum)
nrichness2<-data.frame(richness=rowSums(B40[,-1]>0),design=paste0("GSMC_40B",B40$site))
nrichness<-rbind(nrichness,nrichness1,nrichness2)
###
deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ<-aggregate(deep_succ[,1:(ncol(deep_succ)-1)],by=list(site=deep_succ$site),sum)
nrichness1<-data.frame(richness=rowSums(deep_succ[,-1]>0),design=paste0("deep_SUCC",deep_succ$site))

deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ<-deep_succ[,names(deep_succ) %in% c(ar,names(deep_succ)[ncol(deep_succ)])]
deep_succ<-aggregate(deep_succ[,1:(ncol(deep_succ)-1)],by=list(site=deep_succ$site),sum)
arichness1<-data.frame(richness=rowSums(deep_succ[,-1]>0),design=paste0("deep_SUCC",deep_succ$site))
arichness<-rbind(arichness,nrichness1)
nrichness<-rbind(nrichness,arichness1)

###proportion of artifect
p.a<-merge(nrichness,arichness,by="design")
p.a$ratio<-1-(p.a$richness.x/p.a$richness.y)
###pool
pool<-read.table("fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
row.names(pool)<-pool$OTU
pool<-pool[,-1]
n.r<-pool[row.names(pool) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[rowSums(n.r)>0,]
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata$site<-paste0(metadata$site,"DNA")
metadata$type<-"DNA"
metadata2<-metadata
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata$site<-paste0(metadata$site,"soil")
metadata$type<-"soil"
metadata<-rbind(metadata,metadata2)
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(501:503)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-2)],by=list(design=nr$method,site=nr$site,type=nr$type),sum)
nr <- nr[, c(rep(TRUE,3),(colSums(nr[,-c(1:3)]) != 0))]
arichness<-data.frame(richness=rowSums(nr[,-c(1:3)]>0),design=paste0(nr$design,nr$site,nr$type))
###get the artifect ratio
nr.no<-nr[,names(nr) %in% c(ar,names(nr)[1:3])]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:3)]>0),design=paste0(nr.no$design,nr.no$site,nr$type))
###proportion of artifect
pool.p.a<-merge(nrichness,arichness,by="design",all=T)
pool.p.a$ratio<-1-(pool.p.a$richness.x/pool.p.a$richness.y)

###
pool.p.a$type<-"pooling"
p.a$type<-"unpooled"
a.r<-rbind(pool.p.a,p.a)
a.r2 <- a.r %>%
  mutate(
    site = str_extract(design, "LV|LW|LZ"), 
    pooling = str_extract(design, "DNA|soil")   ,
    design = str_remove_all(design, "LV|LW|LZ|DNA|soil")
  )
a.r2$pooling[is.na(a.r2$pooling)]<-"unpooled"
###
a.r2<-a.r2[!a.r2$design %in% c("Maestre","MDB15","MDB5"),]
a.r2$design[a.r2$design=="Deep"]<-"deep"
a.r2$design[a.r2$design=="GSMc"]<-"GSMc_62"
mod01<-lmer(ratio~ pooling + (1|site) + (1|design),data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
##group into >20 '>20' means PE<1 
a.r2$sub<-NA
a.r2$sub[grepl("GSMc|GSMC",a.r2$method)]<- ">20"
a.r2$sub[!grepl("GSMc|GSMC",a.r2$method)]<- "<20"
mod01<-lmer(ratio~ (1|site)+ (1|design)+ pooling ,data = a.r2[!grepl("<20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
mod01<-lmer(ratio~ (1|site)+ (1|design)+pooling ,data = a.r2[!grepl(">20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
######find rare species
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
b<-read.csv("otu.select.txt",header=F)
a<-a[!duplicated(a$V1),]
a1<-a[a$V1 %in% b$V1,]
a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
###
normal<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/diversity original_not summarize/fungi.rarefy.table2.csv")
n.r<-normal[row.names(normal) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[,!grepl("DNA|soil",names(n.r))]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[rowSums(n.r)>0,]
n.r<-n.r[,!grepl("B",names(n.r))]
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/normal_metadata.csv",header = TRUE,sep = ",")
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(350:353)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-1)],by=list(design=nr$method,site=nr$site),sum)
nr <- nr[, c(rep(TRUE,2),(colSums(nr[,-c(1:2)]) != 0))]
arichness<-data.frame(richness=rowSums(nr[,-c(1:2)]>0),design=paste0(nr$design,nr$site))
###GSMc40AB
select<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/1.accumulative/sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-n.r[grepl(GSMc40A,row.names(n.r)),]
A40$site<-substr(row.names(A40),1,2)
A40<-aggregate(A40[,1:(ncol(A40)-1)],by=list(site=A40$site),sum)
arichness1<-data.frame(richness=rowSums(A40[,-1]>0),design=paste0("GSMC_40A",A40$site))

B40<-n.r[grepl(GSMc40B,row.names(n.r)),]
B40$site<-substr(row.names(B40),1,2)
B40<-aggregate(B40[,1:(ncol(B40)-1)],by=list(site=B40$site),sum)
arichness2<-data.frame(richness=rowSums(B40[,-1]>0),design=paste0("GSMC_40B",B40$site))
arichness<-rbind(arichness,arichness1,arichness2)
###get rare OTU ratio
a.r<-names(nr[,c(rep(TRUE,2),(colSums(nr[,-c(1:2)])/sum(nr[,-c(1:2)]))*100 <= 0.05)])
nr.no<-nr[,names(nr) %in% c(a.r)]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:2)]>0),design=paste0(nr.no$design,nr.no$site))
###GSMc40AB
A40<-normal[,grepl(GSMc40A,names(normal))]
A40<-as.data.frame(t(A40))
A40$site<-substr(row.names(A40),1,2)
A40.no<-A40[,names(A40) %in% c(a.r,names(A40)[ncol(A40)])]
A40.no<-aggregate(A40.no[,1:(ncol(A40.no)-1)],by=list(site=A40.no$site),sum)
nrichness1<-data.frame(richness=rowSums(A40.no[,-1]>0),design=paste0("GSMC_40A",A40.no$site))

B40<-normal[,grepl(GSMc40B,names(normal))]
B40<-as.data.frame(t(B40))
B40$site<-substr(row.names(B40),1,2)
B40.no<-B40[,names(B40) %in% c(a.r,names(B40)[ncol(B40)])]
B40<-aggregate(B40.no[,1:(ncol(B40.no)-1)],by=list(site=B40.no$site),sum)
nrichness2<-data.frame(richness=rowSums(B40[,-1]>0),design=paste0("GSMC_40B",B40$site))
nrichness<-rbind(nrichness,nrichness1,nrichness2)
####deep_succ
deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ<-aggregate(deep_succ[,1:(ncol(deep_succ)-1)],by=list(site=deep_succ$site),sum)
arichness1<-data.frame(richness=rowSums(deep_succ[,-1]>0),design=paste0("deep_SUCC",deep_succ$site))

deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ<-deep_succ[,names(deep_succ) %in% c(a.r,names(deep_succ)[ncol(deep_succ)])]
deep_succ<-aggregate(deep_succ[,1:(ncol(deep_succ)-1)],by=list(site=deep_succ$site),sum)
nrichness1<-data.frame(richness=rowSums(deep_succ[,-1]>0),design=paste0("deep_SUCC",deep_succ$site))
arichness<-rbind(arichness,arichness1)
nrichness<-rbind(nrichness,nrichness1)
###proportion of artifect
p.a<-merge(nrichness,arichness,by="design")
p.a$ratio<-(p.a$richness.x/p.a$richness.y)
###pool
pool<-read.table("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/17.pool similarity/fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
row.names(pool)<-pool$OTU
pool<-pool[,-1]
n.r<-pool[row.names(pool) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata$site<-paste0(metadata$site,"DNA")
metadata$type<-"DNA"
metadata2<-metadata
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata$site<-paste0(metadata$site,"soil")
metadata$type<-"soil"
metadata<-rbind(metadata,metadata2)
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(1095:1097)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-2)],by=list(design=nr$method,site=nr$site,type=nr$type),sum)
nr <- nr[, c(rep(TRUE,3),(colSums(nr[,-c(1:3)]) != 0))]
arichness<-data.frame(richness=rowSums(nr[,-c(1:3)]>0),design=paste0(nr$design,nr$site,nr$type))
###get the artifect ratio
a.r<-names(nr[,c(rep(TRUE,3),(colSums(nr[,-c(1:3)])/sum(nr[,-c(1:3)]))*100 <= 0.05)])
nr.no<-nr[,names(nr) %in% c(a.r)]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:3)]>0),design=paste0(nr.no$design,nr.no$site,nr$type))
###proportion of artifect
pool.p.a<-merge(nrichness,arichness,by="design",all=T)
pool.p.a$ratio<-(pool.p.a$richness.x/pool.p.a$richness.y)
write.csv(pool.p.a,"rare.three.families.pool.p.a.csv")
write.csv(p.a,"rare.three.families.a.csv")
###
pool.p.a$type<-"pooling"
p.a$type<-"unpooled"
a.r<-rbind(pool.p.a,p.a)
a.r2 <- a.r %>%
  mutate(
    site = str_extract(design, "LV|LW|LZ"), 
    pooling = str_extract(design, "DNA|soil")   ,
    design = str_remove_all(design, "LV|LW|LZ|DNA|soil")
  )
a.r2$pooling[is.na(a.r2$pooling)]<-"unpooled"
###
a.r2<-a.r2[!a.r2$design %in% c("Maestre","MDB15","MDB5"),]
a.r2$design[a.r2$design=="Deep"]<-"deep"
a.r2$design[a.r2$design=="GSMc"]<-"GSMc_62"
mod01<-lmer(ratio~ pooling + (1|site) + (1|design),data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 

mod01<-lmer(richness.x~ pooling + (1|site) + (1|design),data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 

##group into >20
a.r2$sub<-NA
a.r2$sub[grepl("GSMc|GSMC",a.r2$design)] <- ">20"
a.r2$sub[!grepl("GSMc|GSMC",a.r2$design)] <- "<20"
a.r2$sub[a.r2$pooling=="soil" & a.r2$design=="Zobel"] <- ">20"
a.r2$sub[a.r2$pooling=="soil" & a.r2$design=="SUCC"] <- ">20"
mod01<-lmer(richness.x~ (1|site)+ (1|design)+ pooling ,data = a.r2[grepl("<20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
mod01<-lmer(richness.x~ (1|site)+ (1|design)+pooling ,data = a.r2[grepl(">20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###find the rare artefact
a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
b<-read.csv("otu.select.txt",header=F)
a<-a[!duplicated(a$V1),]
a1<-a[a$V1 %in% b$V1,]
a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
###get rare OTU ratio
normal<-read.csv("fungi.rarefy.table2.csv")
n.r<-normal[row.names(normal) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[,!grepl("DNA|soil",names(n.r))]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[rowSums(n.r)>0,]
n.r<-n.r[,!grepl("B",names(n.r))]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(350:353)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-1)],by=list(design=nr$method,site=nr$site),sum)
a.r<-names(nr[,c(rep(TRUE,2),(colSums(nr[,-c(1:2)])/sum(nr[,-c(1:2)]))*100 <= 0.05)])
nr.no<-nr[,names(nr) %in% c(a.r)]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:2)]>0),design=paste0(nr.no$design,nr.no$site))
##
a.r2<-read.table("thelephoraceae.txt")
a.r3<-read.table("russula.txt")
a.r4<-read.table("sebacinaceae.txt")
a.r.2<-rbind(a.r2,a.r3,a.r4)
a.r2<-a.r.2[!grepl("f__",a.r.2$V1),]
artefact<-(nr.no[,c(rep(TRUE,2),!names(nr.no[,3:ncol(nr.no)]) %in% a.r2)])
arichness<-data.frame(richness=rowSums(artefact[,-c(1:2)]>0),design=paste0(artefact$design,artefact$site))
###GSMc40AB
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-normal[,grepl(GSMc40A,names(normal))]
A40<-as.data.frame(t(A40))
A40$site<-substr(row.names(A40),1,2)
A40.no<-A40[,names(A40) %in% c(a.r,names(A40)[ncol(A40)])]
A40.no<-aggregate(A40.no[,1:(ncol(A40.no)-1)],by=list(site=A40.no$site),sum)
nrichness1<-data.frame(richness=rowSums(A40.no[,-1]>0),design=paste0("GSMC_40A",A40.no$site))

B40<-normal[,grepl(GSMc40B,names(normal))]
B40<-as.data.frame(t(B40))
B40$site<-substr(row.names(B40),1,2)
B40.no<-B40[,names(B40) %in% c(a.r,names(B40)[ncol(B40)])]
B40.no<-aggregate(B40.no[,1:(ncol(B40.no)-1)],by=list(site=B40.no$site),sum)
nrichness2<-data.frame(richness=rowSums(B40.no[,-1]>0),design=paste0("GSMC_40B",B40.no$site))
nrichness<-rbind(nrichness,nrichness1,nrichness2)
###GSMc40AB artefact
A40<-normal[,grepl(GSMc40A,names(normal))]
A40<-as.data.frame(t(A40))
A40$site<-substr(row.names(A40),1,2)
A40<-A40[,names(A40) %in% c(a.r,names(A40)[ncol(A40)])]
A40.no<-A40[,!names(A40) %in% c(a.r2)]
A40.no<-aggregate(A40.no[,1:(ncol(A40.no)-1)],by=list(site=A40.no$site),sum)
arichness1<-data.frame(richness=rowSums(A40.no[,-1]>0),design=paste0("GSMC_40A",A40.no$site))

B40<-normal[,grepl(GSMc40B,names(normal))]
B40<-as.data.frame(t(B40))
B40$site<-substr(row.names(B40),1,2)
B40<-B40[,names(B40) %in% c(a.r,names(B40)[ncol(B40)])]
B40.no<-B40[,!names(B40) %in% c(a.r2)]
B40.no<-aggregate(B40.no[,1:(ncol(B40.no)-1)],by=list(site=B40.no$site),sum)
arichness2<-data.frame(richness=rowSums(B40.no[,-1]>0),design=paste0("GSMC_40B",B40.no$site))
arichness<-rbind(arichness,arichness1,arichness2)
####deep_succ
deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ.no<-deep_succ[,names(deep_succ) %in% c(a.r,names(deep_succ)[ncol(deep_succ)])]
deep_succ.no<-aggregate(deep_succ.no[,1:(ncol(deep_succ.no)-1)],by=list(site=deep_succ.no$site),sum)
nrichness1<-data.frame(richness=rowSums(deep_succ.no[,-1]>0),design=paste0("deep_SUCC",deep_succ.no$site))

deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ<-deep_succ[,names(deep_succ) %in% c(a.r,names(deep_succ)[ncol(deep_succ)])]
deep_succ.no<-deep_succ[,!names(deep_succ) %in% c(a.r2)]
deep_succ.no<-aggregate(deep_succ.no[,1:(ncol(deep_succ.no)-1)],by=list(site=deep_succ.no$site),sum)
arichness1<-data.frame(richness=rowSums(deep_succ.no[,-1]>0),design=paste0("deep_SUCC",deep_succ.no$site))
arichness<-rbind(arichness,arichness1)
nrichness<-rbind(nrichness,nrichness1)
###proportion of artifect
p.a<-merge(arichness,nrichness,by="design")
p.a$ratio<-(p.a$richness.x/p.a$richness.y)
###pool
pool<-read.table("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/17.pool similarity/fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
row.names(pool)<-pool$OTU
pool<-pool[,-1]
n.r<-pool[row.names(pool) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata$site<-paste0(metadata$site,"DNA")
metadata$type<-"DNA"
metadata2<-metadata
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata$site<-paste0(metadata$site,"soil")
metadata$type<-"soil"
metadata<-rbind(metadata,metadata2)
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(1095:1097)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-2)],by=list(design=nr$method,site=nr$site,type=nr$type),sum)
nr <- nr[, c(rep(TRUE,3),(colSums(nr[,-c(1:3)]) != 0))]
###get the artifect ratio
a.r<-names(nr[,c(rep(TRUE,3),(colSums(nr[,-c(1:3)])/sum(nr[,-c(1:3)]))*100 <= 0.05)])
nr.no<-nr[,names(nr) %in% c(a.r)]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:3)]>0),design=paste0(nr.no$design,nr.no$site,nr$type))
###
artefact<-(nr.no[,c(rep(TRUE,3),!names(nr.no[,4:ncol(nr.no)]) %in% a.r2)])
arichness<-data.frame(richness=rowSums(artefact[,-c(1:3)]>0),design=paste0(artefact$design,artefact$site,nr$type))
###proportion of artifect
pool.p.a<-merge(arichness,nrichness,by="design",all=T)
pool.p.a$ratio<-1-(pool.p.a$richness.x/pool.p.a$richness.y)
###
pool.p.a$type<-"pooling"
p.a$type<-"unpooled"
a.r<-rbind(pool.p.a,p.a)
a.r2 <- a.r %>%
  mutate(
    site = str_extract(design, "LV|LW|LZ"), 
    pooling = str_extract(design, "DNA|soil")   ,
    design = str_remove_all(design, "LV|LW|LZ|DNA|soil")
  )
a.r2$pooling[is.na(a.r2$pooling)]<-"unpooled"
###
a.r2<-a.r2[!a.r2$design %in% c("Maestre","MDB15","MDB5"),]
a.r2$design[a.r2$design=="Deep"]<-"deep"
a.r2$design[a.r2$design=="GSMc"]<-"GSMc_62"
mod01<-lmer(ratio~ pooling + (1|site) + (1|design),data = a.r2)#[grepl("GSMc",a.r2$design),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###unique OTU number and unique artefact
######
a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
b<-read.csv("otu.select.txt",header=F)
a<-a[!duplicated(a$V1),]
a1<-a[a$V1 %in% b$V1,]
a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
#get the unique OTUs
normal<-read.csv("fungi.rarefy.table2.csv")
n.r<-normal[row.names(normal) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[,!grepl("DNA|soil",names(n.r))]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[rowSums(n.r)>0,]
n.r<-n.r[,!grepl("B",names(n.r))]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
n.r<-as.data.frame(t(n.r))
n.r$site<-row.names(n.r)
nr<-merge(n.r,metadata,by="site")
nr<-nr[,-c(350:353)]
nr$site<-substr(nr$site,1,2)
nr<-aggregate(nr[,2:(ncol(nr)-1)],by=list(design=nr$method,site=nr$site),sum)
normal<-nr
row.names(normal)<-paste0(normal$site,normal$design)
normal<-normal[,-c(1:2)]
pool<-read.table("fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
row.names(pool)<-pool$OTU
pool<-pool[,-1]
pool<-pool[row.names(pool) %in% a2$V1,]
soil<-pool[,grepl("soil",names(pool))]
DNA<-pool[,grepl("DNA",names(pool))]
names(soil)<-gsub("soil","",names(soil))
names(DNA)<-gsub("DNA","",names(DNA))
metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata$method[metadata$method=="GSMc_62"]<-"GSMc"
metadata$method[metadata$method=="deep"]<-"Deep"
names(DNA)<-paste0(substr(metadata[match(names(DNA),metadata$site),]$site,1,2),metadata[match(names(DNA),metadata$site),]$method)
names(soil)<-paste0(substr(metadata[match(names(soil),metadata$site),]$site,1,2),metadata[match(names(soil),metadata$site),]$method)
DNA<-as.data.frame(t(DNA))
soil<-as.data.frame(t(soil))

unique.richness.DNA<-data.frame()
unique.richness.soil<-data.frame()
unique.richness.s.n.u<-data.frame()
unique.richness.d.n.u<-data.frame()
for (i in row.names(normal)) {
  d.u<-names(DNA[,DNA[i,]>0])
  if (i %in% rownames(soil)) {
    s.u <- colnames(soil)[which(soil[i, ] > 0)]
  } else {
    s.u <- character(0) 
  }
  
  n.u<-names(normal[,normal[i,]>0])
  DNA.u<-setdiff(d.u,n.u)
  soil.u<-setdiff(s.u,n.u)
  s.n.u<-setdiff(n.u,s.u)
  d.n.u<-setdiff(n.u,d.u)
  unique.richness.DNA1<-data.frame(richness=length(DNA.u),method=i,site=substr(i,1,2))
  unique.richness.DNA<-rbind(unique.richness.DNA,unique.richness.DNA1)
  unique.richness.soil1<-data.frame(richness=length(soil.u),method=i,site=substr(i,1,2))
  unique.richness.soil<-rbind(unique.richness.soil,unique.richness.soil1)
  unique.richness.s.n.u1<-data.frame(richness=length(s.n.u),method=i,site=substr(i,1,2))
  unique.richness.s.n.u<-rbind(unique.richness.s.n.u,unique.richness.s.n.u1)
  unique.richness.d.n.u1<-data.frame(richness=length(d.n.u),method=i,site=substr(i,1,2))
  unique.richness.d.n.u<-rbind(unique.richness.d.n.u,unique.richness.d.n.u1)
}
###GSMc40AB
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(c(select[select$X1Amix=="GSMc40A",]$X,gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
                  gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X)),collapse = "|")
GSMc40B<-paste0(c(select[select$X1Bmix=="GSMc40B",]$X,gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
                  gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X)),collapse = "|")
A40<-n.r[grepl(GSMc40A,row.names(n.r)),]
A40$site<-substr(row.names(A40),1,2)
A40<-aggregate(A40[,1:(ncol(A40)-1)],by=list(site=A40$site),sum)
row.names(A40)<-paste0(A40$site,"GSMc_40A")
B40<-n.r[grepl(GSMc40B,row.names(n.r)),]
B40$site<-substr(row.names(B40),1,2)
B40<-aggregate(B40[,1:(ncol(B40)-1)],by=list(site=B40$site),sum)
row.names(B40)<-paste0(B40$site,"GSMc_40B")
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
deep_succ<-n.r[row.names(n.r) %in% metadata[grepl("Deep|SUCC",metadata$method),]$site,]
deep_succ$site<-substr(row.names(deep_succ),1,2)
deep_succ<-aggregate(deep_succ[,1:(ncol(deep_succ)-1)],by=list(site=deep_succ$site),sum)
row.names(deep_succ)<-paste0(deep_succ$site,"deep_SUCC")
another<-rbind(deep_succ,A40,B40)
another<-another[,-1]
another.unique.richness.DNA<-data.frame()
another.unique.richness.soil<-data.frame()
another.unique.richness.s.n.u<-data.frame()
another.unique.richness.d.n.u<-data.frame()
for (i in row.names(another)) {
  d.u<-names(DNA[,DNA[i,]>0])
  if (i %in% rownames(soil)) {
    s.u <- colnames(soil)[which(soil[i, ] > 0)]
  } else {
    s.u <- character(0) 
  }
  
  n.u<-names(another[,another[i,]>0])
  DNA.u<-setdiff(d.u,n.u)
  soil.u<-setdiff(s.u,n.u)
  s.n.u<-setdiff(n.u,s.u)
  d.n.u<-setdiff(n.u,d.u)
  another.unique.richness.DNA1<-data.frame(richness=length(DNA.u),method=i,site=substr(i,1,2))
  another.unique.richness.DNA<-rbind(another.unique.richness.DNA,another.unique.richness.DNA1)
  another.unique.richness.soil1<-data.frame(richness=length(soil.u),method=i,site=substr(i,1,2))
  another.unique.richness.soil<-rbind(another.unique.richness.soil,another.unique.richness.soil1)
  another.unique.richness.s.n.u1<-data.frame(richness=length(s.n.u),method=i,site=substr(i,1,2))
  another.unique.richness.s.n.u<-rbind(another.unique.richness.s.n.u,another.unique.richness.s.n.u1)
  another.unique.richness.d.n.u1<-data.frame(richness=length(d.n.u),method=i,site=substr(i,1,2))
  another.unique.richness.d.n.u<-rbind(another.unique.richness.d.n.u,another.unique.richness.d.n.u1)
}
###get the artifect
ar<-read.table("artefact.txt")
ar1<-read.table("thelephoraceae.artfect.txt")
ar2<-read.table("sebacinaceae.artfect.txt")
ar<-rbind(ar,ar1,ar2)
ar<-data.frame(V1=gsub("_.*","",ar$V1))
a.r2<-ar$V1
arte.unique.richness.DNA<-data.frame()
arte.unique.richness.soil<-data.frame()
arte.unique.richness.s.n.u<-data.frame()
arte.unique.richness.d.n.u<-data.frame()
for (i in row.names(normal)) {
  d.u<-names(DNA[,DNA[i,]>0])
  if (i %in% rownames(soil)) {
    s.u <- colnames(soil)[which(soil[i, ] > 0)]
  } else {
    s.u <- character(0) 
  }
  n.u<-names(normal[,normal[i,]>0])
  DNA.u1<-setdiff(d.u,n.u)
  soil.u1<-setdiff(s.u,n.u)
  s.n.u1<-setdiff(n.u,s.u)
  d.n.u1<-setdiff(n.u,d.u)
  
  DNA.u<-setdiff(DNA.u1,a.r2)
  soil.u<-setdiff(soil.u1,a.r2)
  s.n.u<-setdiff(s.n.u1,a.r2)
  d.n.u<-setdiff(d.n.u1,a.r2)
  arte.unique.richness.DNA1<-data.frame(richness=length(DNA.u),method=i,site=substr(i,1,2))
  arte.unique.richness.DNA<-rbind(arte.unique.richness.DNA,arte.unique.richness.DNA1)
  arte.unique.richness.soil1<-data.frame(richness=length(soil.u),method=i,site=substr(i,1,2))
  arte.unique.richness.soil<-rbind(arte.unique.richness.soil,arte.unique.richness.soil1)
  arte.unique.richness.s.n.u1<-data.frame(richness=length(s.n.u),method=i,site=substr(i,1,2))
  arte.unique.richness.s.n.u<-rbind(arte.unique.richness.s.n.u,arte.unique.richness.s.n.u1)
  arte.unique.richness.d.n.u1<-data.frame(richness=length(d.n.u),method=i,site=substr(i,1,2))
  arte.unique.richness.d.n.u<-rbind(arte.unique.richness.d.n.u,arte.unique.richness.d.n.u1)
}
###another artefact
arte.another.unique.richness.DNA<-data.frame()
arte.another.unique.richness.soil<-data.frame()
arte.another.unique.richness.s.n.u<-data.frame()
arte.another.unique.richness.d.n.u<-data.frame()
for (i in row.names(another)) {
  d.u<-names(DNA[,DNA[i,]>0])
  if (i %in% rownames(soil)) {
    s.u <- colnames(soil)[which(soil[i, ] > 0)]
  } else {
    s.u <- character(0) 
  }
  n.u<-names(another[,another[i,]>0])
  
  DNA.u1<-setdiff(d.u,n.u)
  soil.u1<-setdiff(s.u,n.u)
  s.n.u1<-setdiff(n.u,s.u)
  d.n.u1<-setdiff(n.u,d.u)
  DNA.u<-setdiff(DNA.u1,a.r2)
  soil.u<-setdiff(soil.u1,a.r2)
  s.n.u<-setdiff(s.n.u1,a.r2)
  d.n.u<-setdiff(d.n.u1,a.r2)
  arte.another.unique.richness.DNA1<-data.frame(richness=length(DNA.u),method=i,site=substr(i,1,2))
  arte.another.unique.richness.DNA<-rbind(arte.another.unique.richness.DNA,arte.another.unique.richness.DNA1)
  arte.another.unique.richness.soil1<-data.frame(richness=length(soil.u),method=i,site=substr(i,1,2))
  arte.another.unique.richness.soil<-rbind(arte.another.unique.richness.soil,arte.another.unique.richness.soil1)
  arte.another.unique.richness.s.n.u1<-data.frame(richness=length(s.n.u),method=i,site=substr(i,1,2))
  arte.another.unique.richness.s.n.u<-rbind(arte.another.unique.richness.s.n.u,arte.another.unique.richness.s.n.u1)
  arte.another.unique.richness.d.n.u1<-data.frame(richness=length(d.n.u),method=i,site=substr(i,1,2))
  arte.another.unique.richness.d.n.u<-rbind(arte.another.unique.richness.d.n.u,arte.another.unique.richness.d.n.u1)
}
###
DNA<-merge(arte.unique.richness.DNA[,1:2],unique.richness.DNA,by="method")
soil<-merge(arte.unique.richness.soil[,1:2],unique.richness.soil,by="method")
DNA.another<-merge(arte.another.unique.richness.DNA[,1:2],another.unique.richness.DNA,by="method")
soil.another<-merge(arte.another.unique.richness.soil[,1:2],another.unique.richness.soil,by="method")
soil<-rbind(soil.another,soil)
DNA<-rbind(DNA.another,DNA)


normal.DNA<-merge(arte.unique.richness.d.n.u[,1:2],unique.richness.d.n.u,by="method")
normal.soil<-merge(arte.unique.richness.s.n.u[,1:2],unique.richness.s.n.u,by="method")
normal.DNA.another<-merge(arte.another.unique.richness.d.n.u[,1:2],another.unique.richness.d.n.u,by="method")
normal.soil.another<-merge(arte.another.unique.richness.s.n.u[,1:2],another.unique.richness.s.n.u,by="method")
normal.soil<-rbind(normal.soil.another,normal.soil)
normal.DNA<-rbind(normal.DNA.another,normal.DNA)

###get the artifect ratio
normal.soil$ratio<-normal.soil$richness.x/normal.soil$richness.y
normal.DNA$ratio<-normal.DNA$richness.x/normal.DNA$richness.y
DNA$ratio<-DNA$richness.x/DNA$richness.y
soil$ratio<-soil$richness.x/soil$richness.y
###
soil$pooling<-"soil"
DNA$pooling<-"DNA"
normal.soil$pooling<-"unpooled"
normal.DNA$pooling<-"unpooled"
###
a.r2<-rbind(DNA,normal.DNA)
a.r3<-rbind(soil,normal.soil)
###
library(marginaleffects)
library(lme4)
library(lmerTest)
library(dplyr)
##soil
a.r3<-a.r3[!is.na(a.r3$ratio),]
a.r3$method<-gsub("LV|LZ|LW","",a.r3$method)
mod01<-lmer(richness.x~ pooling + (1|site) + (1|method),data = a.r3)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
mod01<-lmer(richness.y~ pooling + (1|site) + (1|method),data = a.r3)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 

###DNA
a.r2<-a.r2[!is.na(a.r2$ratio),]
a.r2$method<-gsub("LV|LZ|LW","",a.r2$method)
mod01<-lmer(richness.x~ pooling + (1|site) + (1|method),data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
mod01<-lmer(richness.y~ pooling + (1|site) + (1|method),data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###
mod01<-lmer(richness.y~ pooling + (1|site) + method,data = a.r3)#[grepl("GSMc",a.r2$design),])
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
plot_predictions(mod01,condition = c("method","pooling") )
##group into >20
##soil
a.r3$sub<-NA
a.r3$sub[grepl("GSMc|GSMC|Zobel|SUCC",a.r3$method)]<- ">20"
a.r3$sub[!grepl("GSMc|GSMC|Zobel|SUCC",a.r3$method)]<- "<20"
mod01<-lmer(richness.y~ (1|site)+ (1|method)+ pooling ,data = a.r3[grepl("<20",a.r3$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
mod01<-lmer(richness.y~ (1|site)+ (1|method)+pooling ,data = a.r3[grepl(">20",a.r3$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
##group into >20
a.r2$sub<-NA
a.r2$sub[grepl("GSMc|GSMC",a.r2$method)]<- ">20"
a.r2$sub[!grepl("GSMc|GSMC",a.r2$method)]<- "<20"
mod01<-lmer(richness.y~ (1|site)+ (1|method)+ pooling ,data = a.r2[grepl("<20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
mod01<-lmer(richness.y~ (1|site)+ (1|method)+pooling ,data = a.r2[grepl(">20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 








