##bacteria
library(dplyr)
library(vegan)
library(metagMisc)
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
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
#DNA
DNA<-read.table("bacteria.poolat100.csv",header = T,row.names = 1,sep = ",")
DNA<-DNA[,grepl("DNA",names(DNA))]
names(DNA)<-gsub("DNA","",names(DNA))
metadata.DNA<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata.DNA[,c(1,5)]
bacteria.normal<-as.data.frame(t(DNA))
bacteria.normal3<-decostand(bacteria.normal,method = "total")
dist_normal<- vegdist(bacteria.normal3, method = "bray")
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
eta_squared(mod01)
performance::performance(mod01)
options(scipen = 999)
animal.result<-parameters::model_parameters(mod01)
animal.result$organism<-'animal'
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
plot_predictions(mod01,condition = "pooling")+
  labs(x="Pooling types",y="Similarity")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
###bacteria
bacteria<-rbind(pool.bacteria.between,bacteria.between)
mod01<-lmer(similarity ~  pooling + (1|site),data = bacteria)
performance::performance(mod01)
options(scipen = 999)
bacteria.result<-parameters::model_parameters(mod01)
bacteria.result$organism<-'bacteria'
eta_squared(mod01)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
plot_predictions(mod01,condition = "pooling")+
  labs(x="Pooling types",y="Similarity")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
###fungi
fungi<-rbind(pool.fungi.between,fungi.between)
mod01<-lmer(similarity ~  pooling + (1|site),data = fungi)
performance::performance(mod01)
options(scipen = 999)
fungi.result<-parameters::model_parameters(mod01)
fungi.result$organism<-'fungi'
eta_squared(mod01)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
plot_predictions(mod01,condition = "pooling")+
  labs(x="Pooling types",y="Similarity")+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )
result<-rbind(bacteria.result,
              fungi.result,
              animal.result)
write.csv(result,"similarity.parameters.csv")


