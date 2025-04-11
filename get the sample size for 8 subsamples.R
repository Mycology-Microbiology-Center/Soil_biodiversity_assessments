###check the increasing pattern between the richness and sample size of different organisms
library(iNEXT.3D)
library(metagMisc)
###animal
normal<-read.table("animal.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata.normal2<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
normal.fre2<-merge(normal.fre,metadata.normal2,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$COI_Otu1),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method!="Meastre",]
normal.fre2<-normal.fre2[,1:(ncol(normal.fre2)-2)]
LV<-normal.fre2[grepl("LV",normal.fre2$site),]
LV<-LV[,-1]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
LW<-normal.fre2[grepl("LW",normal.fre2$site),]
LW<-LW[,-1]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-normal.fre2[grepl("LZ",normal.fre2$site),]
LZ<-LZ[,-1]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
###
LZ<-as.data.frame(ifelse(LZ>0,1,0))
LV<-as.data.frame(ifelse(LV>0,1,0))
LW<-as.data.frame(ifelse(LW>0,1,0))
animal<-list(LZ=t(LZ),LV=t(LV),LW=t(LW))
animal.out <- iNEXT3D(animal, diversity = 'TD', q = 0, datatype="incidence_raw",endpoint = 500)
ggiNEXT3D(animal.out, type = 1, facet.var = "Assemblage")
#fungi
normal<-read.table("fungi.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata.normal2<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
normal.fre2<-merge(normal.fre,metadata.normal2,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$'0003b6ab41bde80a3cc58f45791695956d51cb52'),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method!="Meastre",]
normal.fre2<-normal.fre2[,1:(ncol(normal.fre2)-2)]
LV<-normal.fre2[grepl("LV",normal.fre2$site),]
LV<-LV[,-1]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
LW<-normal.fre2[grepl("LW",normal.fre2$site),]
LW<-LW[,-1]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-normal.fre2[grepl("LZ",normal.fre2$site),]
LZ<-LZ[,-1]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
###
LZ<-as.data.frame(ifelse(LZ>0,1,0))
LV<-as.data.frame(ifelse(LV>0,1,0))
LW<-as.data.frame(ifelse(LW>0,1,0))
fungi<-list(LZ=t(LZ),LV=t(LV),LW=t(LW))
fungi.out <- iNEXT3D(fungi, diversity = 'TD', q = 0, datatype="incidence_raw",endpoint = 500)
ggiNEXT3D(fungi.out, type = 1, facet.var = "Assemblage")
#bacteria
normal<-read.table("bacteria.rarefy.table.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata.normal2<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
normal.fre2<-merge(normal.fre,metadata.normal2,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$'0000294ab66fb96003e4afd4237026e0'),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method!="Meastre",]
normal.fre2<-normal.fre2[,1:(ncol(normal.fre2)-2)]
LV<-normal.fre2[grepl("LV",normal.fre2$site),]
LV<-LV[,-1]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
LW<-normal.fre2[grepl("LW",normal.fre2$site),]
LW<-LW[,-1]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-normal.fre2[grepl("LZ",normal.fre2$site),]
LZ<-LZ[,-1]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
###
LZ<-as.data.frame(ifelse(LZ>0,1,0))
LV<-as.data.frame(ifelse(LV>0,1,0))
LW<-as.data.frame(ifelse(LW>0,1,0))
bacteria<-list(LZ=t(LZ),LV=t(LV),LW=t(LW))
bacteria.out <- iNEXT3D(bacteria, diversity = 'TD', q = 0, datatype="incidence_raw",endpoint = 500)
ggiNEXT3D(bacteria.out, type = 1, facet.var = "Assemblage")

####assess the best sampling area for 8 subsamples(4 pairs)
animal<- read.table("animal.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
animal2<-animal[,!grepl("DNA|soil",names(animal))]
animal2<-as.data.frame(t(animal2))
metadata<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
animal2<-animal2[grepl(paste0(metadata[grepl("GSMc",metadata$method),]$site,collapse = "|"),row.names(animal2)),]
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

###
geta.d<-function(site.table,coordinate){
  threshold=2
  select<-unique(gsub("a|b","",row.names(site.table))[duplicated(gsub("a|b","",row.names(site.table)))])
  combn8<-combn(select,4,simplify = F)
  to_remove <- rep(FALSE, length(combn8))  
  for (i in seq_along(combn8)) {
    if (to_remove[i]) next  
    
    for (j in seq_along(combn8)) {
      if (i != j && !to_remove[j]) {  
        if (length(intersect(combn8[[i]], combn8[[j]])) >= threshold) {
          to_remove[j] <- TRUE  
        }
      }
    }
  }
  combn8<-(combn8[!to_remove])
  names(combn8)<-paste0("A",seq_along(combn8))
  rsquare<-data.frame()
  for (i in seq_along(combn8)) {
    table<-(coordinate[grepl(paste0(combn8[[i]],collapse = "|"),coordinate$sample),])
    rsquare1<-data.frame(x=max(table$x)-min(table$x),y=max(table$y)-min(table$y),sample=names(combn8)[i])
    rsquare<-rbind(rsquare,rsquare1)
  }
  ##calculate the diversity
  site.table<-as.data.frame(t(site.table))
  diversity<-data.frame()
  for (i in seq_along(combn8)) {
    table<-(site.table[,grepl(paste0(combn8[[i]],collapse = "|"),names(site.table))])
    diversity1<-data.frame(richness=sum((rowSums(table))>0),sample=names(combn8)[i])
    diversity<-rbind(diversity1,diversity)
  }
  sample8<-merge(diversity,rsquare,by="sample")
  sample8$area<-(((sample8$y*sample8$y)+(sample8$x*sample8$x))/4)*pi
  sample8<-sample8[sample8$area<=2000,]
  return(sample8)
}
animal.LW.diversity<-geta.d(LW,coordinate.LW)
animal.LV.diversity<-geta.d(LV,coordinate.LV)
animal.LZ.diversity<-geta.d(LZ,coordinate.LZ)
animal.LZ.diversity$site<-"LZ"
animal.LV.diversity$site<-"LV"
animal.LW.diversity$site<-"LW"
##bacteria
bacteria<- read.table("bacteria.rarefy.table.csv",header = T,row.names = 1,sep = ",")
bacteria2<-bacteria[,!grepl("DNA|soil",names(bacteria))]
bacteria2<-as.data.frame(t(bacteria2))
metadata<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
bacteria2<-bacteria2[grepl(paste0(metadata[grepl("GSMc",metadata$method),]$site,collapse = "|"),row.names(bacteria2)),]
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
bacteria.LW.diversity<-geta.d(LW,coordinate.LW)
bacteria.LV.diversity<-geta.d(LV,coordinate.LV)
bacteria.LZ.diversity<-geta.d(LZ,coordinate.LZ)
bacteria.LZ.diversity$site<-"LZ"
bacteria.LV.diversity$site<-"LV"
bacteria.LW.diversity$site<-"LW"
#fungi
fungi<- read.table("fungi.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
fungi2<-fungi[,!grepl("DNA|soil",names(fungi))]
fungi2<-as.data.frame(t(fungi2))
metadata<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
fungi2<-fungi2[grepl(paste0(metadata[grepl("GSMc",metadata$method),]$site,collapse = "|"),row.names(fungi2)),]
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
fungi.LW.diversity<-geta.d(LW,coordinate.LW)
fungi.LV.diversity<-geta.d(LV,coordinate.LV)
fungi.LZ.diversity<-geta.d(LZ,coordinate.LZ)
fungi.LZ.diversity$site<-"LZ"
fungi.LV.diversity$site<-"LV"
fungi.LW.diversity$site<-"LW"
###
all<-rbind(fungi.LW.diversity,fungi.LZ.diversity,fungi.LV.diversity)
all2<-all[abs(all$y-all$x)<=5,]
mod<-lmer(richness~area+(1|site),data = all2)
summary(mod)
all<-rbind(bacteria.LW.diversity,bacteria.LZ.diversity,bacteria.LV.diversity)
all2<-all[abs(all$y-all$x)<=5,]
mod<-lmer(richness~area+(1|site),data = all2)
summary(mod)
all<-rbind(animal.LW.diversity,animal.LZ.diversity,animal.LV.diversity)
all2<-all[abs(all$y-all$x)<=5,]
mod<-lmer(richness~area+(1|site),data = all2)
summary(mod)
##get the minimum area
all<-rbind(fungi.LW.diversity,fungi.LZ.diversity,fungi.LV.diversity)
all2<-all[abs(all$y-all$x)<=5,]
all2<-all2[all2$y>19&all2$x>19,]
mod<-lmer(richness~area+(1|site),data = all2)
summary(mod)
all<-rbind(bacteria.LW.diversity,bacteria.LZ.diversity,bacteria.LV.diversity)
all2<-all[abs(all$y-all$x)<=5,]
all2<-all2[all2$y>6.2&all2$x>6.2,]
mod<-lmer(richness~area+(1|site),data = all2)
summary(mod)
all<-rbind(animal.LW.diversity,animal.LZ.diversity,animal.LV.diversity)
all2<-all[abs(all$y-all$x)<=5,]
all2<-all2[all2$y>6.5&all2$x>6.5,]
mod<-lmer(richness~area+(1|site),data = all2)
summary(mod)

###calculate the coverage
###animal
normal<-read.table("animal.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata.normal2<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
normal.fre2<-merge(normal.fre,metadata.normal2,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$COI_Otu1),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method=="GSMc",]
normal.fre2<-normal.fre2[,1:(ncol(normal.fre2)-2)]
LV<-normal.fre2[grepl("LV",normal.fre2$site),]
row.names(LV)<-LV$site
LV<-LV[,-1]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
LW<-normal.fre2[grepl("LW",normal.fre2$site),]
row.names(LW)<-LW$site
LW<-LW[,-1]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-normal.fre2[grepl("LZ",normal.fre2$site),]
row.names(LZ)<-LZ$site
LZ<-LZ[,-1]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
LZ<-as.data.frame(t(LZ))
LW<-as.data.frame(t(LW))
LV<-as.data.frame(t(LV))
###
select<-unique(gsub("a|b","",names(LV))[duplicated(gsub("a|b","",names(LV)))])
LV<-LV[,grepl(paste0(sample(select,10,replace = F),collapse = "|"),names(LV))]
select<-unique(gsub("a|b","",names(LW))[duplicated(gsub("a|b","",names(LW)))])
LW<-LW[,grepl(paste0(sample(select,10,replace = F),collapse = "|"),names(LW))]
select<-unique(gsub("a|b","",names(LZ))[duplicated(gsub("a|b","",names(LZ)))])
LZ<-LZ[,grepl(paste0(sample(select,10,replace = F),collapse = "|"),names(LZ))]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
###
animal<-list(LZ=as.numeric(rowSums(LZ)),LV=as.numeric(rowSums(LV)),LW=as.numeric(rowSums(LW)))
animal.out <- iNEXT3D(animal, diversity = 'TD', q = 0, datatype="abundance")
ggiNEXT3D(animal.out, type = 3, facet.var = "Assemblage")

#fungi
normal<-read.table("fungi.rarefy.table2.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata.normal2<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
normal.fre2<-merge(normal.fre,metadata.normal2,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$'0003b6ab41bde80a3cc58f45791695956d51cb52'),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method=="GSMc",]
normal.fre2<-normal.fre2[,1:(ncol(normal.fre2)-2)]
LV<-normal.fre2[grepl("LV",normal.fre2$site),]
row.names(LV)<-LV$site
LV<-LV[,-1]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
LW<-normal.fre2[grepl("LW",normal.fre2$site),]
row.names(LW)<-LW$site
LW<-LW[,-1]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-normal.fre2[grepl("LZ",normal.fre2$site),]
row.names(LZ)<-LZ$site
LZ<-LZ[,-1]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
LZ<-as.data.frame(t(LZ))
LW<-as.data.frame(t(LW))
LV<-as.data.frame(t(LV))
###
select<-unique(gsub("a|b","",names(LV))[duplicated(gsub("a|b","",names(LV)))])
LV<-LV[,grepl(paste0(sample(select,15,replace = F),collapse = "|"),names(LV))]
select<-unique(gsub("a|b","",names(LW))[duplicated(gsub("a|b","",names(LW)))])
LW<-LW[,grepl(paste0(sample(select,15,replace = F),collapse = "|"),names(LW))]
select<-unique(gsub("a|b","",names(LZ))[duplicated(gsub("a|b","",names(LZ)))])
LZ<-LZ[,grepl(paste0(sample(select,15,replace = F),collapse = "|"),names(LZ))]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
###
fungi<-list(LZ=as.numeric(rowSums(LZ)),LV=as.numeric(rowSums(LV)),LW=as.numeric(rowSums(LW)))
fungi.out <- iNEXT3D(fungi, diversity = 'TD', q = 0, datatype="abundance")
ggiNEXT3D(fungi.out, type = 3, facet.var = "Assemblage")

#bacteria
normal<-read.table("bacteria.rarefy.table.csv",header = T,row.names = 1,sep = ",")
normal.fre<-as.data.frame(t(normal))
normal.fre$site<-row.names(normal.fre) 
normal.fre$site2<-substr(normal.fre$site,1,2) 
metadata.normal2<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata.normal2<-metadata.normal2[,c(1,2)]
normal.fre2<-merge(normal.fre,metadata.normal2,by="site",all=T)
normal.fre2<-normal.fre2[!is.na(normal.fre2$'0000294ab66fb96003e4afd4237026e0'),]
normal.fre2<-normal.fre2[!grepl("B1|B2|B3|B4|B5",normal.fre2$site),]
normal.fre2<-normal.fre2[!is.na(normal.fre2$method),]
normal.fre2<-normal.fre2[normal.fre2$method=="GSMc",]
normal.fre2<-normal.fre2[,1:(ncol(normal.fre2)-2)]
LV<-normal.fre2[grepl("LV",normal.fre2$site),]
row.names(LV)<-LV$site
LV<-LV[,-1]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
LW<-normal.fre2[grepl("LW",normal.fre2$site),]
row.names(LW)<-LW$site
LW<-LW[,-1]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-normal.fre2[grepl("LZ",normal.fre2$site),]
row.names(LZ)<-LZ$site
LZ<-LZ[,-1]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
LZ<-as.data.frame(t(LZ))
LW<-as.data.frame(t(LW))
LV<-as.data.frame(t(LV))
###
select<-unique(gsub("a|b","",names(LV))[duplicated(gsub("a|b","",names(LV)))])
LV<-LV[,grepl(paste0(sample(select,15,replace = F),collapse = "|"),names(LV))]
select<-unique(gsub("a|b","",names(LW))[duplicated(gsub("a|b","",names(LW)))])
LW<-LW[,grepl(paste0(sample(select,15,replace = F),collapse = "|"),names(LW))]
select<-unique(gsub("a|b","",names(LZ))[duplicated(gsub("a|b","",names(LZ)))])
LZ<-LZ[,grepl(paste0(sample(select,15,replace = F),collapse = "|"),names(LZ))]
LW<-LW[rowSums(LW)>0,]
LW<-LW[,colSums(LW)>0]
LZ<-LZ[rowSums(LZ)>0,]
LZ<-LZ[,colSums(LZ)>0]
LV<-LV[rowSums(LV)>0,]
LV<-LV[,colSums(LV)>0]
###
bacteria<-list(LZ=as.numeric(rowSums(LZ)),LV=as.numeric(rowSums(LV)),LW=as.numeric(rowSums(LW)))
bacteria.out <- iNEXT3D(bacteria, diversity = 'TD', q = 0, datatype="abundance")
ggiNEXT3D(bacteria.out, type = 3, facet.var = "Assemblage")







