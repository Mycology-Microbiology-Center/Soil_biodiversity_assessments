# setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
# a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
# b<-read.csv("otu.select.txt",header=F)
# a<-a[!duplicated(a$V1),]
# a1<-a[a$V1 %in% b$V1,]
# a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
# ####
# normal<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/diversity original_not summarize/fungi.rarefy.table2.csv")
# n.r<-normal[row.names(normal) %in% a2$V1,]
# n.r<-n.r[,!grepl("DNA|soil",names(n.r))]
# n.r<-n.r[,!grepl("B",names(n.r))]
# n.r<-n.r[,colSums(n.r)>0]
# n.r<-n.r[rowSums(n.r)>0,]
# a.r<-row.names(n.r[(rowSums(n.r)/sum(n.r))*100 <= 0.05,])
# a2<-a[a$V1 %in% a.r,]
# try<-a2
# ###select for pool
# setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
# a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
# b<-read.csv("otu.select.txt",header=F)
# a<-a[!duplicated(a$V1),]
# a1<-a[a$V1 %in% b$V1,]
# a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
# pool<-read.table("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/17.pool similarity/fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
# row.names(pool)<-pool$OTU
# pool<-pool[,-1]
# n.r<-pool[row.names(pool) %in% a2$V1,]
# n.r<-n.r[,colSums(n.r)>0]
# n.r<-n.r[rowSums(n.r)>0,]
# metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
# metadata$site<-paste0(metadata$site,"DNA")
# metadata$type<-"DNA"
# metadata2<-metadata
# metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
# metadata$site<-paste0(metadata$site,"soil")
# metadata$type<-"soil"
# metadata<-rbind(metadata,metadata2)
# n.r<-as.data.frame(t(n.r))
# n.r$site<-row.names(n.r)
# nr<-merge(n.r,metadata,by="site")
# nr<-nr[,-c(501:503)]
# nr$site<-substr(nr$site,1,2)
# nr<-aggregate(nr[,2:(ncol(nr)-2)],by=list(design=nr$method,site=nr$site,type=nr$type),sum)
# nr <- nr[, c(rep(TRUE,3),(colSums(nr[,-c(1:3)]) != 0))]
# ###get the artifect ratio
# a.r<-names(nr[,(colSums(nr[,-c(1:3)])/sum(nr[,-c(1:3)]))*100<=0.05])
# a2<-a[a$V1 %in% a.r,]
# try2<-a2
# all<-rbind(try,try2)
# all<-all[!duplicated(all$V1),]
# a2<-all
# ###
# f__Russula<-a2[grepl("f__Russula",a2$V2),]
# f__Thelephoraceae<-a2[grepl("f__Thelephoraceae",a2$V2),]
# f__Sebacinaceae<-a2[grepl("f__Sebacinaceae",a2$V2),]
# setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy/correct1")
# write.csv(f__Russula[,1:2],"Russula.artefact2.csv")
# write.csv(f__Thelephoraceae[,1:2],"Thelephoraceae.artefact2.csv")
# write.csv(f__Sebacinaceae[,1:2],"Sebacinaceae.artefact2.csv")
####calculate the proportion
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
b<-read.csv("otu.select.txt",header=F)
a<-a[!duplicated(a$V1),]
a1<-a[a$V1 %in% b$V1,]
a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
###get rare OTU ratio
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
a.r<-names(nr[,c(rep(TRUE,2),(colSums(nr[,-c(1:2)])/sum(nr[,-c(1:2)]))*100 <= 0.05)])
nr.no<-nr[,names(nr) %in% c(a.r)]
nrichness<-data.frame(richness=rowSums(nr.no[,-c(1:2)]>0),design=paste0(nr.no$design,nr.no$site))
##
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
ar<-read.csv("nonartefact.csv",row.names = 1)
a.r2<-ar$x
artefact<-(nr.no[,c(rep(TRUE,2),!names(nr.no[,3:ncol(nr.no)]) %in% a.r2)])
arichness<-data.frame(richness=rowSums(artefact[,-c(1:2)]>0),design=paste0(artefact$design,artefact$site))
###GSMc40AB
select<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/1.accumulative/sheet2_for_pooled.csv")
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
write.csv(pool.p.a,"high.rare.artefact.three.families.pool.p.a.csv")
write.csv(p.a,"rare.artefact.three.families.a.csv")
###
pool.p.a$type<-"pooling"
p.a$type<-"unpooled"
a.r<-rbind(pool.p.a,p.a)
library(dplyr)
library(stringr)
a.r2 <- a.r %>%
  mutate(
    site = str_extract(design, "LV|LW|LZ"), 
    pooling = str_extract(design, "DNA|soil")   ,
    design = str_remove_all(design, "LV|LW|LZ|DNA|soil")
  )
a.r2$pooling[is.na(a.r2$pooling)]<-"unpooled"
###
library(marginaleffects)
library(lme4)
library(lmerTest)
library(dplyr)
a.r2<-a.r2[!a.r2$design %in% c("Maestre","MDB15","MDB5"),]
a.r2$design[a.r2$design=="Deep"]<-"deep"
a.r2$design[a.r2$design=="GSMc"]<-"GSMc_62"
a.r2$design<-gsub("GSMC","GSMc",a.r2$design)
mod01<-lmer(ratio~ pooling + (1|site) + design,data = a.r2)#[grepl("GSMc",a.r2$design),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise"))
library(effectsize)
eta_squared(mod01)
performance::performance(mod01)
a.r2$rare.try<-a.r2$richness.y-a.r2$richness.x
mod01<-lmer((richness.y-richness.x)~ pooling + (1|site) + design,data = a.r2)#[grepl("GSMc",a.r2$design),])
eta_squared(mod01)
performance::performance(mod01)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 

####
# a.r2$sub<-NA
# a.r2$sub[grepl("GSMc|GSMC",a.r2$design)] <- ">20"
# a.r2$sub[!grepl("GSMc|GSMC",a.r2$design)] <- "<20"
# a.r2$sub[a.r2$pooling=="soil" & a.r2$design=="Zobel"] <- ">20"
# a.r2$sub[a.r2$pooling=="soil" & a.r2$design=="SUCC"] <- ">20"
# mod01<-lmer(ratio~ (1|site)+ (1|design)+ pooling ,data = a.r2[grepl("<20",a.r2$sub),])
# b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
# 
# mod01<-lmer(ratio~ (1|site)+ (1|design)+ pooling ,data = a.r2[grepl(">20",a.r2$sub),])
# b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 

###get the ratio
pool.p.a$design2<-pool.p.a$design
pool.p.a$design<-gsub("DNA|soil","",pool.p.a$design)
pool.p.a$design<-gsub("deep","Deep",pool.p.a$design)
p.a$design<-gsub("deep","Deep",p.a$design)
pool.p.a$design<-gsub("GSMc_62","GSMc",pool.p.a$design)
p.a$design<-gsub("GSMC","GSMc",p.a$design)
pool.p.a$design2[grepl("soil",pool.p.a$design2)]<-"soil"
pool.p.a$design2[grepl("DNA",pool.p.a$design2)]<-"DNA"
ratio<-merge(pool.p.a[,-5], p.a[,1:3], by="design",all=T)
ratio$try<-(ratio$richness.x.y/ratio$richness.x.x)
ratio$sub<-NA
ratio$sub[grepl("GSMc|GSMC",ratio$design)] <- ">20"
ratio$sub[!grepl("GSMc|GSMC",ratio$design)] <- "<20"
ratio$sub[ratio$design2=="soil" & grepl("Zobel|SUCC",ratio$design)] <- ">20"
ratio$try<-(ratio$richness.x.y/ratio$richness.x.x)
ratio$site<-substr(ratio$design,nchar(ratio$design)-1,nchar(ratio$design))
get <- ratio %>%
  group_by(site,sub,design2) %>%
  mutate(try2=median(try),mad=mad(try)) %>%
  ungroup()
###
ratio<-merge(pool.p.a[,-5], p.a[,1:4], by="design",all=T)
ratio$try<-(ratio$ratio.x/ratio$ratio.y)
ratio$sub<-NA
ratio$sub[grepl("GSMc|GSMC",ratio$design)] <- ">20"
ratio$sub[!grepl("GSMc|GSMC",ratio$design)] <- "<20"
ratio$sub[ratio$design2=="soil" & grepl("Zobel|SUCC",ratio$design)] <- ">20"
ratio$site<-substr(ratio$design,nchar(ratio$design)-1,nchar(ratio$design))
get <- ratio %>%
  group_by(site,design2) %>%
  mutate(try2=median(try),mad=mad(try)) %>%
  ungroup()
