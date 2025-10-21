######
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
###get the artifect ratio
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
ar<-read.csv("nonartefact.csv",row.names = 1)
ar<-ar$x
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
pool<-read.table("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/17.pool similarity/fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
row.names(pool)<-pool$OTU
pool<-pool[,-1]
n.r<-pool[row.names(pool) %in% a2$V1,]
n.r<-n.r[,colSums(n.r)>0]
n.r<-n.r[rowSums(n.r)>0,]
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
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy/correct1")
write.csv(pool.p.a,"three.families.pool.p.a.csv")
write.csv(p.a,"three.families.a.csv")
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
mod01<-lmer(ratio~ pooling + (1|site) + design,data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
##group into >20
a.r2$sub<-NA
a.r2$sub[grepl("GSMc|GSMC",a.r2$method)]<- ">20"
a.r2$sub[!grepl("GSMc|GSMC",a.r2$method)]<- "<20"
mod01<-lmer(ratio~ (1|site)+ design+ pooling ,data = a.r2[!grepl("<20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise"))
b
mod01<-lmer(ratio~ (1|site)+ design+pooling ,data = a.r2[!grepl(">20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
b
#######get the exact ratio
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

