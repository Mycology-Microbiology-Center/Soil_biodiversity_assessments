######
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
a<-read.csv("OTUs_LULU_blast_1st_tophit.txt",sep="+",header=F)
b<-read.csv("otu.select.txt",header=F)
a<-a[!duplicated(a$V1),]
a1<-a[a$V1 %in% b$V1,]
a2<-a1[grepl("f__Russula|f__Thelephoraceae|f__Sebacinaceae",a1$V2),]
#get the unique OTUs
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
normal<-nr
row.names(normal)<-paste0(normal$site,normal$design)
normal<-normal[,-c(1:2)]
pool<-read.table("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/17.pool similarity/fungi.poolat100.csv",header = T,row.names = 1,sep = ",")
row.names(pool)<-pool$OTU
pool<-pool[,-1]
pool<-pool[row.names(pool) %in% a2$V1,]
soil<-pool[,grepl("soil",names(pool))]
DNA<-pool[,grepl("DNA",names(pool))]
names(soil)<-gsub("soil","",names(soil))
names(DNA)<-gsub("DNA","",names(DNA))
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
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
select<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/1.accumulative/sheet2_for_pooled.csv")
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
metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/normal_metadata.csv",header = TRUE,sep = ",")
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
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/pooled rarefy")
ar<-read.csv("nonartefact.csv",row.names = 1)
a.r2<-ar$x
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
mod01<-lmer(richness.x~ pooling + (1|site) + method,data = a.r3)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) ##unique arterfacts
eta_squared(mod01)
performance::performance(mod01)

mod01<-lmer(richness.y~ pooling + (1|site) + method,data = a.r3)##unique OTUs
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
mod01<-lmer(ratio~ pooling + (1|site) + (1|method),data = a.r3)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###DNA
a.r2<-a.r2[!is.na(a.r2$ratio),]
a.r2$method<-gsub("LV|LZ|LW","",a.r2$method)
mod01<-lmer(richness.x~ pooling + (1|site) + method,data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
mod01<-lmer(richness.y~ pooling + (1|site) + method,data = a.r2)##unique OTUs
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)

mod01<-lmer(ratio~ pooling + (1|site) + (1|method),data = a.r2)
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
###
mod01<-lmer(richness.y~ pooling + (1|site) + method,data = a.r2)#[grepl("GSMc",a.r2$design),])
b<-avg_comparisons(mod01, variables = list(method = "pairwise")) 
plot_predictions(mod01,condition = c("method","pooling") )
##group into >20
##soil
a.r3$sub<-NA
a.r3$sub[grepl("GSMc|GSMC|Zobel|SUCC",a.r3$method)]<- ">20"
a.r3$sub[!grepl("GSMc|GSMC|Zobel|SUCC",a.r3$method)]<- "<20"
mod01<-lmer(richness.y~ (1|site)+ method+ pooling ,data = a.r3[grepl("<20",a.r3$sub),])##PE>1
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
mod01<-lmer(richness.y~ (1|site)+ method+pooling ,data = a.r3[grepl(">20",a.r3$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
##group into >20
a.r2$sub<-NA
a.r2$sub[grepl("GSMc|GSMC",a.r2$method)]<- ">20"
a.r2$sub[!grepl("GSMc|GSMC",a.r2$method)]<- "<20"
mod01<-lmer(richness.y~ (1|site)+ method+ pooling ,data = a.r2[grepl("<20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
mod01<-lmer(richness.y~ (1|site)+ method+pooling ,data = a.r2[grepl(">20",a.r2$sub),])
b<-avg_comparisons(mod01, variables = list(pooling = "pairwise")) 
eta_squared(mod01)
performance::performance(mod01)
###get the ratio
ratio<-merge(DNA, normal.DNA, by="method",all=T)
ratio$try<-(ratio$richness.y.y/ratio$richness.y.x)
ratio$sub<-NA
ratio$sub[grepl("GSMc|GSMC",ratio$method)]<- ">20"
ratio$sub[!grepl("GSMc|GSMC",ratio$method)]<- "<20"
get <- ratio %>%
  group_by(site.y,sub) %>%
  mutate(try2=median(try),mad=mad(try)) %>%
  ungroup()
###unique artefacts
ratio$try<-(ratio$richness.x.y/ratio$richness.x.x)
get <- ratio %>%
  group_by(site.y,sub) %>%
  mutate(try2=median(try),mad=mad(try)) %>%
  ungroup()
###
ratio<-merge(soil, normal.soil, by="method",all=T)
ratio$try<-(ratio$richness.y.y/ratio$richness.y.x)
ratio$sub<-NA
ratio$sub[grepl("GSMc|GSMC|Zobel|SUCC",ratio$method)]<- ">20"
ratio$sub[!grepl("GSMc|GSMC|Zobel|SUCC",ratio$method)]<- "<20"
get <- ratio %>%
  group_by(site.y,sub) %>%
  mutate(try2=median(try),mad=mad(try)) %>%
  ungroup()
###unique artefacts
ratio$try<-(ratio$richness.x.y/ratio$richness.x.x)
get <- ratio %>%
  group_by(site.y,sub) %>%
  mutate(try2=median(try),mad=mad(try)) %>%
  ungroup()