##animal
library(data.table)
animal.DNA<-read.csv("animal_DNA2.csv",header = TRUE,sep = ",",row.names = 1)
animal.DNA<- animal.DNA[animal.DNA$kBP>0.8,]
tax<-animal.DNA[,1:22]
animal.DNA<-animal.DNA[,c(1,23:58)]

names(animal.DNA)<-gsub(".*_","",names(animal.DNA))
names(animal.DNA)<-gsub("mix","",names(animal.DNA))
row.names(animal.DNA)<-animal.DNA$GlobalESV
animal.DNA<-animal.DNA[,c(2:37)]
animal.DNA$Feature.ID<-row.names(animal.DNA)
animal.DNA<-animal.DNA[rowSums(animal.DNA[,c(1:36)])!=0,]

animal.DNA2<-merge(animal.DNA,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.DNA2<-animal.DNA2[!is.na(animal.DNA2$LZ1A),]
animal.DNA2<-animal.DNA2[!grepl("Unassigned",animal.DNA2$Kingdom),]
names(animal.DNA2)[2:37]<-paste0(names(animal.DNA2)[2:37],"DNA")
tax_table.DNA<-animal.DNA2[,c(38:44)]
row.names(tax_table.DNA)<-animal.DNA2$Feature.ID

##soil
animal.soil<-read.csv("animal_soil2.csv",header = TRUE,sep = ",",row.names = 1)
animal.soil<- animal.soil[animal.soil$kBP>0.8,]
animal.soil<-animal.soil[,c(1,23:58)]
names(animal.soil)<-gsub(".*_","",names(animal.soil))
names(animal.soil)<-gsub("mix","",names(animal.soil))

tax<-read.table("soil_taxonomy.csv",header = T, sep=",")
row.names(animal.soil)<-animal.soil$GlobalESV
animal.soil<-animal.soil[,c(2:37)]
animal.soil$Feature.ID<-row.names(animal.soil)
animal.soil<-animal.soil[rowSums(animal.soil[,c(1:36)])!=0,]

animal.soil2<-merge(animal.soil,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.soil2<-animal.soil2[!is.na(animal.soil2$LZ1A),]
animal.soil2<-animal.soil2[!grepl("Unassigned",animal.soil2$Kingdom),]
names(animal.soil2)[2:37]<-paste0(names(animal.soil2)[2:37],"soil")
tax_table.soil<-animal.soil2[,c(38:45)]
row.names(tax_table.soil)<-animal.soil2$Feature.ID


mer<-function(x,y){
  merge(x,y,by="Feature.ID",all=T)
}

animal<-list(animal.DNA2,animal.soil2)
animal2<-Reduce(mer,animal)

row.names(animal2)<-animal2$Feature.ID
animal2<-animal2[,grepl("LZ|LV|LW",names(animal2))]

animal2[is.na(animal2)]<-0

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")
metadata<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
metadata2<-rbind(rbind(metadata.DNA,metadata.soil),metadata[,-5])

##DarkDiv
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
metadata.normal<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")
DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-animal.normal[,grepl(DarkDiv1,names(animal.normal))]
DarkDiv<-DarkDiv[,grepl("LZ",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))

d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
library(vegan)
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-animal.normal[,grepl(DarkDiv1,names(animal.normal))]
DarkDiv<-DarkDiv[,grepl("LV",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-animal.normal[,grepl(DarkDiv1,names(animal.normal))]
DarkDiv<-DarkDiv[,grepl("LW",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
DarkDiv<-shannon.normal2
DarkDiv$method<-"DarkDiv"
##Zobel
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
Zobel1<-paste0(metadata.normal$site[metadata.normal$method=="Zobel"],collapse="|")
Zobel<-animal.normal[,grepl(Zobel1,names(animal.normal))]
Zobel<-Zobel[,grepl("LZ",names(Zobel))]
Zobel<-as.data.frame(t(Zobel))

d<-data.frame(try=colSums(Zobel))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(Zobel),sequencing.depth=sum(Zobel),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
Zobel1<-paste0(metadata.normal$site[metadata.normal$method=="Zobel"],collapse="|")
Zobel<-animal.normal[,grepl(Zobel1,names(animal.normal))]
Zobel<-Zobel[,grepl("LV",names(Zobel))]
Zobel<-as.data.frame(t(Zobel))
d<-data.frame(try=colSums(Zobel))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(Zobel),sequencing.depth=sum(Zobel),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
Zobel1<-paste0(metadata.normal$site[metadata.normal$method=="Zobel"],collapse="|")
Zobel<-animal.normal[,grepl(Zobel1,names(animal.normal))]
Zobel<-Zobel[,grepl("LW",names(Zobel))]
Zobel<-as.data.frame(t(Zobel))
d<-data.frame(try=colSums(Zobel))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(Zobel),sequencing.depth=sum(Zobel),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
Zobel<-shannon.normal2
Zobel$method<-"Zobel"
##LUCAS
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-animal.normal[,grepl(LUCAS1,names(animal.normal))]
LUCAS<-LUCAS[,grepl("LZ",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))

d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-animal.normal[,grepl(LUCAS1,names(animal.normal))]
LUCAS<-LUCAS[,grepl("LV",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))
d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-animal.normal[,grepl(LUCAS1,names(animal.normal))]
LUCAS<-LUCAS[,grepl("LW",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))
d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
LUCAS<-shannon.normal2
LUCAS$method<-"LUCAS"
##SUCC
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-animal.normal[,grepl(SUCC1,names(animal.normal))]
SUCC<-SUCC[,grepl("LZ",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))

d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-animal.normal[,grepl(SUCC1,names(animal.normal))]
SUCC<-SUCC[,grepl("LV",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))
d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-animal.normal[,grepl(SUCC1,names(animal.normal))]
SUCC<-SUCC[,grepl("LW",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))
d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
SUCC<-shannon.normal2
SUCC$method<-"SUCC"
#deep
###get the unpooled information
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-animal.normal2[,grepl(Deep1,names(animal.normal2))]
Deep1<-Deep1[,grepl("LZ",names(Deep1))]

animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)

Deep2<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
animal.normal$OTU<-row.names(animal.normal)
Deep2<-animal.normal[,grepl("LZA1|LZA7|OTU",names(animal.normal))]
Deep1<-select(Deep1,-c("LZA1","LZA7"))
Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,Deep2,by="OTU")
Deep<-select(Deep,-"OTU")
Deep<-Deep[,!grepl("B",names(Deep))]
deep<-Deep[,grepl("LZ",names(Deep))]
deep<-as.data.frame(t(deep))
d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-animal.normal2[,grepl(Deep1,names(animal.normal2))]
Deep1<-Deep1[,grepl("LV",names(Deep1))]

animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
Deep2<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
animal.normal$OTU<-row.names(animal.normal)
Deep2<-animal.normal[,grepl("LVA1|OTU",names(animal.normal))]
Deep1<-select(Deep1,-c("LVA1"))
Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,Deep2,by="OTU")
Deep<-select(Deep,-"OTU")

Deep<-Deep[,grepl("LV",names(Deep))]
deep<-Deep[,!grepl("B",names(Deep))]

deep<-as.data.frame(t(deep))
b<-as.data.frame(renyi(deep,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=1,sequencing.depth=rowSums(deep),sample=row.names(as.data.frame(renyi(deep,scales=1,hill = TRUE))))

d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-animal.normal2[,grepl(Deep1,names(animal.normal2))]
Deep1<-Deep1[,grepl("LW",names(Deep1))]

Deep<-Deep1
Deep<-as.data.frame(Deep[,!grepl("B",names(Deep))])
deep<-as.data.frame(t(Deep))
b<-as.data.frame(renyi(deep,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=1,sequencing.depth=rowSums(deep),sample=row.names(as.data.frame(renyi(deep,scales=1,hill = TRUE))))
shannon.normal$type<-"unpooled"

##
shannon.normal$site<-"LW"
deep<-shannon.normal
deep$method<-"deep"
##GSMc_62
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-animal.normal[,grepl(GSMc_621,names(animal.normal))]
GSMc_62<-GSMc_62[,grepl("LZ",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")

shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-animal.normal[,grepl(GSMc_621,names(animal.normal))]
GSMc_62<-GSMc_62[,grepl("LV",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")
shannon.normal$type<-"unpooled"
#
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-animal.normal[,grepl(GSMc_621,names(animal.normal))]
GSMc_62<-GSMc_62[,grepl("LW",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")
shannon.normal$type<-"unpooled"
shannon.normal$site<-"LW"
shannon.normal2$type<-"unpooled"
shannon.normal<-rbind(shannon.normal,shannon.normal2)
##
GSMc_62<-shannon.normal
GSMc_62$method<-"GSMc_62"
##GSMc_40A
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-animal.normal[,grepl(GSMc40A,names(animal.normal))]
GSMc_40A<-GSMc_40A[,grepl("LZ",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LV",select$X)
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-animal.normal[,grepl(GSMc40A,names(animal.normal))]
GSMc_40A<-GSMc_40A[,grepl("LV",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LW",select$X)
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-animal.normal[,grepl(GSMc40A,names(animal.normal))]
GSMc_40A<-GSMc_40A[,grepl("LW",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
GSMc_40A<-shannon.normal2
GSMc_40A$method<-"GSMc_40A"
##GSMc_40B
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
select<-read.csv("sheet2_for_pooled.csv")
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-animal.normal[,grepl(GSMc40B,names(animal.normal))]
GSMc_40B<-GSMc_40B[,grepl("LZ",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LV",select$X)
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-animal.normal[,grepl(GSMc40B,names(animal.normal))]
GSMc_40B<-GSMc_40B[,grepl("LV",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LW",select$X)
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-animal.normal[,grepl(GSMc40B,names(animal.normal))]
GSMc_40B<-GSMc_40B[,grepl("LW",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
GSMc_40B<-shannon.normal2
GSMc_40B$method<-"GSMc_40B"
##deep_SUCC
###get the unpooled information
##deep_SUCC
SUCC0<-paste(c(metadata.normal$site[metadata.normal$method=="SUCC"],"OTU"),collapse = "|")
SUCC1<-animal.normal[,grepl(SUCC0,names(animal.normal))]
SUCC1$OTU<-row.names(SUCC1)
###get the unrarefied SUCC
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

SUCC2<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC2<-animal.normal2[,grepl(SUCC2,names(animal.normal2))]
SUCC2<-SUCC2[,grepl("LW39|LV32",names(SUCC2))]
SUCC2$OTU<-row.names(SUCC2)
SUCC1<-merge(SUCC1,SUCC2,by="OTU",all=T)
SUCC1[is.na(SUCC1)]<-0
#####
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-animal.normal2[,grepl(Deep1,names(animal.normal2))]
Deep1<-Deep1[,grepl("LZ",names(Deep1))]

Deep2<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
animal.normal$OTU<-row.names(animal.normal)
Deep2<-animal.normal[,grepl("LZA1|LZA7|OTU",names(animal.normal))]
Deep1<-select(Deep1,-c("LZA1","LZA7"))
Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,Deep2,by="OTU")

Deep<-Deep[,!grepl("B",names(Deep))]
deep1<-Deep
try2<-merge(deep1,SUCC1,by="OTU",all=T)
try2[is.na(try2)]<-0
ds2<-try2
###for deep in LV
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-animal.normal2[,grepl(Deep1,names(animal.normal2))]
Deep1<-Deep1[,grepl("LV",names(Deep1))]

Deep2<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
animal.normal$OTU<-row.names(animal.normal)
Deep2<-animal.normal[,grepl("LVA1|OTU",names(animal.normal))]
Deep1<-select(Deep1,-c("LVA1"))
Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,Deep2,by="OTU")

try2<-merge(ds2,Deep,by="OTU",all=T)
try2[is.na(try2)]<-0
ds2<-try2

###for deep in LW
animal.normal1<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
tax<-animal.normal1[,1:22]
animal.normal1<- animal.normal1[animal.normal1$kBP>0.8,]
animal.normal1<-animal.normal1[,c(1,23:355)]
names(animal.normal1)<-gsub(".*_","",names(animal.normal1))
names(animal.normal1)<-gsub("mix","",names(animal.normal1))

row.names(animal.normal1)<-animal.normal1$GlobalESV
animal.normal1<-animal.normal1[,c(2:334)]
animal.normal1$Feature.ID<-row.names(animal.normal1)

animal.normal1<-animal.normal1[colSums(animal.normal1[,c(1:333)])>100]
animal.normal1<-animal.normal1[rowSums(animal.normal1[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal1,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
row.names(animal.normal2)<-animal.normal2$Feature.ID

Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-animal.normal2[,grepl(Deep1,names(animal.normal2))]
Deep1<-Deep1[,grepl("LW",names(Deep1))]
Deep1$OTU<-row.names(Deep1)

try2<-merge(ds2,Deep1,by="OTU",all=T)
try2[is.na(try2)]<-0
ds2<-try2
Deep<-ds2[,!grepl("B",names(ds2))]
ds2<-Deep
group<-c("LZA8","LZ38", "LZA1", "LZ31", "LZA7", "LZ37",
         "LVA2","LV32", "LVA4","LV34", "LVA5","LV35",
         "LVA8","LV38", "LVA9", "LV39","LVA1","LV31", "LWA9", "LW39")
library(dplyr)
a<-select(ds2,group)
for(i in 1:(length(group)/2)){
  a[,(length(group)+i)]<-rowSums(a[,group[((2*i)-1):(2*i)]])
  names(a)[ncol(a)]<-paste0(substr(group[((2*i)-1)],1,2),i)
}
deep_SUCC<-a[,c((length(group)+1):ncol(a))]
deep_SUCC1<-deep_SUCC[,grepl("LZ",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
deep_SUCC1<-deep_SUCC[,grepl("LV",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##LW
deep_SUCC1<-deep_SUCC[,grepl("LW",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
b<-as.data.frame(renyi(deep_SUCC1,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=1,sequencing.depth=rowSums(deep_SUCC1),sample=row.names(as.data.frame(renyi(deep_SUCC1,scales=1,hill = TRUE))))
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
deep_SUCC<-shannon.normal2
deep_SUCC$method<-"deep_SUCC"

animal.shannon<-rbind(deep,deep_SUCC,GSMc_40A,GSMc_40B,GSMc_62,DarkDiv,SUCC,Zobel,LUCAS)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
write.csv(animal.shannon,"unpooled.animal.shannon.exp.atsummarysequencing.depth.csv")

##bacteria
library(data.table)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/bacteria")
bacteria.soil<-read.csv("bacteria_soil.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.soil)<-gsub("mix","",names(bacteria.soil))
bacteria.soil$Feature.ID<-row.names(bacteria.soil)
bacteria.soil<-bacteria.soil[rowSums(bacteria.soil[,c(1:36)])!=0,]
bacteria.soil<-bacteria.soil[bacteria.soil$Confidence > 0.8,]
names(bacteria.soil)[1:36]<-paste0(names(bacteria.soil)[1:36],"soil")

##DNA
library(data.table)
bacteria.DNA<-read.csv("bacteria_DNA.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.DNA)<-gsub("mix","",names(bacteria.DNA))
bacteria.DNA$Feature.ID<-row.names(bacteria.DNA)
bacteria.DNA<-bacteria.DNA[rowSums(bacteria.DNA[,c(1:36)])!=0,]
bacteria.DNA<-bacteria.DNA[bacteria.DNA$Confidence > 0.8,]
names(bacteria.DNA)[1:36]<-paste0(names(bacteria.DNA)[1:36],"DNA")
mer<-function(x,y){
  merge(x,y,by="Feature.ID",all=T)
}

bacteria<-list(bacteria.DNA,bacteria.soil)
bacteria2<-Reduce(mer,bacteria)
row.names(bacteria2)<-bacteria2$Feature.ID
bacteria2<-bacteria2[,grepl("LZ|LV|LW",names(bacteria2))]
bacteria2[is.na(bacteria2)]<-0

metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")
metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")
metadata2<-rbind(rbind(metadata.DNA,metadata.soil))

###
##deep_SUCC
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
deep0<-(metadata.normal$site[metadata.normal$method=="Deep"])
deep0<-deep0[!grepl("B",deep0)]
deep0<-paste(deep0,collapse = "|")
SUCC0<-paste(c(metadata.normal$site[metadata.normal$method=="SUCC"],"OTU"),collapse = "|")
SUCC1<-bacteria.normal[,grepl(SUCC0,names(bacteria.normal))]
deep1<-bacteria.normal[,grepl(deep0,names(bacteria.normal))]
deep1$OTU<-row.names(deep1)
SUCC1$OTU<-row.names(SUCC1)

try2<-merge(deep1,SUCC1,by="OTU",all=T)
try2[is.na(try2)]<-0
ds2<-try2

group<-c("LZ32", "LZA2", "LZ33", "LZA3", "LZ37", "LZA7", "LW33", "LWA3", "LW32", "LWA2", "LV32", "LVA2",
         "LV33", "LVA3", "LV39", "LVA9", "LZ31", "LZA1", "LZ34", "LZA4",
         "LZ35", "LZA5", "LZ36", "LZA6", "LZ38", "LZA8", "LZ39", "LZA9", "LW31", "LWA1", "LW34", "LWA4",
         "LW35", "LWA5", "LW36", "LWA6", "LW37", "LWA7", "LW38", "LWA8",
         "LW39", "LWA9", "LV31", "LVA1", "LV34", "LVA4", "LV35", "LVA5", "LV36", "LVA6", "LV37", "LVA7", "LV38", "LVA8")
library(dplyr)
a<-select(ds2,group)
for(i in 1:(length(group)/2)){
  a[,(length(group)+i)]<-rowSums(a[,group[((2*i)-1):(2*i)]])
  names(a)[ncol(a)]<-paste0(substr(group[((2*i)-1)],1,2),i)
}
deep_SUCC<-a[,c((length(group)+1):ncol(a))]
###
deep_SUCC1<-deep_SUCC[,grepl("LZ",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))

d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")

shannon.normal$type<-"unpooled"
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
#LV
deep_SUCC1<-deep_SUCC[,grepl("LV",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
deep_SUCC1<-deep_SUCC[,grepl("LW",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")

shannon.normal$type<-"unpooled"
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
deep_SUCC<-shannon.normal2
deep_SUCC$method<-"deep_SUCC"
##DarkDiv
###get the unpooled information
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)

DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-bacteria.normal[,grepl(DarkDiv1,names(bacteria.normal))]
DarkDiv<-DarkDiv[,grepl("LZ",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
library(vegan)
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)

DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-bacteria.normal[,grepl(DarkDiv1,names(bacteria.normal))]
DarkDiv<-DarkDiv[,grepl("LV",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)

DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-bacteria.normal[,grepl(DarkDiv1,names(bacteria.normal))]
DarkDiv<-DarkDiv[,grepl("LW",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
DarkDiv<-shannon.normal2
DarkDiv$method<-"DarkDiv"
##LUCAS
###get the unpooled information
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-bacteria.normal[,grepl(LUCAS1,names(bacteria.normal))]
LUCAS<-LUCAS[,grepl("LZ",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))

d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-bacteria.normal[,grepl(LUCAS1,names(bacteria.normal))]
LUCAS<-LUCAS[,grepl("LV",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))
d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-bacteria.normal[,grepl(LUCAS1,names(bacteria.normal))]
LUCAS<-LUCAS[,grepl("LW",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))
d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
LUCAS<-shannon.normal2
LUCAS$method<-"LUCAS"
##SUCC
###get the unpooled information


SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-bacteria.normal[,grepl(SUCC1,names(bacteria.normal))]
SUCC<-SUCC[,grepl("LZ",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))

d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information


SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-bacteria.normal[,grepl(SUCC1,names(bacteria.normal))]
SUCC<-SUCC[,grepl("LV",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))
d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-bacteria.normal[,grepl(SUCC1,names(bacteria.normal))]
SUCC<-SUCC[,grepl("LW",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))
d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
SUCC<-shannon.normal2
SUCC$method<-"SUCC"
##deep
###get the unpooled information
deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
deep<-bacteria.normal[,grepl(deep1,names(bacteria.normal))]
deep<-deep[,grepl("LZ",names(deep))]
deep<-deep[,!grepl("B",names(deep))]
deep<-as.data.frame(t(deep))

d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
deep<-bacteria.normal[,grepl(deep1,names(bacteria.normal))]
deep<-deep[,grepl("LV",names(deep))]
deep<-deep[,!grepl("B",names(deep))]
deep<-as.data.frame(t(deep))
d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
deep<-bacteria.normal[,grepl(deep1,names(bacteria.normal))]
deep<-deep[,grepl("LW",names(deep))]
deep<-deep[,!grepl("B",names(deep))]
deep<-as.data.frame(t(deep))
d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
deep<-shannon.normal2
deep$method<-"deep"
##GSMc_62
###get the unpooled information
GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-bacteria.normal[,grepl(GSMc_621,names(bacteria.normal))]
GSMc_62<-GSMc_62[,grepl("LZ",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")

shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information


GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-bacteria.normal[,grepl(GSMc_621,names(bacteria.normal))]
GSMc_62<-GSMc_62[,grepl("LV",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")
shannon.normal$type<-"unpooled"
#
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information


GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-bacteria.normal[,grepl(GSMc_621,names(bacteria.normal))]
GSMc_62<-GSMc_62[,grepl("LW",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")
shannon.normal$type<-"unpooled"
shannon.normal$site<-"LW"
shannon.normal2$type<-"unpooled"
shannon.normal<-rbind(shannon.normal,shannon.normal2)
##
GSMc_62<-shannon.normal
GSMc_62$method<-"GSMc_62"
##GSMc_40A
###get the unpooled information

select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-bacteria.normal[,grepl(GSMc40A,names(bacteria.normal))]
GSMc_40A<-GSMc_40A[,grepl("LZ",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information

select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LV",select$X)
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-bacteria.normal[,grepl(GSMc40A,names(bacteria.normal))]
GSMc_40A<-GSMc_40A[,grepl("LV",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information


select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LW",select$X)
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-bacteria.normal[,grepl(GSMc40A,names(bacteria.normal))]
GSMc_40A<-GSMc_40A[,grepl("LW",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
GSMc_40A<-shannon.normal2
GSMc_40A$method<-"GSMc_40A"
##GSMc_40B
###get the unpooled information

select<-read.csv("sheet2_for_pooled.csv")
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-bacteria.normal[,grepl(GSMc40B,names(bacteria.normal))]
GSMc_40B<-GSMc_40B[,grepl("LZ",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information

select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LV",select$X)
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-bacteria.normal[,grepl(GSMc40B,names(bacteria.normal))]
GSMc_40B<-GSMc_40B[,grepl("LV",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information


select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LW",select$X)
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-bacteria.normal[,grepl(GSMc40B,names(bacteria.normal))]
GSMc_40B<-GSMc_40B[,grepl("LW",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
GSMc_40B<-shannon.normal2
GSMc_40B$method<-"GSMc_40B"
##

bacteria.shannon<-rbind(deep,deep_SUCC,GSMc_40A,GSMc_40B,GSMc_62,DarkDiv,SUCC,Zobel,LUCAS)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/2.compare mean and total richness of unpooled methods")
write.csv(bacteria.shannon,"unpooled.bacteria.shannon.exp.atsummarysequencing.depth.csv")

##fungi
library(data.table)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/fungi")
fungi.soil<-read.csv("fungi_soil.csv",header = TRUE,sep = ",",row.names = 1)
fungi.soil<-fungi.soil[fungi.soil$Confidence > 0.8,]
names(fungi.soil)<-gsub("mix","",names(fungi.soil))

fungi.soil<-fungi.soil[,c(1:35)]
fungi.soil$Feature.ID<-row.names(fungi.soil)
fungi.soil<-fungi.soil[rowSums(fungi.soil[,c(1:35)])!=0,]

tax<-read.table("original result/taxonomy.tsv",header = TRUE,sep = "\t")
fungi.soil2<-merge(fungi.soil,tax,by="Feature.ID",all=T)
fungi.soil2<-fungi.soil2[!is.na(fungi.soil2$LZ9B),]
names(fungi.soil2)[2:36]<-paste0(names(fungi.soil2)[2:36],"soil")
tax_table.soil<-fungi.soil2[,c(37:38)]
row.names(tax_table.soil)<-fungi.soil2$Feature.ID
tax_table.soil$Feature.ID<-row.names(tax_table.soil)
##DNA
library(data.table)
fungi.DNA<-read.csv("fungi_DNA.csv",header = TRUE,sep = ",",row.names = 1)
names(fungi.DNA)<-gsub("mix","",names(fungi.DNA))
fungi.DNA<-fungi.DNA[fungi.DNA$Confidence > 0.8,]

fungi.DNA$Feature.ID<-row.names(fungi.DNA)
fungi.DNA<-fungi.DNA[rowSums(fungi.DNA[,c(1:36)])!=0,]
names(fungi.DNA)[1:36]<-paste0(names(fungi.DNA)[1:36],"DNA")
tax_table.DNA<-fungi.DNA[,c(37:39)]
row.names(tax_table.DNA)<-fungi.DNA$Feature.ID

mer<-function(x,y){
  merge(x,y,by="Feature.ID",all=T)
}

fungi<-list(fungi.DNA,fungi.soil2)
fungi2<-Reduce(mer,fungi)

row.names(fungi2)<-fungi2$Feature.ID
fungi2<-fungi2[,grepl("LZ|LV|LW",names(fungi2))]

fungi2[is.na(fungi2)]<-0

metadata<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))
##DarkDiv
###get the unpooled information
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-fungi.normal[,grepl(DarkDiv1,names(fungi.normal))]
DarkDiv<-DarkDiv[,grepl("LZ",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
library(vegan)
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-fungi.normal[,grepl(DarkDiv1,names(fungi.normal))]
DarkDiv<-DarkDiv[,grepl("LV",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
DarkDiv1<-paste0(metadata.normal$site[metadata.normal$method=="DarkDiv"],collapse="|")
DarkDiv<-fungi.normal[,grepl(DarkDiv1,names(fungi.normal))]
DarkDiv<-DarkDiv[,grepl("LW",names(DarkDiv))]
DarkDiv<-as.data.frame(t(DarkDiv))
d<-data.frame(try=colSums(DarkDiv))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(DarkDiv),sequencing.depth=sum(DarkDiv),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
DarkDiv<-shannon.normal2
DarkDiv$method<-"DarkDiv"
##LUCAS
###get the unpooled information
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-fungi.normal[,grepl(LUCAS1,names(fungi.normal))]
LUCAS<-LUCAS[,grepl("LZ",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))

d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
fungi.normal1<-read.csv("fungi_normal.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal1<-fungi.normal1[fungi.normal1$Confidence > 0.8,]
fungi.normal1$Feature.ID<-row.names(fungi.normal1)
fungi.normal1<-fungi.normal1[rowSums(fungi.normal1[,c(1:331)])!=0,]
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS1<-fungi.normal1[,grepl(LUCAS1,names(fungi.normal1))]
LUCAS1<-LUCAS1[,grepl("LV",names(LUCAS1))]

fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
LUCAS2<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
fungi.normal$OTU<-row.names(fungi.normal)
LUCAS2<-fungi.normal[,grepl("LV65|OTU",names(fungi.normal))]
LUCAS1<-select(LUCAS1,-"LV65")
LUCAS1$OTU<-row.names(LUCAS1)
LUCAS<-merge(LUCAS1,LUCAS2,by="OTU")
LUCAS<-select(LUCAS,-"OTU")
LUCAS<-LUCAS[,grepl("LV",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))
d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
LUCAS1<-paste0(metadata.normal$site[metadata.normal$method=="LUCAS"],collapse="|")
LUCAS<-fungi.normal[,grepl(LUCAS1,names(fungi.normal))]
LUCAS<-LUCAS[,grepl("LW",names(LUCAS))]
LUCAS<-as.data.frame(t(LUCAS))
d<-data.frame(try=colSums(LUCAS))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(LUCAS),sequencing.depth=sum(LUCAS),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
LUCAS<-shannon.normal2
LUCAS$method<-"LUCAS"
##SUCC
###get the unpooled information
SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-fungi.normal[,grepl(SUCC1,names(fungi.normal))]
SUCC<-SUCC[,grepl("LZ",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))

d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information


SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-fungi.normal[,grepl(SUCC1,names(fungi.normal))]
SUCC<-SUCC[,grepl("LV",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))
d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
SUCC1<-paste0(metadata.normal$site[metadata.normal$method=="SUCC"],collapse="|")
SUCC<-fungi.normal[,grepl(SUCC1,names(fungi.normal))]
SUCC<-SUCC[,grepl("LW",names(SUCC))]
SUCC<-as.data.frame(t(SUCC))
d<-data.frame(try=colSums(SUCC))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(SUCC),sequencing.depth=sum(SUCC),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
SUCC<-shannon.normal2
SUCC$method<-"SUCC"
##deep
###get the unpooled information
deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
deep<-fungi.normal[,grepl(deep1,names(fungi.normal))]
deep<-deep[,grepl("LZ",names(deep))]
deep<-deep[,!grepl("B",names(deep))]
deep<-as.data.frame(t(deep))

d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
fungi.normal1<-read.csv("fungi_normal.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal1<-fungi.normal1[fungi.normal1$Confidence > 0.8,]
fungi.normal1$Feature.ID<-row.names(fungi.normal1)
fungi.normal1<-fungi.normal1[rowSums(fungi.normal1[,c(1:331)])!=0,]
Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-fungi.normal1[,grepl(Deep1,names(fungi.normal1))]
Deep1<-Deep1[,grepl("LV",names(Deep1))]

fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
Deep2<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
fungi.normal$OTU<-row.names(fungi.normal)
Deep2<-fungi.normal[,grepl("LVA8|LVA9|OTU",names(fungi.normal))]
Deep1<-select(Deep1,-c("LVA8","LVA9"))
Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,Deep2,by="OTU")
Deep<-select(Deep,-"OTU")

deep<-Deep[,grepl("LV",names(Deep))]
deep<-deep[,!grepl("B",names(deep))]
deep<-as.data.frame(t(deep))

d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
deep<-fungi.normal[,grepl(deep1,names(fungi.normal))]
deep<-deep[,grepl("LW",names(deep))]
deep<-deep[,!grepl("B",names(deep))]
deep<-as.data.frame(t(deep))
d<-data.frame(try=colSums(deep))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep),sequencing.depth=sum(deep),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
deep<-shannon.normal2
deep$method<-"deep"
##GSMc_62
###get the unpooled information
GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-fungi.normal[,grepl(GSMc_621,names(fungi.normal))]
GSMc_62<-GSMc_62[,grepl("LZ",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")

shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information


GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-fungi.normal[,grepl(GSMc_621,names(fungi.normal))]
GSMc_62<-GSMc_62[,grepl("LV",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")
shannon.normal$type<-"unpooled"
#
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information


GSMc_621<-paste0(metadata.normal$site[metadata.normal$method=="GSMc"],collapse="|")
GSMc_62<-fungi.normal[,grepl(GSMc_621,names(fungi.normal))]
GSMc_62<-GSMc_62[,grepl("LW",names(GSMc_62))]
GSMc_62<-GSMc_62[,!grepl("B",names(GSMc_62))]
GSMc_62<-as.data.frame(t(GSMc_62))
d<-data.frame(try=colSums(GSMc_62))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_62),sequencing.depth=sum(GSMc_62),sample="try")
shannon.normal$type<-"unpooled"
shannon.normal$site<-"LW"
shannon.normal2$type<-"unpooled"
shannon.normal<-rbind(shannon.normal,shannon.normal2)
##
GSMc_62<-shannon.normal
GSMc_62$method<-"GSMc_62"
##GSMc_40A
###get the unpooled information

select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-fungi.normal[,grepl(GSMc40A,names(fungi.normal))]
GSMc_40A<-GSMc_40A[,grepl("LZ",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information

select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LV",select$X)
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-fungi.normal[,grepl(GSMc40A,names(fungi.normal))]
GSMc_40A<-GSMc_40A[,grepl("LV",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LW",select$X)
GSMc40A<-paste0(select[select$X1Amix=="GSMc40A",]$X,collapse = "|")
GSMc_40A<-fungi.normal[,grepl(GSMc40A,names(fungi.normal))]
GSMc_40A<-GSMc_40A[,grepl("LW",names(GSMc_40A))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
d<-data.frame(try=colSums(GSMc_40A))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40A),sequencing.depth=sum(GSMc_40A),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
GSMc_40A<-shannon.normal2
GSMc_40A$method<-"GSMc_40A"
##GSMc_40B
###get the unpooled information
select<-read.csv("sheet2_for_pooled.csv")
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-fungi.normal[,grepl(GSMc40B,names(fungi.normal))]
GSMc_40B<-GSMc_40B[,grepl("LZ",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LV",select$X)
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-fungi.normal[,grepl(GSMc40B,names(fungi.normal))]
GSMc_40B<-GSMc_40B[,grepl("LV",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
select<-read.csv("sheet2_for_pooled.csv")
select$X<-gsub("LZ","LW",select$X)
GSMc40B<-paste0(select[select$X1Bmix=="GSMc40B",]$X,collapse = "|")
GSMc_40B<-fungi.normal[,grepl(GSMc40B,names(fungi.normal))]
GSMc_40B<-GSMc_40B[,grepl("LW",names(GSMc_40B))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
d<-data.frame(try=colSums(GSMc_40B))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(GSMc_40B),sequencing.depth=sum(GSMc_40B),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
GSMc_40B<-shannon.normal2
GSMc_40B$method<-"GSMc_40B"
#deep_SUCC
##deep_SUCC
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
SUCC0<-paste(c(metadata.normal$site[metadata.normal$method=="SUCC"],"OTU"),collapse = "|")
SUCC1<-fungi.normal[,grepl(SUCC0,names(fungi.normal))]
SUCC1$OTU<-row.names(SUCC1)

fungi.normal1<-read.csv("fungi_normal.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal1<-fungi.normal1[fungi.normal1$Confidence > 0.8,]
fungi.normal1$Feature.ID<-row.names(fungi.normal1)
fungi.normal1<-fungi.normal1[rowSums(fungi.normal1[,c(1:331)])!=0,]
Deep1<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
Deep1<-fungi.normal1[,grepl(Deep1,names(fungi.normal1))]
Deep1<-Deep1[,grepl("LV",names(Deep1))]

Deep2<-paste0(metadata.normal$site[metadata.normal$method=="Deep"],collapse="|")
fungi.normal$OTU<-row.names(fungi.normal)
Deep2<-fungi.normal[,grepl("LVA8|LVA9|OTU",names(fungi.normal))]
Deep1<-select(Deep1,-c("LVA8","LVA9"))
Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,Deep2,by="OTU")
Deep<-Deep[,!grepl("B",names(Deep))]
deep1<-Deep
try2<-merge(deep1,SUCC1,by="OTU",all=T)
try2[is.na(try2)]<-0
ds2<-try2

Deep1<-metadata.normal$site[metadata.normal$method=="Deep"]
Deep1<-Deep1[!grepl("B",Deep1)]
Deep1<-paste0(Deep1,collapse="|")
Deep<-fungi.normal[,grepl(Deep1,names(fungi.normal))]
Deep1<-Deep[,grepl("LW|LZ",names(Deep))]

Deep1$OTU<-row.names(Deep1)
Deep<-merge(Deep1,ds2,by="OTU")
Deep<-Deep[,!grepl("B",names(Deep))]

ds2<-Deep
group<-c("LZ32", "LZA2",  "LZ37", "LZA7", "LW33", "LWA3", "LW32", "LWA2", 
         "LV39", "LVA9",  
         "LZ36", "LZA6", "LZ38", "LZA8", "LZ39", "LZA9", "LW31", "LWA1", "LW34", "LWA4",
         "LW35", "LWA5", "LW36", "LWA6",  "LW38", "LWA8",
         "LW39", "LWA9", "LV31", "LVA1",  "LV35", "LVA5", "LV36", "LVA6",  "LV38", "LVA8")
library(dplyr)
a<-select(ds2,group)
for(i in 1:(length(group)/2)){
  a[,(length(group)+i)]<-rowSums(a[,group[((2*i)-1):(2*i)]])
  names(a)[ncol(a)]<-paste0(substr(group[((2*i)-1)],1,2),i)
}
deep_SUCC<-a[,c((length(group)+1):ncol(a))]

###
deep_SUCC1<-deep_SUCC[,grepl("LZ",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")
shannon.normal$type<-"unpooled"

##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
#LV
deep_SUCC1<-deep_SUCC[,grepl("LV",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
deep_SUCC1<-deep_SUCC[,grepl("LW",names(deep_SUCC))]
deep_SUCC1<-as.data.frame(t(deep_SUCC1))
d<-data.frame(try=colSums(deep_SUCC1))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(deep_SUCC1),sequencing.depth=sum(deep_SUCC1),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
deep_SUCC<-shannon.normal2
deep_SUCC$method<-"deep_SUCC"
##Zobel
###get the unpooled information
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
Zobel1<-paste0(metadata.normal$site[metadata.normal$method=="Zobel"],collapse="|")
Zobel<-fungi.normal[,grepl(Zobel1,names(fungi.normal))]
Zobel<-Zobel[,grepl("LZ",names(Zobel))]
Zobel<-as.data.frame(t(Zobel))
d<-data.frame(try=colSums(Zobel))
d<-as.data.frame(t(d))
library(vegan)
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(Zobel),sequencing.depth=sum(Zobel),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal2<-shannon.normal
shannon.normal2$site<-"LZ"
##LV
###get the unpooled information
Zobel1<-paste0(metadata.normal$site[metadata.normal$method=="Zobel"],collapse="|")
Zobel<-fungi.normal[,grepl(Zobel1,names(fungi.normal))]
Zobel<-Zobel[,grepl("LV",names(Zobel))]
Zobel<-as.data.frame(t(Zobel))
d<-data.frame(try=colSums(Zobel))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(Zobel),sequencing.depth=sum(Zobel),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LV"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
#LW
###get the unpooled information
Zobel1<-paste0(metadata.normal$site[metadata.normal$method=="Zobel"],collapse="|")
Zobel<-fungi.normal[,grepl(Zobel1,names(fungi.normal))]
Zobel<-Zobel[,grepl("LW",names(Zobel))]
Zobel<-as.data.frame(t(Zobel))
d<-data.frame(try=colSums(Zobel))
d<-as.data.frame(t(d))
b<-as.data.frame(renyi(d,scales=1,hill = TRUE))
names(b)<-"shannon"
shannon.normal<-data.frame(shannon=b$shannon,number=nrow(Zobel),sequencing.depth=sum(Zobel),sample="try")
shannon.normal$type<-"unpooled"
##
shannon.normal$site<-"LW"
shannon.normal2<-rbind(shannon.normal,shannon.normal2)
##
Zobel<-shannon.normal2
Zobel$method<-"Zobel"
fungi.shannon<-rbind(deep,deep_SUCC,GSMc_40A,GSMc_40B,GSMc_62,DarkDiv,SUCC,Zobel,LUCAS)
write.csv(fungi.shannon,"unpooled.fungi.shannon.exp.atsummarysequencing.depth.csv")
