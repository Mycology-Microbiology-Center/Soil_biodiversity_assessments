###get the rarefied table
##bacteria
library(data.table)
library(metagMisc)
library(phyloseq)
bacteria.soil<-read.csv("bacteria_soil.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.soil)<-gsub("mix","",names(bacteria.soil))
bacteria.soil$Feature.ID<-row.names(bacteria.soil)
bacteria.soil<-bacteria.soil[rowSums(bacteria.soil[,c(1:36)])!=0,]

bacteria.soil<-bacteria.soil[bacteria.soil$Confidence > 0.8,]
names(bacteria.soil)[1:36]<-paste0(names(bacteria.soil)[1:36],"soil")
tax_table.soil<-bacteria.soil[,c(37:39)]
##normal
bacteria.normal<-read.csv("bacteria_normal.csv",header = TRUE,sep = ",",row.names = 1)
bacteria.normal<-bacteria.normal[bacteria.normal$Confidence > 0.8,]
tax_table.normal<-bacteria.normal[,c(335:336)]
row.names(tax_table.normal)<-bacteria.normal$OTU.ID
tax_table.normal$Feature.ID<-row.names(tax_table.normal)

row.names(bacteria.normal)<-bacteria.normal$OTU.ID
bacteria.normal<-bacteria.normal[,2:334]
bacteria.normal2<-bacteria.normal[,colSums(bacteria.normal)>100]
names(bacteria.normal2)<-gsub(".*_","",names(bacteria.normal2))

##DNA
bacteria.DNA<-read.csv("bacteria_DNA.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.DNA)<-gsub("mix","",names(bacteria.DNA))
bacteria.DNA$Feature.ID<-row.names(bacteria.DNA)
bacteria.DNA<-bacteria.DNA[rowSums(bacteria.DNA[,c(1:36)])!=0,]
bacteria.DNA<-bacteria.DNA[bacteria.DNA$Confidence > 0.8,]
tax_table.DNA<-bacteria.DNA[,c(37:39)]
names(bacteria.DNA)[1:36]<-paste0(names(bacteria.DNA)[1:36],"DNA")
mer<-function(x,y){
  merge(x,y,by="Feature.ID",all=T)
}

bacteria.normal2$Feature.ID<-row.names(bacteria.normal2)
bacteria<-list(bacteria.normal2,bacteria.DNA,bacteria.soil)
bacteria2<-Reduce(mer,bacteria)

row.names(bacteria2)<-bacteria2$Feature.ID
bacteria2<-bacteria2[,grepl("LZ|LV|LW",names(bacteria2))]

bacteria2[is.na(bacteria2)]<-0

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil),metadata[,-5])
tax<-rbind(tax_table.soil,tax_table.DNA,tax_table.normal)
tax<-tax[!duplicated(tax$Feature.ID),]
row.names(tax)<-tax$Feature.ID
row.names(metadata2)<-metadata2$site
###
bacteria2<-select(bacteria2,-names(sort(colSums(bacteria2))[1:3]))
###
metadata2<-metadata2[match(names(bacteria2),metadata2$site),]
bacteria <- phyloseq(
  otu_table = otu_table(as.matrix(bacteria2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata2)
)
min <- round(min(sample_sums(bacteria)))

d2 <- prune_taxa(taxa_sums(bacteria) > 0, bacteria)
bacteria.normal.rarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,multithread=T,seeds = NULL)
##
for(i in 1:1000){
  bacterianormal.rarefy_diversity[[i]]<-estimate_richness(bacteria.normal.rarefy[[i]], measures = c("Observed","Shannon"))
}

write.table(bacterianormal.rarefy_diversity,"bacteria.rarefy.diversity.csv",sep=",")
bacteria.rarefy.table<-as.data.frame(bacteria.normal.rarefy[[300]]@otu_table)
write.table(bacteria.rarefy.table,"bacteria.rarefy.table.csv",sep=",")

#fungi
fungi.soil<-read.csv("fungi_soil.csv",header = TRUE,sep = ",",row.names = 1)
fungi.soil<-fungi.soil[fungi.soil$Confidence > 0.8,]
names(fungi.soil)<-gsub("mix","",names(fungi.soil))

fungi.soil<-fungi.soil[,c(1:35)]
fungi.soil$Feature.ID<-row.names(fungi.soil)
fungi.soil<-fungi.soil[rowSums(fungi.soil[,c(1:35)])!=0,]

tax<-read.table("taxonomy.tsv",header = TRUE,sep = "\t")
fungi.soil2<-merge(fungi.soil,tax,by="Feature.ID",all=T)
fungi.soil2<-fungi.soil2[!is.na(fungi.soil2$LZ9B),]
names(fungi.soil2)[2:36]<-paste0(names(fungi.soil2)[2:36],"soil")
tax_table.soil<-fungi.soil2[,c(37:38)]
row.names(tax_table.soil)<-fungi.soil2$Feature.ID
tax_table.soil$Feature.ID<-row.names(tax_table.soil)
#normal
fungi.normal<-read.csv("fungi_normal.csv",header = TRUE,sep = ",",row.names = 1)
fungi.normal<-fungi.normal[fungi.normal$Confidence > 0.8,]

fungi.normal$Feature.ID<-row.names(fungi.normal)
fungi.normal<-fungi.normal[rowSums(fungi.normal[,c(1:331)])!=0,]
tax_table.normal<-fungi.normal[,c(332:334)]
row.names(tax_table.normal)<-fungi.normal$Feature.ID

##DNA
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

fungi.normal$Feature.ID<-row.names(fungi.normal)
fungi<-list(fungi.normal,fungi.DNA,fungi.soil2)
fungi2<-Reduce(mer,fungi)

row.names(fungi2)<-fungi2$Feature.ID
fungi2<-fungi2[,grepl("LZ|LV|LW",names(fungi2))]

fungi2[is.na(fungi2)]<-0

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata<-read.csv("normal_metadata2.csv",header = TRUE,sep = ",")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil),metadata[,-5])
tax<-rbind(tax_table.soil,tax_table.DNA,tax_table.normal)
tax<-tax[!duplicated(tax$Feature.ID),]
row.names(metadata2)<-metadata2$site
###
fungi2<-select(fungi2,-names(sort(colSums(fungi2))[1:61]))
###
metadata2<-metadata2[match(names(fungi2),metadata2$site),]
fungi <- phyloseq(
  otu_table = otu_table(as.matrix(fungi2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata2)
)
min <- round(min(sample_sums(fungi)))

d2 <- prune_taxa(taxa_sums(fungi) > 0, fungi)
fungi.normal.rarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,seeds = NULL)
##
for(i in 1:1000){
  funginormal.rarefy_diversity[[i]]<-estimate_richness(fungi.normal.rarefy[[i]], measures = c("Observed","Shannon"))
}

write.table(funginormal.rarefy_diversity,"fungi.rarefy.diversity2.csv",sep=",")
fungi.rarefy.table<-as.data.frame(fungi.normal.rarefy[[300]]@otu_table)
write.table(fungi.rarefy.table,"fungi.rarefy.table2.csv",sep=",")

###animal
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

#normal
animal.normal<-read.csv("animal_normal2.csv",header = TRUE,sep = ",",row.names = 1)
animal.normal<- animal.normal[animal.normal$kBP>0.8,]
tax<-animal.normal[,1:22]

animal.normal<-animal.normal[,c(1,23:355)]
names(animal.normal)<-gsub(".*_","",names(animal.normal))
names(animal.normal)<-gsub("mix","",names(animal.normal))

row.names(animal.normal)<-animal.normal$GlobalESV
animal.normal<-animal.normal[,c(2:334)]
animal.normal$Feature.ID<-row.names(animal.normal)

animal.normal<-animal.normal[colSums(animal.normal[,c(1:333)])>100]
animal.normal<-animal.normal[rowSums(animal.normal[,c(1:310)])!=0,]

animal.normal2<-merge(animal.normal,tax,by.x="Feature.ID",by.y = "COI_GlobalESV",all=T)
animal.normal2<-animal.normal2[!is.na(animal.normal2$LVB4),]
animal.normal2<-animal.normal2[!grepl("Unassigned",animal.normal2$Kingdom),]
tax_table.normal<-animal.normal2[,c(312:320)]
row.names(tax_table.normal)<-animal.normal2$Feature.ID

mer<-function(x,y){
  merge(x,y,by="Feature.ID",all=T)
}

animal<-list(animal.normal2,animal.DNA2,animal.soil2)
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
tax_table.soil$COI<-row.names(tax_table.soil)
tax_table.DNA$COI<-row.names(tax_table.DNA)
tax_table.normal$COI<-row.names(tax_table.normal)

tax<-rbind(tax_table.soil[,c(1:3,9)],tax_table.DNA[,c(1,4,7,8)],tax_table.normal[,c(1,4,7,10)])
tax<-tax[!duplicated(tax$COI),]
row.names(tax)<-tax$COI
row.names(metadata2)<-metadata2$site
###
animal2<-select(animal2,-names(sort(colSums(animal2))[1:54]))
###
metadata2<-metadata2[match(names(animal2),metadata2$site),]
animal <- phyloseq(
  otu_table = otu_table(as.matrix(animal2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata2)
)
min <- round(min(sample_sums(animal)))

d2 <- prune_taxa(taxa_sums(animal) > 0, animal)
animal.normal.rarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,multithread=T,seeds = NULL)
##
animal.normal.rarefy_diversity<-list()
for(i in 1:1000){
  animal.normal.rarefy_diversity[[i]]<-estimate_richness(animal.normal.rarefy[[i]], measures = c("Observed","Shannon"))
}
write.table(animalnormal.rarefy_diversity,"animal.rarefy.diversity2.csv",sep=",")
animal.rarefy.table<-as.data.frame(animal.normal.rarefy[[300]]@otu_table)
write.table(animal.rarefy.table,"animal.rarefy.table2.csv",sep=",")
###rarefy the pooled table separately
##bacteria
bacteria.soil<-read.csv("bacteria_soil.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.soil)<-gsub("mix","",names(bacteria.soil))
bacteria.soil$Feature.ID<-row.names(bacteria.soil)
bacteria.soil<-bacteria.soil[rowSums(bacteria.soil[,c(1:36)])!=0,]

bacteria.soil<-bacteria.soil[bacteria.soil$Confidence > 0.8,]
names(bacteria.soil)[1:36]<-paste0(names(bacteria.soil)[1:36],"soil")
tax_table.soil<-bacteria.soil[,c(37:39)]

##DNA
bacteria.DNA<-read.csv("bacteria_DNA.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.DNA)<-gsub("mix","",names(bacteria.DNA))
bacteria.DNA$Feature.ID<-row.names(bacteria.DNA)
bacteria.DNA<-bacteria.DNA[rowSums(bacteria.DNA[,c(1:36)])!=0,]
bacteria.DNA<-bacteria.DNA[bacteria.DNA$Confidence > 0.8,]
tax_table.DNA<-bacteria.DNA[,c(37:39)]
names(bacteria.DNA)[1:36]<-paste0(names(bacteria.DNA)[1:36],"DNA")
mer<-function(x,y){
  merge(x,y,by="Feature.ID",all=T)
}

bacteria<-list(bacteria.DNA,bacteria.soil)
bacteria2<-Reduce(mer,bacteria)

row.names(bacteria2)<-bacteria2$Feature.ID
bacteria2<-bacteria2[,grepl("LZ|LV|LW",names(bacteria2))]

bacteria2[is.na(bacteria2)]<-0

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))
tax<-rbind(tax_table.soil,tax_table.DNA)
tax<-tax[!duplicated(tax$Feature.ID),]
row.names(tax)<-tax$Feature.ID
row.names(metadata2)<-metadata2$site
###
metadata2<-metadata2[match(names(bacteria2),metadata2$site),]
bacteria <- phyloseq(
  otu_table = otu_table(as.matrix(bacteria2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata2)
)
min <- round(min(sample_sums(bacteria)))

d2 <- prune_taxa(taxa_sums(bacteria) > 0, bacteria)
bacteria.normal.rarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,seeds = NULL)
##
bacterianormal.rarefy_diversity<-list()
for(i in 1:1000){
  bacterianormal.rarefy_diversity[[i]]<-estimate_richness(bacteria.normal.rarefy[[i]], measures = c("Observed","Shannon"))
}

write.table(bacterianormal.rarefy_diversity,"bacteria.pooled.rarefy.csv",sep=",")
bacteria.rarefy.table<-as.data.frame(bacteria.normal.rarefy[[300]]@otu_table)
write.table(bacteria.rarefy.table,"bacteria.pooled.rarefy.table.csv",sep=",")
#fungi
fungi.soil<-read.csv("fungi_soil.csv",header = TRUE,sep = ",",row.names = 1)
fungi.soil<-fungi.soil[fungi.soil$Confidence > 0.8,]
names(fungi.soil)<-gsub("mix","",names(fungi.soil))

fungi.soil<-fungi.soil[,c(1:35)]
fungi.soil$Feature.ID<-row.names(fungi.soil)
fungi.soil<-fungi.soil[rowSums(fungi.soil[,c(1:35)])!=0,]

tax<-read.table("taxonomy.tsv",header = TRUE,sep = "\t")
fungi.soil2<-merge(fungi.soil,tax,by="Feature.ID",all=T)
fungi.soil2<-fungi.soil2[!is.na(fungi.soil2$LZ9B),]
names(fungi.soil2)[2:36]<-paste0(names(fungi.soil2)[2:36],"soil")
tax_table.soil<-fungi.soil2[,c(37:38)]
row.names(tax_table.soil)<-fungi.soil2$Feature.ID
tax_table.soil$Feature.ID<-row.names(tax_table.soil)
##DNA
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

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))
tax<-rbind(tax_table.soil,tax_table.DNA)
tax<-tax[!duplicated(tax$Feature.ID),]
row.names(metadata2)<-metadata2$site

###
metadata2<-metadata2[match(names(fungi2),metadata2$site),]
fungi <- phyloseq(
  otu_table = otu_table(as.matrix(fungi2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata2)
)
min <- round(min(sample_sums(fungi)))

d2 <- prune_taxa(taxa_sums(fungi) > 0, fungi)
fungi.normal.rarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,seeds = NULL)
##
funginormal.rarefy_diversity<-list()
for(i in 1:1000){
  funginormal.rarefy_diversity[[i]]<-estimate_richness(fungi.normal.rarefy[[i]], measures = c("Observed","Shannon"))
}

write.table(funginormal.rarefy_diversity,"fungi.pooled.rarefy.csv",sep=",")
fungi.rarefy.table<-as.data.frame(fungi.normal.rarefy[[300]]@otu_table)
write.table(fungi.rarefy.table,"fungi.pooled.rarefy.table.csv",sep=",")
###animal
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

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))
tax_table.soil$COI<-row.names(tax_table.soil)
tax_table.DNA$COI<-row.names(tax_table.DNA)

tax<-rbind(tax_table.soil[,c(1:3,9)],tax_table.DNA[,c(1,4,7,8)])
tax<-tax[!duplicated(tax$COI),]
row.names(tax)<-tax$COI
row.names(metadata2)<-metadata2$site
###
animal2<-select(animal2,-names(sort(colSums(animal2))[1:7]))
###
metadata2<-metadata2[match(names(animal2),metadata2$site),]
animal <- phyloseq(
  otu_table = otu_table(as.matrix(animal2), taxa_are_rows = TRUE),
  tax_table = tax_table(as.matrix(tax)),
  sample_data = sample_data(metadata2)
)
min <- round(min(sample_sums(animal)))

d2 <- prune_taxa(taxa_sums(animal) > 0, animal)
animal.normal.rarefy<-phyloseq_mult_raref(d2,trimOTUs = F, replace = T, MinSizeThreshold = min, SampSize = min, iter = 1000,seeds = NULL)
##
animal.normal.rarefy_diversity<-list()
for(i in 1:1000){
  animal.normal.rarefy_diversity[[i]]<-estimate_richness(animal.normal.rarefy[[i]], measures = c("Observed","Shannon"))
}

write.table(animal.normal.rarefy_diversity,"animal.pooled.rarefy.csv",sep=",")
animal.rarefy.table<-as.data.frame(animal.normal.rarefy[[300]]@otu_table)
write.table(animal.rarefy.table,"animal.pooled.rarefy.table.csv",sep=",")
###get the diversity of individual samples
##fungi
SE<-function(x) {
  sd(x)/sqrt(length(x))
}
get.d<-function(fungi.normal.quantile,i){
  fungi.normal.quantile.r<-fungi.normal.quantile[,grepl(i,names(fungi.normal.quantile))]
  row.names(fungi.normal.quantile.r)<-fungi.normal.quantile$site
  fungi.normal.quantile.r.sd<-apply(fungi.normal.quantile.r,1,sd)
  fungi.normal.quantile.r.mean<-apply(fungi.normal.quantile.r,1,mean)
  fungi.normal.quantile.r.SE<-apply(fungi.normal.quantile.r,1,SE)
  fungi.normal.quantile.diversityresult<-data.frame(estimate=fungi.normal.quantile.r.mean,sd=fungi.normal.quantile.r.sd,se=fungi.normal.quantile.r.SE)
  row.names(fungi.normal.quantile.diversityresult)<-row.names(fungi.normal.quantile)
  return(fungi.normal.quantile.diversityresult)
}
fungi.normal.quantile<-read.csv("fungi.rarefy.diversity2.csv" ,header = T,row.names = 1,sep = ",")
fungi.normal.quantile.diversityresult<-get.d(fungi.normal.quantile,"Shannon")
write.table(fungi.normal.quantile.diversityresult,"fungi.shannon.csv",sep=",")
##richness
fungi.normal.quantile.diversityresult<-get.d(fungi.normal.quantile,"Obser")
write.table(fungi.normal.quantile.diversityresult,"fungi.richness.csv",sep=",")
##animal
animal.normal.quantile<-read.csv("animal.rarefy.diversity2.csv" ,header = T,row.names = 1,sep = ",")
animal.normal.quantile.diversityresult<-get.d(animal.normal.quantile,"Shannon")
write.table(animal.normal.quantile.diversityresult,"animal.shannon.csv",sep=",")
##richness
animal.normal.quantile.diversityresult<-get.d(animal.normal.quantile,"Obser")
write.table(animal.normal.quantile.diversityresult,"animal.richness.csv",sep=",")
##bacteria
bacteria.normal.quantile<-read.csv("bacteria.rarefy.diversity.csv" ,header = T,row.names = 1,sep = ",")
bacteria.normal.quantile.diversityresult<-get.d(bacteria.normal.quantile,"Shannon")
write.table(bacteria.normal.quantile.diversityresult,"bacteria.shannon.csv",sep=",")
##richness
bacteria.normal.quantile.diversityresult<-get.d(bacteria.normal.quantile,"Obser")
write.table(bacteria.normal.quantile.diversityresult,"bacteria.richness.csv",sep=",")
###for pooled
#bacteria
bacteria.normal.quantile<-read.csv("bacteria.pooled.rarefy.csv",header = T,sep = ",")
bacteria.normal.quantile.diversityresult<-get.d(bacteria.normal.quantile,"Shannon")
write.table(bacteria.normal.quantile.diversityresult,"bacteria.pooled.shannon.csv",sep=",")
##richness
bacteria.normal.quantile.diversityresult<-get.d(bacteria.normal.quantile,"Obser")
write.table(bacteria.normal.quantile.diversityresult,"bacteria.pooled.richness.csv",sep=",")
#animal
animal.normal.quantile<-read.csv("animal.pooled.rarefy.csv",header = T,sep = ",")
animal.normal.quantile.diversityresult<-get.d(animal.normal.quantile,"Shannon")
write.table(animal.normal.quantile.diversityresult,"animal.pooled.shannon.csv",sep=",")
##richness
animal.normal.quantile.diversityresult<-get.d(animal.normal.quantile,"Obser")
write.table(animal.normal.quantile.diversityresult,"animal.pooled.richness.csv",sep=",")
#fungi
fungi.normal.quantile<-read.csv("fungi.pooled.rarefy.csv",header = T,sep = ",")
fungi.normal.quantile.diversityresult<-get.d(fungi.normal.quantile,"Shannon")
write.table(fungi.normal.quantile.diversityresult,"fungi.pooled.shannon.csv",sep=",")
##richness
fungi.normal.quantile.diversityresult<-get.d(fungi.normal.quantile,"Obser")
write.table(fungi.normal.quantile.diversityresult,"fungi.pooled.richness.csv",sep=",")
###get the table with sampling designs rarefied at summary sequencing depth
###100% table
##bacteria
bacteria.soil<-read.csv("bacteria_soil.csv",header = TRUE,sep = ",",row.names = 1)
names(bacteria.soil)<-gsub("mix","",names(bacteria.soil))
bacteria.soil$Feature.ID<-row.names(bacteria.soil)
bacteria.soil<-bacteria.soil[rowSums(bacteria.soil[,c(1:36)])!=0,]
bacteria.soil<-bacteria.soil[bacteria.soil$Confidence > 0.8,]
names(bacteria.soil)[1:36]<-paste0(names(bacteria.soil)[1:36],"soil")

##DNA
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

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")
metadata2<-rbind(rbind(metadata.DNA,metadata.soil))

###get the unpooled information
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
a<-colSums(bacteria.normal)[1]
##
raref<-data.frame(rarefaction=colSums(bacteria.normal))
raref$site<-row.names(raref)
raref$site2<-substr(row.names(raref),1,2)
raref<-raref[!grepl("soil|DNA",raref$site),]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata<-metadata[,c(1,2)]
raref2<-merge(metadata,raref,by="site",all=T)
raref2<-raref2[!is.na(raref2$rarefaction),]
raref2<-raref2[!is.na(raref2$method),]
raref3<-raref2 %>%
  group_by(site2,method) %>%
  summarise(summarise.rarefa=sum(rarefaction),.groups="keep")
a<-unique(raref3$summarise.rarefa)
##rarefy table
metadata2$method[metadata2$method=="deep"]<-"Deep"
metadata2$method[metadata2$method=="GSMc_62"]<-"GSMc"
rarefy.t<-list()
for (i in 1:length(a)) {
  d1<-rrarefy(t(bacteria2),as.numeric(a[i]))
  try<-paste0(raref3[raref3$summarise.rarefa==a[i],]$method,raref3[raref3$summarise.rarefa==a[i],]$site2)
  try<-metadata2[paste0(metadata2$method,substr(metadata2$site,1,2)) %in% try,]$site
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,names(d1) %in% try]
  d1<-d1[,!grepl("LZ1ADNA|LZ1BDNA|LZ1Asoil|LZ1Bsoil|LV1ADNA|LV1BDNA|LV1Asoil|LV1Bsoil|LW1ADNA|LW1BDNA|LW1Asoil|LW1Bsoil|LVDDNA|LWDDNA|LZDDNA|LVDsoil|LWDsoil|LZDsoil",names(d1))]
  rarefy.t[[i]]<-d1
}
###GSMc40AB
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-c(select[select$X1Amix=="GSMc40A",]$X,
           gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
           gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X))
GSMc40A<-paste0(GSMc40A,collapse = "|")
a<-raref[grepl(GSMc40A,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40A"],collapse="|")
rarefy.t1<-list()
for (i in 1:3) {
  d1<-rrarefy(t(bacteria2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(GSMC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t1[[i]]<-d1
}
###
GSMc40B<-c(select[select$X1Bmix=="GSMc40B",]$X,
           gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
           gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X))
GSMc40B<-paste0(GSMc40B,collapse = "|")
a<-raref[grepl(GSMc40B,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40B"],collapse="|")
rarefy.t2<-list()
for (i in 1:3) {
  d1<-rrarefy(t(bacteria2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(GSMC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t2[[i]]<-d1
}

##deep_SUCC
deep_SUCC1<-paste0(metadata2$site[metadata2$method=="deep_SUCC"],collapse="|")
a<-raref3[grepl("Deep|SUCC",raref3$method),]
a<-a%>%
  group_by(site2) %>%
  summarise(summarise.rarefa2=sum(summarise.rarefa))
rarefy.t3<-list()
for (i in 1:3) {
  d1<-rrarefy(t(bacteria2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(deep_SUCC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t3[[i]]<-d1
}

all<-c(rarefy.t,rarefy.t1,rarefy.t2,rarefy.t3)

for (i in 1:17) {
  try<- all[[i]]
  try$OTU<-row.names(try)
  all[[i]]<-try
}

me<-function(x,y){
  merge(x,y,by="OTU",all=T)
}

all2<-Reduce(me,all)
write.csv(all2,"bacteria.poolat100.csv")
##fungi
###100% table
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

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))

###get the unpooled information
fungi.normal <- read.csv("fungi.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
a<-colSums(fungi.normal)[1]
##
raref<-data.frame(rarefaction=colSums(fungi.normal))
raref$site<-row.names(raref)
raref$site2<-substr(row.names(raref),1,2)
raref<-raref[!grepl("soil|DNA",raref$site),]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata<-metadata[,c(1,2)]
raref2<-merge(metadata,raref,by="site",all=T)
raref2<-raref2[!is.na(raref2$rarefaction),]
raref2<-raref2[!is.na(raref2$method),]
raref3<-raref2 %>%
  group_by(site2,method) %>%
  summarise(summarise.rarefa=sum(rarefaction),.groups="keep")
a<-unique(raref3$summarise.rarefa)
##rarefy table
metadata2$method[metadata2$method=="deep"]<-"Deep"
metadata2$method[metadata2$method=="GSMc_62"]<-"GSMc"
rarefy.t<-list()
for (i in 1:length(a)) {
  d1<-rrarefy(t(fungi2),as.numeric(a[i]))
  try<-paste0(raref3[raref3$summarise.rarefa==a[i],]$method,raref3[raref3$summarise.rarefa==a[i],]$site2)
  try<-metadata2[paste0(metadata2$method,substr(metadata2$site,1,2)) %in% try,]$site
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,names(d1) %in% try]
  d1<-d1[,!grepl("LZ1ADNA|LZ1BDNA|LZ1Asoil|LZ1Bsoil|LV1ADNA|LV1BDNA|LV1Asoil|LV1Bsoil|LW1ADNA|LW1BDNA|LW1Asoil|LW1Bsoil|LVDDNA|LWDDNA|LZDDNA|LVDsoil|LWDsoil|LZDsoil",names(d1))]
  rarefy.t[[i]]<-d1
}
###GSMc40AB
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-c(select[select$X1Amix=="GSMc40A",]$X,
           gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
           gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X))
GSMc40A<-paste0(GSMc40A,collapse = "|")
a<-raref[grepl(GSMc40A,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40A"],collapse="|")
rarefy.t1<-list()
for (i in 1:3) {
  d1<-rrarefy(t(fungi2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(GSMC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t1[[i]]<-d1
}
###
GSMc40B<-c(select[select$X1Bmix=="GSMc40B",]$X,
           gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
           gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X))
GSMc40B<-paste0(GSMc40B,collapse = "|")
a<-raref[grepl(GSMc40B,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40B"],collapse="|")
rarefy.t2<-list()
for (i in 1:3) {
  d1<-rrarefy(t(fungi2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(GSMC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t2[[i]]<-d1
}

##deep_SUCC
deep_SUCC1<-paste0(metadata2$site[metadata2$method=="deep_SUCC"],collapse="|")
a<-raref3[grepl("Deep|SUCC",raref3$method),]
a<-a%>%
  group_by(site2) %>%
  summarise(summarise.rarefa2=sum(summarise.rarefa))
rarefy.t3<-list()
for (i in 1:3) {
  d1<-rrarefy(t(fungi2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(deep_SUCC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t3[[i]]<-d1
}

all<-c(rarefy.t,rarefy.t1,rarefy.t2,rarefy.t3)

for (i in 1:18) {
  try<- all[[i]]
  try$OTU<-row.names(try)
  all[[i]]<-try
}

me<-function(x,y){
  merge(x,y,by="OTU")
}
all2<-Reduce(me,all)
write.csv(all2,"fungi.poolat100.csv")
###100% table
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

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))

###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
a<-colSums(animal.normal)[1]
##
raref<-data.frame(rarefaction=colSums(animal.normal))
raref$site<-row.names(raref)
raref$site2<-substr(row.names(raref),1,2)
raref<-raref[!grepl("soil|DNA",raref$site),]
metadata<-read.csv("normal_metadata.csv",header = TRUE,sep = ",")
metadata<-metadata[,c(1,2)]
raref2<-merge(metadata,raref,by="site",all=T)
raref2<-raref2[!is.na(raref2$rarefaction),]
raref2<-raref2[!is.na(raref2$method),]
raref3<-raref2 %>%
  group_by(site2,method) %>%
  summarise(summarise.rarefa=sum(rarefaction),.groups="keep")
a<-unique(raref3$summarise.rarefa)
##rarefy table
metadata2$method[metadata2$method=="deep"]<-"Deep"
metadata2$method[metadata2$method=="GSMc_62"]<-"GSMc"
rarefy.t<-list()
for (i in 1:length(a)) {
  d1<-rrarefy(t(animal2),as.numeric(a[i]))
  try<-paste0(raref3[raref3$summarise.rarefa==a[i],]$method,raref3[raref3$summarise.rarefa==a[i],]$site2)
  try<-metadata2[paste0(metadata2$method,substr(metadata2$site,1,2)) %in% try,]$site
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,names(d1) %in% try]
  d1<-d1[,!grepl("LZ1ADNA|LZ1BDNA|LZ1Asoil|LZ1Bsoil|LV1ADNA|LV1BDNA|LV1Asoil|LV1Bsoil|LW1ADNA|LW1BDNA|LW1Asoil|LW1Bsoil|LVDDNA|LWDDNA|LZDDNA|LVDsoil|LWDsoil|LZDsoil",names(d1))]
  rarefy.t[[i]]<-d1
}
###GSMc40AB
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-c(select[select$X1Amix=="GSMc40A",]$X,
           gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
           gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X))
GSMc40A<-paste0(GSMc40A,collapse = "|")
a<-raref[grepl(GSMc40A,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40A"],collapse="|")
rarefy.t1<-list()
for (i in 1:3) {
  d1<-rrarefy(t(animal2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(GSMC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t1[[i]]<-d1
}
###
GSMc40B<-c(select[select$X1Bmix=="GSMc40B",]$X,
           gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
           gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X))
GSMc40B<-paste0(GSMc40B,collapse = "|")
a<-raref[grepl(GSMc40B,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40B"],collapse="|")
rarefy.t2<-list()
for (i in 1:3) {
  d1<-rrarefy(t(animal2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(GSMC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t2[[i]]<-d1
}

##deep_SUCC
deep_SUCC1<-paste0(metadata2$site[metadata2$method=="deep_SUCC"],collapse="|")
a<-raref3[grepl("Deep|SUCC",raref3$method),]
a<-a%>%
  group_by(site2) %>%
  summarise(summarise.rarefa2=sum(summarise.rarefa))
rarefy.t3<-list()
for (i in 1:3) {
  d1<-rrarefy(t(animal2),as.numeric(a[i,2]))
  d1<-as.data.frame(d1)
  d1<-as.data.frame(t(d1))
  d1<-d1[,grepl(deep_SUCC1,names(d1))]
  d1<-d1[,grepl(as.character(a[i,1]),names(d1))]
  rarefy.t3[[i]]<-d1
}

all<-c(rarefy.t,rarefy.t1,rarefy.t2,rarefy.t3)

for (i in 1:19) {
  try<- all[[i]]
  try$OTU<-row.names(try)
  all[[i]]<-try
}

me<-function(x,y){
  merge(x,y,by="OTU")
}

all2<-Reduce(me,all)
write.csv(all2,"animal.poolat100.csv")


