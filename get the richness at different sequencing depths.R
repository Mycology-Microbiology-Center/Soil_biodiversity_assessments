##bacteria
library(data.table)
library(tidyr)
library(purrr)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
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
##DarkDiv
DarkDiv1<-paste0(metadata2$site[metadata2$method=="DarkDiv"],collapse="|")
DarkDiv<-bacteria2[,grepl(DarkDiv1,names(bacteria2))]
DarkDiv<-as.data.frame(t(DarkDiv))
p<-rarecurve(DarkDiv, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(DarkDiv),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
DarkDiv<-b
###get the unpooled information
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
a<-colSums(bacteria.normal)[1]
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
DarkDiv1<-bacteria.total.richenss[bacteria.total.richenss$method=="DarkDiv",c(1,3)]
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
###
a<-raref3[grepl("DarkDiv",raref3$method),]
##
DarkDiv1<-merge(DarkDiv1,a[,c(1,3)],by.x="site",by.y="site2")
DarkDiv1$mean<-round(DarkDiv1$mean)
names(DarkDiv1)<-c("Group","value","n_seqs")
DarkDiv1$site<-DarkDiv1$Group
DarkDiv1$type<-paste0(DarkDiv1$site,"unpooled")
DarkDiv1$type1<-"unpooled"
DarkDiv1$method<-"DarkDiv"
all<-rbind(DarkDiv1,b)

##
DarkDiv2<-DarkDiv1[,2:4]
names(DarkDiv2)<-paste0("unpooled",names(DarkDiv2))
DarkDiv50<-merge(b[abs(b$n_seqs-(unique(DarkDiv1$n_seqs))*0.50)<10,],DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv75<-merge(b[abs(b$n_seqs-(unique(DarkDiv1$n_seqs))*0.75)<10,],DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv25<-merge(b[abs(b$n_seqs-(unique(DarkDiv1$n_seqs))*0.25)<10,],DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv2<-merge(b[abs(b$n_seqs-unique(DarkDiv1$n_seqs))<10,],DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv2$proportion<-"100%"
DarkDiv25$proportion<-"25%"
DarkDiv50$proportion<-"50%"
DarkDiv75$proportion<-"75%"
DarkDiv2<-rbind(DarkDiv2,DarkDiv25,DarkDiv50,DarkDiv75)
DarkDiv2$ratio<-DarkDiv2$value/DarkDiv2$unpooledvalue
##GSMC62
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_62"],collapse="|")
GSMC62<-bacteria2[,grepl(GSMC1,names(bacteria2))]
GSMC62<-as.data.frame(t(GSMC62))
p<-rarecurve(GSMC62, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMC62),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMC62<-b
###get the unpooled information
bacteria.normal <- read.csv("bacteria.rarefy.table.csv",header = TRUE,sep = ",",row.names = 1)
a<-colSums(bacteria.normal)[1]
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
GSMC1<-bacteria.total.richenss[bacteria.total.richenss$method=="GSMc",c(1,3)]
###
a<-raref3[grepl("GSMc",raref3$method),]
##
GSMC1<-merge(GSMC1,a[,c(1,3)],by.x="site",by.y="site2")
#
GSMC1$mean<-round(GSMC1$mean)
names(GSMC1)<-c("Group","value","n_seqs")
GSMC1$site<-GSMC1$Group
GSMC1$type<-paste0(GSMC1$site,"unpooled")
GSMC1$type1<-"unpooled"
GSMC1$method<-"GSMC62"
all<-rbind(GSMC1,b)
# Assuming `p` is the plot object returned by plot_predictions()
out<-ggplot(aes(x=n_seqs,y=value,group=type,color=site),data=b[b$n_seqs<75000,])+
  geom_point(data = GSMC1, aes(x = n_seqs, y = value, color = site), size = 3, alpha = 0.8) +  # Enhance points with size and transparency
  geom_line(aes(linetype = type1)) +  # Enhance points with size and transparency
  labs(
    x = "Sequencing depth",
    y = "number of OTUs",
    color = "Sites") +  # Add meaningful labels and titles
  geom_vline(
    xintercept = GSMC1$n_seqs, 
    linetype = "dashed", 
    color = "grey", 
    size = 0.1
  ) +
  geom_hline(
    yintercept = GSMC1$value, 
    linetype = "dashed", 
    color = "grey", 
    size = 0.1
  ) +
  scale_y_continuous(breaks = c(0,500,1000,1500,c(GSMC1$value),(max(b$value)),3000))+
  scale_x_continuous(breaks = c(0,20000,40000,60000,(GSMC1$n_seqs),max(b$n_seqs),70000))+
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

cbPalette <- rev(c("#E69F99","#009E73","#1188D3"))
cbPalette <-setNames(cbPalette,c("LV","LZ","LW"))
out + scale_colour_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette)  
##
GSMC2<-GSMC1[,2:4]
names(GSMC2)<-paste0("unpooled",names(GSMC2))
GSMC25<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")

GSMC50<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")
GSMC75<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")
GSMC2<-merge(rbind(b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")
GSMC25$proportion<-"25%"
GSMC50$proportion<-"50%"
GSMC75$proportion<-"75%"
GSMC2$proportion<-"100%"
GSMC2<-rbind(GSMC2,GSMC25,GSMC50,GSMC75)
GSMC2$ratio<-GSMC2$value/GSMC2$unpooledvalue
###GSM40A
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40A"],collapse="|")
GSMc_40A<-bacteria2[,grepl(GSMC1,names(bacteria2))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
p<-rarecurve(GSMc_40A, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMc_40A),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMc_40A<-b
##
raref<-data.frame(rarefaction=colSums(bacteria.normal))
raref$site<-row.names(raref)
raref$site2<-substr(row.names(raref),1,2)
raref<-raref[!grepl("soil|DNA",raref$site),]
select<-read.csv("sheet2_for_pooled.csv")
GSMc40A<-c(select[select$X1Amix=="GSMc40A",]$X,
           gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
           gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X))
GSMc40A<-paste0(GSMc40A,collapse = "|")
###
a<-raref[grepl(GSMc40A,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
GSMC1<-bacteria.total.richenss[bacteria.total.richenss$method=="GSMc40A",c(1,3)]
GSMC1<-merge(GSMC1,a,by.x="site",by.y="site2")
GSMC1$mean<-round(GSMC1$mean)
names(GSMC1)<-c("Group","value","n_seqs")
GSMC1$site<-GSMC1$Group
GSMC1$type<-paste0(GSMC1$site,"unpooled")
GSMC1$type1<-"unpooled"
all<-rbind(GSMC1,b)

##
GSMC3<-GSMC1[,2:4]
names(GSMC3)<-paste0("unpooled",names(GSMC3))

GSMC25<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")

GSMC50<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")
GSMC75<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")
GSMC3<-merge(rbind(b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")
GSMC25$proportion<-"25%"
GSMC50$proportion<-"50%"
GSMC75$proportion<-"75%"
GSMC3$proportion<-"100%"
GSMC3<-rbind(GSMC3,GSMC25,GSMC50,GSMC75)
GSMC3$ratio<-GSMC3$value/GSMC3$unpooledvalue

##GSMc_40B
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40B"],collapse="|")
GSMc_40B<-bacteria2[,grepl(GSMC1,names(bacteria2))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
p<-rarecurve(GSMc_40B, step = 20,label = F)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMc_40B),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMc_40B<-b
##
GSMc40B<-c(select[select$X1Bmix=="GSMc40B",]$X,
           gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
           gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X))
GSMc40B<-paste0(GSMc40B,collapse = "|")
###
a<-raref[grepl(GSMc40B,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))

###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
GSMC1<-bacteria.total.richenss[bacteria.total.richenss$method=="GSMc40B",c(1,3)]
GSMC1<-merge(GSMC1,a,by.x="site",by.y="site2")
GSMC1$mean<-round(GSMC1$mean)
names(GSMC1)<-c("Group","value","n_seqs")
GSMC1$site<-GSMC1$Group
GSMC1$type<-paste0(GSMC1$site,"unpooled")
GSMC1$type1<-"unpooled"
GSMC1$method<-"GSMc40B"
all<-rbind(GSMC1,b)

###
GSMC4<-GSMC1[,2:4]
names(GSMC4)<-paste0("unpooled",names(GSMC4))

GSMC25<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")

GSMC50<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")
GSMC75<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")
GSMC4<-merge(rbind(b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")
GSMC25$proportion<-"25%"
GSMC50$proportion<-"50%"
GSMC75$proportion<-"75%"
GSMC4$proportion<-"100%"
GSMC4<-rbind(GSMC4,GSMC25,GSMC50,GSMC75)
GSMC4$ratio<-GSMC4$value/GSMC4$unpooledvalue
##LUCAS
LUCAS1<-paste0(metadata2$site[metadata2$method=="LUCAS"],collapse="|")
LUCAS<-bacteria2[,grepl(LUCAS1,names(bacteria2))]
LUCAS<-as.data.frame(t(LUCAS))
p<-rarecurve(LUCAS, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(LUCAS),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
LUCAS<-b
###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
LUCAS1<-bacteria.total.richenss[bacteria.total.richenss$method=="LUCAS",c(1,3)]
#
a<-raref3[grepl("LUCAS",raref3$method),]
##
LUCAS1<-merge(LUCAS1,a[,c(1,3)],by.x="site",by.y="site2")
LUCAS1$mean<-round(LUCAS1$mean)
names(LUCAS1)<-c("Group","value","n_seqs")
LUCAS1$site<-LUCAS1$Group
LUCAS1$type<-paste0(LUCAS1$site,"unpooled")
LUCAS1$type1<-"unpooled"
LUCAS1$method<-"LUCAS"
all<-rbind(LUCAS1,b)
###
LUCAS2<-LUCAS1[,2:4]
names(LUCAS2)<-paste0("unpooled",names(LUCAS2))
LUCAS50<-merge(rbind(b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.50)<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")
LUCAS75<-merge(rbind(b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.75)<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")
LUCAS25<-merge(rbind(b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.25)<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")

LUCAS2<-merge(rbind(b[abs(b$n_seqs-LUCAS1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-LUCAS1$n_seqs[1])<10,]$site),],
                    b[abs(b$n_seqs-LUCAS1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-LUCAS1$n_seqs[2])<10,]$site),],
                    b[abs(b$n_seqs-LUCAS1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-LUCAS1$n_seqs[3])<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")
LUCAS25$proportion<-"25%"
LUCAS50$proportion<-"50%"
LUCAS75$proportion<-"75%"
LUCAS2$proportion<-"100%"
LUCAS2<-rbind(LUCAS2,LUCAS25,LUCAS50,LUCAS75)
LUCAS2$ratio<-LUCAS2$value/LUCAS2$unpooledvalue
##Zobel
Zobel1<-paste0(metadata2$site[metadata2$method=="Zobel"],collapse="|")
Zobel<-bacteria2[,grepl(Zobel1,names(bacteria2))]
Zobel<-as.data.frame(t(Zobel))
p<-rarecurve(Zobel, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(Zobel),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
Zobel<-b
###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
Zobel1<-bacteria.total.richenss[bacteria.total.richenss$method=="Zobel",c(1,3)]
a<-raref3[grepl("Zobel",raref3$method),]
Zobel1<-merge(Zobel1,a[,c(1,3)],by.x="site",by.y="site2")
Zobel1$mean<-round(Zobel1$mean)
names(Zobel1)<-c("Group","value","n_seqs")
Zobel1$site<-Zobel1$Group
Zobel1$type<-paste0(Zobel1$site,"unpooled")
Zobel1$type1<-"unpooled"
Zobel1$method<-"Zobel"
all<-rbind(Zobel1,b)
#
Zobel2<-Zobel1[,2:4]
names(Zobel2)<-paste0("unpooled",names(Zobel2))
Zobel50<-merge(rbind(b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.50)<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")
Zobel75<-merge(rbind(b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.75)<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")
Zobel25<-merge(rbind(b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.25)<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")

Zobel2<-merge(rbind(b[abs(b$n_seqs-Zobel1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-Zobel1$n_seqs[1])<10,]$site),],
                    b[abs(b$n_seqs-Zobel1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-Zobel1$n_seqs[2])<10,]$site),],
                    b[abs(b$n_seqs-Zobel1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-Zobel1$n_seqs[3])<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")
Zobel25$proportion<-"25%"
Zobel50$proportion<-"50%"
Zobel75$proportion<-"75%"
Zobel2$proportion<-"100%"
Zobel2<-rbind(Zobel2,Zobel25,Zobel50,Zobel75)
Zobel2$ratio<-Zobel2$value/Zobel2$unpooledvalue

##deep
deep1<-paste0(metadata2$site[metadata2$method=="deep"],collapse="|")
deep<-bacteria2[,grepl(deep1,names(bacteria2))]
deep<-as.data.frame(t(deep))
p<-rarecurve(deep, step = 20,label = F)

library(tidyr)
library(purrr)
library(stringr)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(deep),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))

deep<-b
###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
deep1<-bacteria.total.richenss[bacteria.total.richenss$method=="deep",c(1,3)]
a<-raref3[grepl("Deep",raref3$method),]
deep1<-merge(deep1,a[,c(1,3)],by.x="site",by.y="site2")
deep1$mean<-round(deep1$mean)
names(deep1)<-c("Group","value","n_seqs")
deep1$site<-deep1$Group
deep1$type<-paste0(deep1$site,"unpooled")
deep1$type1<-"unpooled"
all<-rbind(deep1,b)

##
deep2<-deep1[,2:4]
names(deep2)<-paste0("unpooled",names(deep2))
deep50<-merge(rbind(b[abs(b$n_seqs-(deep1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(deep1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(deep1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(deep1$n_seqs[3])*0.50)<10,]$site),]),deep2,by.x="site",by.y="unpooledsite")
deep75<-merge(rbind(b[abs(b$n_seqs-(deep1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(deep1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(deep1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(deep1$n_seqs[3])*0.75)<10,]$site),]),deep2,by.x="site",by.y="unpooledsite")
deep25<-merge(rbind(b[abs(b$n_seqs-(deep1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(deep1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(deep1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(deep1$n_seqs[3])*0.25)<10,]$site),]),deep2,by.x="site",by.y="unpooledsite")
deep2<-merge((rbind(b[abs(b$n_seqs-deep1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-deep1$n_seqs[1])<10,]$site),],
                    b[abs(b$n_seqs-deep1$n_seqs[2])<12,][grepl("LW",b[abs(b$n_seqs-deep1$n_seqs[2])<12,]$site),],
                    b[abs(b$n_seqs-deep1$n_seqs[3])<12,][grepl("LZ",b[abs(b$n_seqs-deep1$n_seqs[3])<12,]$site),]))[c(1:3,5,7,9),],deep2,by.x="site",by.y="unpooledsite")

deep25$proportion<-"25%"
deep50$proportion<-"50%"
deep75$proportion<-"75%"
deep2$proportion<-"100%"
deep2<-rbind(deep2,deep25,deep50,deep75)
deep2$ratio<-deep2$value/deep2$unpooledvalue

##deep_SUCC
deep_SUCC1<-paste0(metadata2$site[metadata2$method=="deep_SUCC"],collapse="|")
deep_SUCC<-bacteria2[,grepl(deep_SUCC1,names(bacteria2))]
deep_SUCC<-as.data.frame(t(deep_SUCC))
p<-rarecurve(deep_SUCC, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(deep_SUCC),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
deep_SUCC<-b
###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
deep_SUCC1<-bacteria.total.richenss[bacteria.total.richenss$method=="deep_SUCC",c(1,3)]
##
a<-raref3[grepl("Deep|SUCC",raref3$method),]
a<-a%>%
  group_by(site2) %>%
  summarise(summarise.rarefa2=sum(summarise.rarefa))
deep_SUCC1<-merge(deep_SUCC1,a,by.x="site",by.y="site2")
##
deep_SUCC1$mean<-round(deep_SUCC1$mean)
names(deep_SUCC1)<-c("Group","value","n_seqs")
deep_SUCC1$site<-deep_SUCC1$Group
deep_SUCC1$type<-paste0(deep_SUCC1$site,"unpooled")
deep_SUCC1$type1<-"unpooled"
all<-rbind(deep_SUCC1,b)
##
deep_SUCC2<-deep_SUCC1[,2:4]
names(deep_SUCC2)<-paste0("unpooled",names(deep_SUCC2))
deep_SUCC50<-merge(rbind(b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.50)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.50)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[3])*0.50)<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC75<-merge(rbind(b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.75)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.75)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[3])*0.75)<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC25<-merge(rbind(b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.25)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.25)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[3])*0.25)<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC2<-merge(rbind(b[abs(b$n_seqs-deep_SUCC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-deep_SUCC1$n_seqs[1])<10,]$site),],
                        b[abs(b$n_seqs-deep_SUCC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-deep_SUCC1$n_seqs[2])<10,]$site),],
                        b[abs(b$n_seqs-deep_SUCC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-deep_SUCC1$n_seqs[3])<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC25$proportion<-"25%"
deep_SUCC50$proportion<-"50%"
deep_SUCC75$proportion<-"75%"
deep_SUCC2$proportion<-"100%"
deep_SUCC2<-rbind(deep_SUCC2,deep_SUCC25,deep_SUCC50,deep_SUCC75)
deep_SUCC2$ratio<-deep_SUCC2$value/deep_SUCC2$unpooledvalue  
##SUCC
SUCC1<-paste0(metadata2$site[metadata2$method=="SUCC"],collapse="|")
SUCC<-bacteria2[,grepl(SUCC1,names(bacteria2))]
SUCC<-as.data.frame(t(SUCC))
p<-rarecurve(SUCC, step = 20,label = F)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(SUCC),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
SUCC<-b
##
###get the unpooled information
bacteria.total.richenss <- read.csv("bacteria.totalrichness.csv",header = TRUE,sep = ",",row.names = 1)
SUCC1<-bacteria.total.richenss[bacteria.total.richenss$method=="SUCC",c(1,3)]
a<-raref3[grepl("SUCC",raref3$method),]
SUCC1<-merge(SUCC1,a[,c(1,3)],by.x="site",by.y="site2")
SUCC1$mean<-round(SUCC1$mean)
names(SUCC1)<-c("Group","value","n_seqs")
SUCC1$site<-SUCC1$Group
SUCC1$type<-paste0(SUCC1$site,"unpooled")
SUCC1$type1<-"unpooled"
all<-rbind(SUCC1,b)

##
SUCC2<-SUCC1[,2:4]
names(SUCC2)<-paste0("unpooled",names(SUCC2))
SUCC50<-merge(rbind(b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.50)<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC75<-merge(rbind(b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.75)<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC25<-merge(rbind(b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.25)<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC2<-merge(rbind(b[abs(b$n_seqs-SUCC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-SUCC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-SUCC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-SUCC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-SUCC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-SUCC1$n_seqs[3])<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC25$proportion<-"25%"
SUCC50$proportion<-"50%"
SUCC75$proportion<-"75%"
SUCC2$proportion<-"100%"
SUCC2<-rbind(SUCC2,SUCC25,SUCC50,SUCC75)
SUCC2$ratio<-SUCC2$value/SUCC2$unpooledvalue  
##
deep2$method<-"deep"
LUCAS2$method<-"LUCAS"
DarkDiv2$method<-"DarkDiv"
Zobel2$method<-"Zobel"
GSMC2$method<-"GSMc_62"
GSMC3$method<-"GSMc_40A"
GSMC4$method<-"GSMc_40B"
deep_SUCC2$method<-"deep_SUCC"
SUCC2$method<-"SUCC"
all<-rbind(SUCC2,DarkDiv2,Zobel2,GSMC2,GSMC3,GSMC4,deep2,deep_SUCC2,LUCAS2)
write.csv(all,"proportion.bacteria.richness.summarydepthforpooled.csv")
##
SUCC$method<-"SUCC"
DarkDiv$method<-"DarkDiv"
Zobel$method<-"Zobel"
GSMC62$method<-"GSMC62"
GSMc_40A$method<-"GSMc40A"
GSMc_40B$method<-"GSMc40B"
deep$method<-"deep"
deep_SUCC$method<-"deep_SUCC"
LUCAS$method<-"LUCAS"
all1<-rbind(SUCC,DarkDiv,Zobel,GSMC62,GSMc_40A,GSMc_40B,deep,deep_SUCC,LUCAS)
all1<-all1[all1$n_seqs>39290,]
##
all2<-all %>%
  group_by(type1,method) %>%
  summarize(mean.ratio=mean(ratio),sd=sd(ratio),.groups = "keep")

all3<-all %>%
  group_by(type1) %>%
  summarize(mean.ratio=mean(ratio),sd=sd(ratio))
##
library(marginaleffects)
library(lme4)
library(lmerTest)
library(ggplot2)
library(rcompanion)
library(dplyr)
all2<-all[,c(1,6,9,10)]
##DNA
DNA<-all2[grepl("DNA",all2$type1),]
mod01<-lmer(ratio ~ method+ (1|site),data = DNA)
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
  labs(x="method",y="Pooling effect")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##soil
soil<-all2[grepl("soil",all2$type1),]
mod01<-lmer(ratio ~ method+ (1|site),data = soil)
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
  labs(x="method",y="Pooling effect")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  scale_y_continuous(breaks = c(1,1.4,1.8,2.2))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))  
##animal
library(data.table)
setwd("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/animal_table")
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
library(data.table)
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
##DarkDiv
DarkDiv1<-paste0(metadata2$site[metadata2$method=="DarkDiv"],collapse="|")
DarkDiv<-animal2[,grepl(DarkDiv1,names(animal2))]
DarkDiv<-as.data.frame(t(DarkDiv))
p<-rarecurve(DarkDiv, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(DarkDiv),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
DarkDiv<-b
#
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
DarkDiv1<-animal.total.richenss[animal.total.richenss$method=="DarkDiv",c(1,3)]
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
###
a<-raref3[grepl("DarkDiv",raref3$method),]
##
DarkDiv1<-merge(DarkDiv1,a[,c(1,3)],by.x="site",by.y="site2")

###get the unpooled information
DarkDiv1$mean<-round(DarkDiv1$mean)
names(DarkDiv1)<-c("Group","value","n_seqs")
DarkDiv1$site<-DarkDiv1$Group
DarkDiv1$type<-paste0(DarkDiv1$site,"unpooled")
DarkDiv1$type1<-"unpooled"
all<-rbind(DarkDiv1,b)

###
DarkDiv2<-DarkDiv1[,2:4]
names(DarkDiv2)<-paste0("unpooled",names(DarkDiv2))
DarkDiv25<-merge(rbind(b[abs(b$n_seqs-(DarkDiv1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(DarkDiv1$n_seqs[1])*0.25)<10,]$site),],
                       b[abs(b$n_seqs-(DarkDiv1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(DarkDiv1$n_seqs[2])*0.25)<10,]$site),],
                       b[abs(b$n_seqs-(DarkDiv1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(DarkDiv1$n_seqs[3])*0.25)<10,]$site),]),DarkDiv2,by.x="site",by.y="unpooledsite")

DarkDiv50<-merge(rbind(b[abs(b$n_seqs-(DarkDiv1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(DarkDiv1$n_seqs[1])*0.50)<10,]$site),],
                       b[abs(b$n_seqs-(DarkDiv1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(DarkDiv1$n_seqs[2])*0.50)<10,]$site),],
                       b[abs(b$n_seqs-(DarkDiv1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(DarkDiv1$n_seqs[3])*0.50)<10,]$site),]),DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv75<-merge(rbind(b[abs(b$n_seqs-(DarkDiv1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(DarkDiv1$n_seqs[1])*0.75)<10,]$site),],
                       b[abs(b$n_seqs-(DarkDiv1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(DarkDiv1$n_seqs[2])*0.75)<10,]$site),],
                       b[abs(b$n_seqs-(DarkDiv1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(DarkDiv1$n_seqs[3])*0.75)<10,]$site),]),DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv2<-merge(rbind(b[abs(b$n_seqs-DarkDiv1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-DarkDiv1$n_seqs[1])<10,]$site),],
                      b[abs(b$n_seqs-DarkDiv1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-DarkDiv1$n_seqs[2])<10,]$site),],
                      b[abs(b$n_seqs-DarkDiv1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-DarkDiv1$n_seqs[3])<10,]$site),]),DarkDiv2,by.x="site",by.y="unpooledsite")
DarkDiv2$proportion<-"100%"
DarkDiv25$proportion<-"25%"
DarkDiv50$proportion<-"50%"
DarkDiv75$proportion<-"75%"
DarkDiv2<-rbind(DarkDiv2,DarkDiv25,DarkDiv50,DarkDiv75)
DarkDiv2$ratio<-DarkDiv2$value/DarkDiv2$unpooledvalue
##GSMC62
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_62"],collapse="|")
GSMC62<-animal2[,grepl(GSMC1,names(animal2))]
GSMC62<-as.data.frame(t(GSMC62))
p<-rarecurve(GSMC62, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMC62),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMC62<-b
###get the unpooled information
animal.normal <- read.csv("animal.rarefy.table2.csv",header = TRUE,sep = ",",row.names = 1)
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
GSMC1<-animal.total.richenss[animal.total.richenss$method=="GSMc",c(1,3)]
a<-raref3[grepl("GSMc",raref3$method),]
##
GSMC1<-merge(GSMC1,a[,c(1,3)],by.x="site",by.y="site2")

GSMC1$mean<-round(GSMC1$mean)
names(GSMC1)<-c("Group","value","n_seqs")
GSMC1$site<-GSMC1$Group
GSMC1$type<-paste0(GSMC1$site,"unpooled")
GSMC1$type1<-"unpooled"
all<-rbind(GSMC1,b)
##
GSMC2<-GSMC1[,2:4]
names(GSMC2)<-paste0("unpooled",names(GSMC2))
GSMC25<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")

GSMC50<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")
GSMC75<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")
GSMC2<-merge(rbind(b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,]$site),]),GSMC2,by.x="site",by.y="unpooledsite")
GSMC25$proportion<-"25%"
GSMC50$proportion<-"50%"
GSMC75$proportion<-"75%"
GSMC2$proportion<-"100%"
GSMC2<-rbind(GSMC2,GSMC25,GSMC50,GSMC75)

GSMC2$ratio<-GSMC2$value/GSMC2$unpooledvalue
###GSM40A
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40A"],collapse="|")
GSMc_40A<-animal2[,grepl(GSMC1,names(animal2))]
GSMc_40A<-as.data.frame(t(GSMc_40A))
p<-rarecurve(GSMc_40A, step = 20,label = F)

library(tidyr)
library(purrr)
library(stringr)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMc_40A),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMc_40A<-b
##
raref<-data.frame(rarefaction=colSums(animal.normal))
raref$site<-row.names(raref)
raref$site2<-substr(row.names(raref),1,2)
raref<-raref[!grepl("soil|DNA",raref$site),]
select<-read.csv("C:/Users/meirong/Desktop/PhD project/downstream/richness and shannon/forth/pooled and unpooled comparison/1.accumulative/sheet2_for_pooled.csv")
GSMc40A<-c(select[select$X1Amix=="GSMc40A",]$X,
           gsub("LZ","LV",select[select$X1Amix=="GSMc40A",]$X),
           gsub("LZ","LW",select[select$X1Amix=="GSMc40A",]$X))
GSMc40A<-paste0(GSMc40A,collapse = "|")
###
a<-raref[grepl(GSMc40A,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
GSMC1<-animal.total.richenss[animal.total.richenss$method=="GSMc40A",c(1,3)]
GSMC1<-merge(GSMC1,a,by.x="site",by.y="site2")
GSMC1$mean<-round(GSMC1$mean)
names(GSMC1)<-c("Group","value","n_seqs")
GSMC1$site<-GSMC1$Group
GSMC1$type<-paste0(GSMC1$site,"unpooled")
GSMC1$type1<-"unpooled"
all<-rbind(GSMC1,b)

##
GSMC3<-GSMC1[,2:4]
names(GSMC3)<-paste0("unpooled",names(GSMC3))
GSMC25<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")

GSMC50<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")
GSMC75<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")
GSMC3<-merge(rbind(b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,]$site),]),GSMC3,by.x="site",by.y="unpooledsite")
GSMC25$proportion<-"25%"
GSMC50$proportion<-"50%"
GSMC75$proportion<-"75%"
GSMC3$proportion<-"100%"
GSMC3<-rbind(GSMC3,GSMC25,GSMC50,GSMC75)

GSMC3$ratio<-GSMC3$value/GSMC3$unpooledvalue

##GSMc_40B
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40B"],collapse="|")
GSMc_40B<-animal2[,grepl(GSMC1,names(animal2))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
p<-rarecurve(GSMc_40B, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMc_40B),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMc_40B<-b
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
GSMC1<-animal.total.richenss[animal.total.richenss$method=="GSMc40B",c(1,3)]
##
GSMc40B<-c(select[select$X1Bmix=="GSMc40B",]$X,
           gsub("LZ","LV",select[select$X1Bmix=="GSMc40B",]$X),
           gsub("LZ","LW",select[select$X1Bmix=="GSMc40B",]$X))
GSMc40B<-paste0(GSMc40B,collapse = "|")
###
a<-raref[grepl(GSMc40B,raref$site),]
a<-a %>%
  group_by(site2) %>%
  summarise(summarise.rarefaction=sum(rarefaction))
##
GSMC1<-merge(GSMC1,a,by.x="site",by.y="site2")
GSMC1$mean<-round(GSMC1$mean)
names(GSMC1)<-c("Group","value","n_seqs")
GSMC1$site<-GSMC1$Group
GSMC1$type<-paste0(GSMC1$site,"unpooled")
GSMC1$type1<-"unpooled"
all<-rbind(GSMC1,b)

##
GSMC4<-GSMC1[,2:4]
names(GSMC4)<-paste0("unpooled",names(GSMC4))
GSMC25<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.25)<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")

GSMC50<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.50)<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")
GSMC75<-merge(rbind(b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(GSMC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(GSMC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(GSMC1$n_seqs[3])*0.75)<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")
GSMC4<-merge(rbind(b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-GSMC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-GSMC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-GSMC1$n_seqs[3])<10,]$site),]),GSMC4,by.x="site",by.y="unpooledsite")
GSMC25$proportion<-"25%"
GSMC50$proportion<-"50%"
GSMC75$proportion<-"75%"
GSMC4$proportion<-"100%"
GSMC4<-rbind(GSMC4,GSMC25,GSMC50,GSMC75)

GSMC4$ratio<-GSMC4$value/GSMC4$unpooledvalue
##LUCAS
LUCAS1<-paste0(metadata2$site[metadata2$method=="LUCAS"],collapse="|")
LUCAS<-animal2[,grepl(LUCAS1,names(animal2))]
LUCAS<-as.data.frame(t(LUCAS))
p<-rarecurve(LUCAS, step = 20,label = F)

library(tidyr)
library(purrr)
library(stringr)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(LUCAS),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
LUCAS<-b
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
LUCAS1<-animal.total.richenss[animal.total.richenss$method=="LUCAS",c(1,3)]
a<-raref3[grepl("LUCAS",raref3$method),]
##
LUCAS1<-merge(LUCAS1,a[,c(1,3)],by.x="site",by.y="site2")
LUCAS1$mean<-round(LUCAS1$mean)
names(LUCAS1)<-c("Group","value","n_seqs")
LUCAS1$site<-LUCAS1$Group
LUCAS1$type<-paste0(LUCAS1$site,"unpooled")
LUCAS1$type1<-"unpooled"
all<-rbind(LUCAS1,b)

##
LUCAS2<-LUCAS1[,2:4]
names(LUCAS2)<-paste0("unpooled",names(LUCAS2))

LUCAS50<-merge(rbind(b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.50)<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")
LUCAS75<-merge(rbind(b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.75)<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")
LUCAS25<-merge(rbind(b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(LUCAS1$n_seqs[1])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(LUCAS1$n_seqs[2])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(LUCAS1$n_seqs[3])*0.25)<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")

LUCAS2<-merge(rbind(b[abs(b$n_seqs-LUCAS1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-LUCAS1$n_seqs[1])<10,]$site),],
                    b[abs(b$n_seqs-LUCAS1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-LUCAS1$n_seqs[2])<10,]$site),],
                    b[abs(b$n_seqs-LUCAS1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-LUCAS1$n_seqs[3])<10,]$site),]),LUCAS2,by.x="site",by.y="unpooledsite")
LUCAS25$proportion<-"25%"
LUCAS50$proportion<-"50%"
LUCAS75$proportion<-"75%"
LUCAS2$proportion<-"100%"
LUCAS2<-rbind(LUCAS2,LUCAS25,LUCAS50,LUCAS75)
LUCAS2$ratio<-LUCAS2$value/LUCAS2$unpooledvalue
LUCAS2$method<-"LUCAS"
##Zobel
Zobel1<-paste0(metadata2$site[metadata2$method=="Zobel"],collapse="|")
Zobel<-animal2[,grepl(Zobel1,names(animal2))]
Zobel<-as.data.frame(t(Zobel))
p<-rarecurve(Zobel, step = 20,label = F)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(Zobel),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
Zobel<-b
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
Zobel1<-animal.total.richenss[animal.total.richenss$method=="Zobel",c(1,3)]
a<-raref3[grepl("Zobel",raref3$method),]
Zobel1<-merge(Zobel1,a[,c(1,3)],by.x="site",by.y="site2")
Zobel1$mean<-round(Zobel1$mean)
names(Zobel1)<-c("Group","value","n_seqs")
Zobel1$site<-Zobel1$Group
Zobel1$type<-paste0(Zobel1$site,"unpooled")
Zobel1$type1<-"unpooled"
all<-rbind(Zobel1,b)

##
Zobel2<-Zobel1[,2:4]
names(Zobel2)<-paste0("unpooled",names(Zobel2))

Zobel50<-merge(rbind(b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.50)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.50)<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")
Zobel75<-merge(rbind(b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.75)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.75)<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")
Zobel25<-merge(rbind(b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(Zobel1$n_seqs[1])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(Zobel1$n_seqs[2])*0.25)<10,]$site),],
                     b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(Zobel1$n_seqs[3])*0.25)<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")

Zobel2<-merge(rbind(b[abs(b$n_seqs-Zobel1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-Zobel1$n_seqs[1])<10,]$site),],
                    b[abs(b$n_seqs-Zobel1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-Zobel1$n_seqs[2])<10,]$site),],
                    b[abs(b$n_seqs-Zobel1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-Zobel1$n_seqs[3])<10,]$site),]),Zobel2,by.x="site",by.y="unpooledsite")
Zobel25$proportion<-"25%"
Zobel50$proportion<-"50%"
Zobel75$proportion<-"75%"
Zobel2$proportion<-"100%"
Zobel2<-rbind(Zobel2,Zobel25,Zobel50,Zobel75)
Zobel2$ratio<-Zobel2$value/Zobel2$unpooledvalue

##deep
deep1<-paste0(metadata2$site[metadata2$method=="deep"],collapse="|")
deep<-animal2[,grepl(deep1,names(animal2))]
deep<-as.data.frame(t(deep))
p<-rarecurve(deep, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(deep),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
deep<-b
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
deep1<-animal.total.richenss[animal.total.richenss$method=="deep",c(1,3)]
a<-raref3[grepl("Deep",raref3$method),]
deep1<-merge(deep1,a[,c(1,3)],by.x="site",by.y="site2")
deep1$mean<-round(deep1$mean)
names(deep1)<-c("Group","value","n_seqs")
deep1$site<-deep1$Group
deep1$type<-paste0(deep1$site,"unpooled")
deep1$type1<-"unpooled"
all<-rbind(deep1,b)

##
deep2<-deep1[,2:4]
names(deep2)<-paste0("unpooled",names(deep2))
deep50<-merge(rbind(b[abs(b$n_seqs-(deep1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(deep1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[2])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(deep1$n_seqs[2])*0.50)<10,]$site),]),deep2,by.x="site",by.y="unpooledsite")
deep75<-merge(rbind(b[abs(b$n_seqs-(deep1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(deep1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[2])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(deep1$n_seqs[2])*0.75)<10,]$site),]),deep2,by.x="site",by.y="unpooledsite")
deep25<-merge(rbind(b[abs(b$n_seqs-(deep1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(deep1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(deep1$n_seqs[2])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(deep1$n_seqs[2])*0.25)<10,]$site),]),deep2,by.x="site",by.y="unpooledsite")


deep2<-merge((rbind(b[grepl("LV",b$site),][grepl("soil",b[grepl("LV",b$site),]$type1),][which.max(b[grepl("LV",b$site),][grepl("soil",b[grepl("LV",b$site),]$type1),]$n_seqs),],
                    b[grepl("LV",b$site),][grepl("DNA",b[grepl("LV",b$site),]$type1),][which.max(b[grepl("LV",b$site),][grepl("DNA",b[grepl("LV",b$site),]$type1),]$n_seqs),],
                    b[abs(b$n_seqs-deep1$n_seqs[2])<10,][grepl("LZ",b[abs(b$n_seqs-deep1$n_seqs[2])<10,]$site),])),deep2,by.x="site",by.y="unpooledsite")


deep25$proportion<-"25%"
deep50$proportion<-"50%"
deep75$proportion<-"75%"
deep2$proportion<-"100%"
deep2<-rbind(deep2,deep25,deep50,deep75)
deep2$ratio<-deep2$value/deep2$unpooledvalue
##deep_SUCC
deep_SUCC1<-paste0(metadata2$site[metadata2$method=="deep_SUCC"],collapse="|")
deep_SUCC<-animal2[,grepl(deep_SUCC1,names(animal2))]
deep_SUCC<-as.data.frame(t(deep_SUCC))
p<-rarecurve(deep_SUCC, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(deep_SUCC),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
deep_SUCC<-b
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
deep_SUCC1<-animal.total.richenss[animal.total.richenss$method=="deep_SUCC",c(1,3)]
##
a<-raref3[grepl("Deep|SUCC",raref3$method),]
a<-a%>%
  group_by(site2) %>%
  summarise(summarise.rarefa2=sum(summarise.rarefa))
deep_SUCC1<-merge(deep_SUCC1,a,by.x="site",by.y="site2")

deep_SUCC1$mean<-round(deep_SUCC1$mean)
names(deep_SUCC1)<-c("Group","value","n_seqs")
deep_SUCC1$site<-deep_SUCC1$Group
deep_SUCC1$type<-paste0(deep_SUCC1$site,"unpooled")
deep_SUCC1$type1<-"unpooled"
all<-rbind(deep_SUCC1,b)

##
deep_SUCC2<-deep_SUCC1[,2:4]
names(deep_SUCC2)<-paste0("unpooled",names(deep_SUCC2))
deep_SUCC50<-merge(rbind(b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.50)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.50)<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC75<-merge(rbind(b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.75)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.75)<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC25<-merge(rbind(b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[1])*0.25)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.25)<10,]$site),],
                         b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(deep_SUCC1$n_seqs[2])*0.25)<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC2<-merge(rbind(b[abs(b$n_seqs-deep_SUCC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-deep_SUCC1$n_seqs[1])<10,]$site),],
                        b[abs(b$n_seqs-deep_SUCC1$n_seqs[2])<10,][grepl("LZ",b[abs(b$n_seqs-deep_SUCC1$n_seqs[2])<10,]$site),]),deep_SUCC2,by.x="site",by.y="unpooledsite")
deep_SUCC25$proportion<-"25%"
deep_SUCC50$proportion<-"50%"
deep_SUCC75$proportion<-"75%"
deep_SUCC2$proportion<-"100%"
deep_SUCC2<-rbind(deep_SUCC2,deep_SUCC25,deep_SUCC50,deep_SUCC75)

deep_SUCC2$ratio<-deep_SUCC2$value/deep_SUCC2$unpooledvalue  

##SUCC
SUCC1<-paste0(metadata2$site[metadata2$method=="SUCC"],collapse="|")
SUCC<-animal2[,grepl(SUCC1,names(animal2))]
SUCC<-as.data.frame(t(SUCC))
p<-rarecurve(SUCC, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(SUCC),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
SUCC<-b
###get the unpooled information
animal.total.richenss <- read.csv("animal.totalrichness.csv",header = TRUE,sep = ",")
SUCC1<-animal.total.richenss[animal.total.richenss$method=="SUCC",c(1,3)]
a<-raref3[grepl("SUCC",raref3$method),]
SUCC1<-merge(SUCC1,a[,c(1,3)],by.x="site",by.y="site2")
SUCC1$mean<-round(SUCC1$mean)
names(SUCC1)<-c("Group","value","n_seqs")
SUCC1$site<-SUCC1$Group
SUCC1$type<-paste0(SUCC1$site,"unpooled")
SUCC1$type1<-"unpooled"
all<-rbind(SUCC1,b)

SUCC2<-SUCC1[,2:4]
names(SUCC2)<-paste0("unpooled",names(SUCC2))

SUCC50<-merge(rbind(b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.50)<10,][grepl("LV",b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.50)<10,][grepl("LW",b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.50)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.50)<10,][grepl("LZ",b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.50)<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC75<-merge(rbind(b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.75)<10,][grepl("LV",b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.75)<10,][grepl("LW",b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.75)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.75)<10,][grepl("LZ",b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.75)<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC25<-merge(rbind(b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.25)<10,][grepl("LV",b[abs(b$n_seqs-(SUCC1$n_seqs[1])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.25)<10,][grepl("LW",b[abs(b$n_seqs-(SUCC1$n_seqs[2])*0.25)<10,]$site),],
                    b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.25)<10,][grepl("LZ",b[abs(b$n_seqs-(SUCC1$n_seqs[3])*0.25)<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC2<-merge(rbind(b[abs(b$n_seqs-SUCC1$n_seqs[1])<10,][grepl("LV",b[abs(b$n_seqs-SUCC1$n_seqs[1])<10,]$site),],
                   b[abs(b$n_seqs-SUCC1$n_seqs[2])<10,][grepl("LW",b[abs(b$n_seqs-SUCC1$n_seqs[2])<10,]$site),],
                   b[abs(b$n_seqs-SUCC1$n_seqs[3])<10,][grepl("LZ",b[abs(b$n_seqs-SUCC1$n_seqs[3])<10,]$site),]),SUCC2,by.x="site",by.y="unpooledsite")
SUCC25$proportion<-"25%"
SUCC50$proportion<-"50%"
SUCC75$proportion<-"75%"
SUCC2$proportion<-"100%"
SUCC2<-rbind(SUCC2,SUCC25,SUCC50,SUCC75)
SUCC2$ratio<-SUCC2$value/SUCC2$unpooledvalue  
##
deep2$method<-"deep"
LUCAS2$method<-"LUCAS"
DarkDiv2$method<-"DarkDiv"
Zobel2$method<-"Zobel"
GSMC2$method<-"GSMc_62"
GSMC3$method<-"GSMc_40A"
GSMC4$method<-"GSMc_40B"
deep_SUCC2$method<-"deep_SUCC"
SUCC2$method<-"SUCC"
all<-rbind(SUCC2,DarkDiv2,Zobel2,GSMC2,GSMC3,GSMC4,deep2,deep_SUCC2,LUCAS2)
write.csv(all,"proportion.animal.richness.summarydepthforpooled.csv")
##
SUCC$method<-"SUCC"
DarkDiv$method<-"DarkDiv"
Zobel$method<-"Zobel"
GSMC62$method<-"GSMC62"
GSMc_40A$method<-"GSMc40B"
GSMc_40B$method<-"GSMc40B"
deep$method<-"deep"
deep_SUCC$method<-"deep_SUCC"
LUCAS$method<-"LUCAS"
all1<-rbind(SUCC,DarkDiv,Zobel,GSMC62,GSMc_40A,GSMc_40B,deep,deep_SUCC,LUCAS)
all1<-all1[all1$n_seqs>15962,]


##
all2<-all %>%
  group_by(type1,method) %>%
  summarize(mean.ratio=mean(ratio),sd=sd(ratio),.groups = "keep")

all3<-all %>%
  group_by(type1) %>%
  summarize(mean.ratio=mean(ratio),sd=sd(ratio))
all2<-all[,c(1,6,9,10)]
##DNA
DNA<-all2[grepl("DNA",all2$type1),]
mod01<-lmer(ratio ~ method+ (1|site),data = DNA)
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
  labs(x="method",y="Pooling effect")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##soil
soil<-all2[grepl("soil",all2$type1),]
mod01<-lmer(ratio ~ method+ (1|site),data = soil)
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
  labs(x="method",y="Pooling effect")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))  
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

metadata<-read.csv("soil.metadata.csv",header = TRUE,sep = ",")
metadata.DNA<-metadata
metadata.DNA$site<-paste0(metadata.DNA$site,"DNA")

metadata.soil<-metadata
metadata.soil$site<-paste0(metadata.soil$site,"soil")

metadata2<-rbind(rbind(metadata.DNA,metadata.soil))

##DarkDiv
DarkDiv1<-paste0(metadata2$site[metadata2$method=="DarkDiv"],collapse="|")
DarkDiv<-fungi2[,grepl(DarkDiv1,names(fungi2))]
DarkDiv<-as.data.frame(t(DarkDiv))
p<-rarecurve(DarkDiv, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(DarkDiv),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
DarkDiv<-b

##GSMC62
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_62"],collapse="|")
GSMC62<-fungi2[,grepl(GSMC1,names(fungi2))]
GSMC62<-as.data.frame(t(GSMC62))
p<-rarecurve(GSMC62, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMC62),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMC62<-b

##GSMc_40B
GSMC1<-paste0(metadata2$site[metadata2$method=="GSMc_40B"],collapse="|")
GSMc_40B<-fungi2[,grepl(GSMC1,names(fungi2))]
GSMc_40B<-as.data.frame(t(GSMc_40B))
p<-rarecurve(GSMc_40B, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(GSMc_40B),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,5,nchar(b$Group)))
b$type1<-substr(b$Group,5,nchar(b$Group))
GSMc_40B<-b

##LUCAS
LUCAS1<-paste0(metadata2$site[metadata2$method=="LUCAS"],collapse="|")
LUCAS<-fungi2[,grepl(LUCAS1,names(fungi2))]
LUCAS<-as.data.frame(t(LUCAS))
p<-rarecurve(LUCAS, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(LUCAS),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
LUCAS<-b

##Zobel
Zobel1<-paste0(metadata2$site[metadata2$method=="Zobel"],collapse="|")
Zobel<-fungi2[,grepl(Zobel1,names(fungi2))]
Zobel<-as.data.frame(t(Zobel))
p<-rarecurve(Zobel, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(Zobel),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
Zobel<-b
##deep
deep1<-paste0(metadata2$site[metadata2$method=="deep"],collapse="|")
deep<-fungi2[,grepl(deep1,names(fungi2))]
deep<-as.data.frame(t(deep))
p<-rarecurve(deep, step = 20,label = F)

library(tidyr)
library(purrr)
library(stringr)
b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(deep),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
deep<-b


##deep_SUCC
deep_SUCC1<-paste0(metadata2$site[metadata2$method=="deep_SUCC"],collapse="|")
deep_SUCC<-fungi2[,grepl(deep_SUCC1,names(fungi2))]
deep_SUCC<-as.data.frame(t(deep_SUCC))
p<-rarecurve(deep_SUCC, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(deep_SUCC),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
deep_SUCC<-b

##SUCC
SUCC1<-paste0(metadata2$site[metadata2$method=="SUCC"],collapse="|")
SUCC<-fungi2[,grepl(SUCC1,names(fungi2))]
SUCC<-as.data.frame(t(SUCC))
p<-rarecurve(SUCC, step = 20,label = F)

b<-map_dfr(p,bind_rows) %>%
  bind_cols(Group=rownames(SUCC),.) %>%
  pivot_longer(-Group)%>%
  drop_na()%>%
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>%
  select(-name)

b$site<-substr(b$Group,1,2)
b$type<-paste0(b$site,substr(b$Group,4,nchar(b$Group)))
b$type1<-substr(b$Group,4,nchar(b$Group))
SUCC<-b
####
deep2$method<-"deep"
LUCAS2$method<-"LUCAS"
DarkDiv2$method<-"DarkDiv"
Zobel2$ratio<-Zobel2$value/Zobel2$unpooledvalue  
Zobel2$method<-"Zobel"
GSMC2$method<-"GSMc_62"
GSMC3$method<-"GSMc_40A"
GSMC4$method<-"GSMc_40B"
deep_SUCC2$method<-"deep_SUCC"
SUCC2$method<-"SUCC"
all<-rbind(SUCC2,DarkDiv2,Zobel2,GSMC2,GSMC3,GSMC4,deep2,deep_SUCC2,LUCAS2)
write.csv(all,"proportion.fungi.richness.summarydepthforpooled.csv")
##
SUCC$method<-"SUCC"
DarkDiv$method<-"DarkDiv"
Zobel$method<-"Zobel"
GSMC62$method<-"GSMC62"
GSMc_40A$method<-"GSMc40A"
GSMc_40B$method<-"GSMc40B"
deep$method<-"deep"
deep_SUCC$method<-"deep_SUCC"
LUCAS$method<-"LUCAS"
all1<-rbind(SUCC,DarkDiv,Zobel,GSMC62,GSMc_40A,GSMc_40B,deep,deep_SUCC,LUCAS)
all1<-all1[all1$n_seqs>14409,]
##
all2<-all %>%
  group_by(type1,method) %>%
  summarize(mean.ratio=mean(ratio),sd=sd(ratio),.groups = "keep")

all3<-all %>%
  group_by(type1) %>%
  summarize(mean.ratio=mean(ratio),sd=sd(ratio))

all2<-all[,c(1,6,9,10)]
##DNA
DNA<-all2[grepl("DNA",all2$type1),]
mod01<-lmer(ratio ~ method+ (1|site),data = DNA)
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
  labs(x="method",y="Pooling effect")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5))
##soil
soil<-all2[grepl("soil",all2$type1),]
mod01<-lmer(ratio ~ method+ (1|site),data = soil)
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
  labs(x="method",y="Pooling effect")+
  theme_light() +
  scale_x_discrete(limits=as.character(c( "GSMc_62","GSMc_40A","GSMc_40B","Zobel","DarkDiv","LUCAS","SUCC","deep_SUCC","deep")))+ 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )+
  geom_hline(yintercept=1, linetype="dashed", color = "grey")+
  geom_text(data = y.site, aes(x = area , y = ymax, label = letter,hjust=-0.5)) 