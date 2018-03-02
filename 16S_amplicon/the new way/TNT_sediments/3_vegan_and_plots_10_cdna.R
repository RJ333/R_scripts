#vegan analysis and plots
##for vegan
##otus in spalten/proben als rownames 
tnt10_tvegan_cdna<-read.csv(file.choose(),row.names=1,sep=";")

#transpose back for merging
tnt10_vegan_cdna<-t(tnt10_tvegan_cdna)
head(tnt10_vegan_cdna)

#2 dimensionen 
nmds_tnt10_cdna_notrans<-metaMDS(tnt10_vegan_cdna,try=100,autotransform=FALSE)
nmds_tnt10_cdna<-nmds_tnt10_cdna_notrans
#extract data
names(nmds_tnt10_cdna)
nmds_tnt10_cdna_coord<-as.data.frame(nmds_tnt10_cdna$points)
str(nmds_tnt10_cdna_coord)
row.names(nmds_tnt10_cdna_coord)<-row.names(tnt10_vegan_cdna)
plot(nmds_tnt10_cdna_coord)
text(nmds_tnt10_cdna_coord,labels=row.names(nmds_tnt10_cdna_coord),cex=0.5) #gives the sample ids, must be called when the plot is already open
#merge with meta data
meta_nmds_tnt10_cdna<-merge(nmds_tnt10_cdna_coord,meta_data_tnt10,by="row.names",all.x=TRUE)
head(meta_nmds_tnt10_cdna)

#make "row.names" to actual row.names
row.names(meta_nmds_tnt10_cdna)<-meta_nmds_tnt10_cdna$Row.names
meta_nmds_tnt10_cdna<-meta_nmds_tnt10_cdna[,-1]
head(meta_nmds_tnt10_cdna)

ggplot() + 
  geom_point(data=subset(meta_nmds_tnt10_cdna,cruise=="L17_07"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=2.5,colour="blue") +
  geom_point(data=subset(meta_nmds_tnt10_cdna,cruise=="L17_08"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=4,colour="red")+
  geom_point(data=subset(meta_nmds_tnt10_cdna,cruise=="L16_14"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=4,colour="green")+
  geom_point(data=subset(meta_nmds_tnt10_cdna,station=="Mo7"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel)),size=4,colour="black")+
  geom_text(data=meta_nmds_tnt10_cdna,aes(x=MDS1,y=MDS2,label=station,size=3,vjust=-1))+
  coord_equal()+
  theme_bw()+
  theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	theme(legend.position='none')
  #facet_wrap(~nucleic_acid)
ggsave(file="tnt10_nmds_cdna.png", width=14, height=14)

#anosim 
anosim_treatment_water_cdna<-anosim(cast_tnt_otu10_table3,meta_nmds_water_cdna$treatment,permutations=999,distance="bray")
plot(anosim_treatment)


#cca
str(meta_water_ccdna)
meta_water_ccdna$days<-as.factor(meta_water_ccdna$days)
str(meta_water_ccdna)
cca_water_ccdna<-cca(water_ccdna_subset_nonzero~disturbance + days,data=meta_water_ccdna)

colfunc <- colorRampPalette(c("yellow2","red","purple"))
colfunc(16)
plot(cca_water_ccdna,main="CCA_water_ccdna", type="p",display="sites")
with(meta_water_ccdna, points(cca_water_ccdna,col=colfunc(16)[meta_water_ccdna$days],pch=c(16,17)[meta_water_ccdna$treatment]))
with(meta_water_ccdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_ccdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(cca_water_ccdna,groups=interaction(meta_water_ccdna$days,meta_water_ccdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_ccdna,labels=row.names(meta_water_ccdna),cex=0.7)

names(cca_water_ccdna)
summary(cca_water_ccdna)

#dca
dca_water_ccdna<-decorana(water_ccdna_subset_nonzero)
plot(dca_water_ccdna,type="p",display="sites",main="dca_water_ccdna")
with(meta_water_ccdna, points(dca_water_ccdna,col=colfunc(16)[meta_water_ccdna$days],pch=c(16,17)[meta_water_ccdna$treatment]))
with(meta_water_ccdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_ccdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(dca_water_ccdna,groups=interaction(meta_water_ccdna$days,meta_water_ccdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_ccdna,labels=row.names(meta_water_ccdna),cex=0.7)

names(dca_water_ccdna)
summary(dca_water_ccdna)

#rda
rda_water_ccdna<-rda(water_ccdna_subset_nonzero~days+disturbance,meta_water_ccdna)
plot(rda_water_ccdna,type="p",display="sites",main="rda_water_ccdna")
with(meta_water_ccdna, points(rda_water_ccdna,col=colfunc(16)[meta_water_ccdna$days],pch=c(16,17)[meta_water_ccdna$treatment]))
with(meta_water_ccdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_ccdna, legend("bottomleft", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(rda_water_ccdna,groups=interaction(meta_water_ccdna$days,meta_water_ccdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
names(rda_water_ccdna)
summary(rda_water_ccdna)

#pca
pca_water_ccdna<-rda(water_ccdna_subset_nonzero)
plot(pca_water_ccdna,main="PCA_water_ccdna", type="p",display="sites")
with(meta_water_ccdna, points(pca_water_ccdna,col=colfunc(16)[meta_water_ccdna$days],pch=c(16,17)[meta_water_ccdna$treatment]))
with(meta_water_ccdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_ccdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(pca_water_ccdna,groups=interaction(meta_water_ccdna$days,meta_water_ccdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_ccdna,labels=row.names(meta_water_ccdna),cex=0.7)

#pca
pca_water_ccdna_scaled<-rda(water_ccdna_subset_nonzero,scale = TRUE)
plot(pca_water_ccdna_scaled,main="PCA_water_ccdna_scaled", type="p",display="sites")
with(meta_water_ccdna, points(pca_water_ccdna_scaled,col=colfunc(16)[meta_water_ccdna$days],pch=c(16,17)[meta_water_ccdna$treatment]))
with(meta_water_ccdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_ccdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(pca_water_ccdna_scaled,groups=interaction(meta_water_ccdna$disturbance,meta_water_ccdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
text(pca_water_ccdna_scaled,labels=row.names(meta_water_ccdna),cex=0.7)