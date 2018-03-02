#vegan analysis and plots
##for vegan
##otus in spalten/proben als rownames 
tnt10_tvegan_dna<-read.csv(file.choose(),row.names=1,sep=";")

#transpose back for merging
tnt10_vegan_dna<-t(tnt10_tvegan_dna)
head(tnt10_vegan_dna)

#2 dimensionen
nmds_tnt10_dna_notrans<-metaMDS(tnt10_vegan_dna,try=100,autotransform=FALSE)

nmds_tnt10_dna<-nmds_tnt10_dna_notrans
#extract data
names(nmds_tnt10_dna)
nmds_tnt10_dna_coord<-as.data.frame(nmds_tnt10_dna$points)
str(nmds_tnt10_dna_coord)
row.names(nmds_tnt10_dna_coord)<-row.names(tnt10_vegan_dna)
plot(nmds_tnt10_dna_coord)
text(nmds_tnt10_dna_coord,labels=row.names(nmds_tnt10_dna_coord),cex=0.5) #gives the sample ids, must be called when the plot is already open
#merge with meta data
meta_nmds_tnt10_dna<-merge(nmds_tnt10_dna_coord,meta_data_tnt10,by="row.names",all.x=TRUE)
head(meta_nmds_tnt10_dna)
 
#make "row.names" to actual row.names
row.names(meta_nmds_tnt10_dna)<-meta_nmds_tnt10_dna$Row.names
meta_nmds_tnt10_dna<-meta_nmds_tnt10_dna[,-1]
head(meta_nmds_tnt10_dna)

ggplot() + 
  geom_point(data=subset(meta_nmds_tnt10_dna,cruise=="L17_07"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=2.5,colour="blue") +
  geom_point(data=subset(meta_nmds_tnt10_dna,cruise=="L17_08"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=4,colour="red")+
  geom_point(data=subset(meta_nmds_tnt10_dna,cruise=="L16_14"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=4,colour="green")+
  geom_point(data=subset(meta_nmds_tnt10_dna,station=="Mo7"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel)),size=4,colour="yellow")+
  geom_text(data=meta_nmds_tnt10_dna,aes(x=MDS1,y=MDS2,label=station,size=3,vjust=-1))+
  coord_equal()
  #facet_wrap(~nucleic_acid)
ggsave(file="tnt10_nmds_dna.png", width=12, height=12)

#anosim 
anosim_treatment_water_dna<-anosim(cast_tnt_otu10_table3,meta_nmds_water_dna$treatment,permutations=999,distance="bray")
plot(anosim_treatment)


#cca
str(meta_water_cdna)
meta_water_cdna$days<-as.factor(meta_water_cdna$days)
str(meta_water_cdna)
cca_water_cdna<-cca(water_cdna_subset_nonzero~disturbance + days,data=meta_water_cdna)

colfunc <- colorRampPalette(c("yellow2","red","purple"))
colfunc(16)
plot(cca_water_cdna,main="CCA_water_cdna", type="p",display="sites")
with(meta_water_cdna, points(cca_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(cca_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_cdna,labels=row.names(meta_water_cdna),cex=0.7)

names(cca_water_cdna)
summary(cca_water_cdna)

#dca
dca_water_cdna<-decorana(water_cdna_subset_nonzero)
plot(dca_water_cdna,type="p",display="sites",main="dca_water_cdna")
with(meta_water_cdna, points(dca_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(dca_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_cdna,labels=row.names(meta_water_cdna),cex=0.7)

names(dca_water_cdna)
summary(dca_water_cdna)

#rda
rda_water_cdna<-rda(water_cdna_subset_nonzero~days+disturbance,meta_water_cdna)
plot(rda_water_cdna,type="p",display="sites",main="rda_water_cdna")
with(meta_water_cdna, points(rda_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("bottomleft", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(rda_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
names(rda_water_cdna)
summary(rda_water_cdna)

#pca
pca_water_cdna<-rda(water_cdna_subset_nonzero)
plot(pca_water_cdna,main="PCA_water_cdna", type="p",display="sites")
with(meta_water_cdna, points(pca_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(pca_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_cdna,labels=row.names(meta_water_cdna),cex=0.7)

#pca
pca_water_cdna_scaled<-rda(water_cdna_subset_nonzero,scale = TRUE)
plot(pca_water_cdna_scaled,main="PCA_water_cdna_scaled", type="p",display="sites")
with(meta_water_cdna, points(pca_water_cdna_scaled,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(pca_water_cdna_scaled,groups=interaction(meta_water_cdna$disturbance,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
text(pca_water_cdna_scaled,labels=row.names(meta_water_cdna),cex=0.7)