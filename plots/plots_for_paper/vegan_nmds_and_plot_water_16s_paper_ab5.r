#####################NMDS mit Bray Curtis als Distanzmaß

###################nur cDNA
#zero counts entfernen
water5_cdna_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start_cdna_ab5, meta in sample name enthalten, nur wasser
separation<-colSums(water5_cdna_subset2)!=0
water5_cdna_subset_nonzero<-water5_cdna_subset2[,separation]
ncol(water5_cdna_subset2)
ncol(water5_cdna_subset_nonzero)
names(water5_cdna_subset_nonzero)
head(water5_cdna_subset_nonzero)

#2 dimensionen
nmds_water5_cdna<-metaMDS(water5_cdna_subset_nonzero,try=100,autotransform=FALSE)
#extract data
names(nmds_water5_cdna)
nmds_water5_cdna_coord<-as.data.frame(nmds_water5_cdna$points)
str(nmds_water5_cdna_coord)
row.names(nmds_water5_cdna_coord)<-row.names(water5_cdna_subset_nonzero)
plot(nmds_water5_cdna_coord)
text(nmds_water5_cdna_coord,labels=row.names(nmds_water5_cdna_coord),cex=0.7) #gives the sample ids, must be called when the plot is already open

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
meta_nmds_water5_cdna<-merge(nmds_water5_cdna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_nmds_water5_cdna)
#make "row.names" to actual row.names
row.names(meta_nmds_water5_cdna)<-meta_nmds_water5_cdna$Row.names
meta_nmds_water5_cdna<-meta_nmds_water5_cdna[,-1]
head(meta_nmds_water5_cdna)

##################nur DNA
#zero counts entfernen
water5_dna_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start_dna_ab5, meta in sample name enthalten, nur wasser
separation<-colSums(water5_dna_subset2)!=0
water5_dna_subset_nonzero<-water5_dna_subset2[,separation]
ncol(water5_dna_subset2)
ncol(water5_dna_subset_nonzero)
names(water5_dna_subset_nonzero)
head(water5_dna_subset_nonzero)

#2 dimensionen
nmds_water5_dna<-metaMDS(water5_dna_subset_nonzero,try=100,autotransform=FALSE)
#extract data
names(nmds_water5_dna)
nmds_water5_dna_coord<-as.data.frame(nmds_water5_dna$points)
str(nmds_water5_dna_coord)
row.names(nmds_water5_dna_coord)<-row.names(water5_dna_subset_nonzero)
plot(nmds_water5_dna_coord)
text(nmds_water5_dna_coord,labels=row.names(nmds_water5_dna_coord),cex=0.7) #gives the sample ids, must be called when the plot is already open

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
meta_nmds_water5_dna<-merge(nmds_water5_dna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_nmds_water5_dna)
#make "row.names" to actual row.names
row.names(meta_nmds_water5_dna)<-meta_nmds_water5_dna$Row.names
meta_nmds_water5_dna<-meta_nmds_water5_dna[,-1]
head(meta_nmds_water5_dna)

#combining dna and cdna, marking samples for highlighting
write.csv(meta_nmds_water5_dna,file="meta_nmds_water5_dna.csv")
write.csv(meta_nmds_water5_cdna,file="meta_nmds_water5_cdna.csv")

combined_nmds5<-read.csv(file.choose(),sep=";",row.names=1)

#plotting combined

path_both5<-combined_nmds5
path_both5$time<-as.numeric(as.character(path_both5$time))
path_both5$days<-as.numeric(as.character(path_both5$days))
path_both5<-path_both5[order(path_both5$days),]

##gefärbte punkte
ggplot() + 
  geom_point(data=subset(combined_nmds5,group=="normal"),aes(x=MDS1,y=MDS2,shape=treatment,colour=days),size=2.5) +
  geom_point(data=subset(combined_nmds5,group=="bloom"),aes(x=MDS1,y=MDS2,shape=treatment,colour=days),size=4,colour="red")+
  geom_point(data=subset(combined_nmds5,group=="second"),aes(x=MDS1,y=MDS2,shape=treatment,colour=days),size=4,colour="green")+
  geom_path(data=path_both5,aes(x=MDS1,y=MDS2,group=treatment,colour=as.numeric(as.character(days))),size=1)+
  geom_text(data=combined_nmds5,aes(x=MDS1,y=MDS2,label=days,size=5,vjust=-0.7))+
  coord_equal() +
  facet_wrap(~nucleic_acid)

##gefärbte linien
ggplot() + 
  geom_path(data=path_both5,aes(x=MDS1,y=MDS2,group=treatment,colour=group),size=2)+
  geom_point(data=combined_nmds5,aes(x=MDS1,y=MDS2,shape=treatment),size=2.5)+
  geom_text(data=combined_nmds5,aes(x=MDS1,y=MDS2,label=days,vjust=-0.7),size=3)+
  scale_colour_discrete(limits=levels(path_both5$group))+
  coord_equal() +
  facet_wrap(~nucleic_acid)
  #theme_bw()
ggsave(file="nmds_water5_paper_test.jpg", width=10, height=5)

	

####test ordisurf
#check order
meta_nmds_water5_dna_ordered<-read.csv(file.choose(),sep=";",row.names=1)
nmds_water5_dna$points
meta_nmds_water5_dna_ordered
ordisurf(nmds_water5_dna,meta_nmds_water5_dna_ordered$glyph,main="hallo",col="blue")
text(nmds_water5_dna_coord,labels=row.names(nmds_water5_dna_coord),cex=0.7)
stressplot(nmds_water5_dna)

ordiplot(nmds_water5_dna,type="n")
orditorp(nmds_water5_dna,display="sites",air=0.01,cex=1)
ordicluster(nmds_water5_dna,hclust(vegdist(water5_dna_subset_nonzero,"bray")))



dca_water5_cdna<-decorana(water5_cdna_subset_nonzero)
plot(dca_water5_cdna,type="p", display="sites")
text(dca_water5_cdna,labels=row.names(water5_cdna_subset_nonzero),cex=0.7)

rda_water5_cdna<-rda(water5_cdna_subset_nonzero,try=100,autotransform=FALSE)
plot(rda_water5_cdna,type="p", display="sites")
text(rda_water5_cdna,labels=row.names(water5_cdna_subset_nonzero),cex=0.7)