#####################NMDS mit Bray Curtis als Distanzmaß
#zero counts entfernen
water_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start, meta in sample name enthalten, nur wasser
separation<-colSums(water_subset2)!=0
water_subset_nonzero<-water_subset2[,separation]
ncol(water_subset2)
ncol(water_subset_nonzero)
names(water_subset_nonzero)
head(water_subset_nonzero)

#2 dimensionen
nmds_water<-metaMDS(water_subset_nonzero,try=100,autotransform=FALSE)
#extract data
names(nmds_water)
nmds_water_coord<-as.data.frame(nmds_water$points)
str(nmds_water_coord)
row.names(nmds_water_coord)<-row.names(water_subset_nonzero)
plot(nmds_water_coord)
text(nmds_water_coord,labels=row.names(nmds_water_coord),cex=0.7) #gives the sample ids, must be called when the plot is already open

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
meta_nmds_water<-merge(nmds_water_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_nmds_water)
#make "row.names" to actual row.names
row.names(meta_nmds_water)<-meta_nmds_water$Row.names
meta_nmds_water<-meta_nmds_water[,-1]
head(meta_nmds_water)

###################nur cDNA
#zero counts entfernen
water_cdna_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start_cdna, meta in sample name enthalten, nur wasser
separation<-colSums(water_cdna_subset2)!=0
water_cdna_subset_nonzero<-water_cdna_subset2[,separation]
ncol(water_cdna_subset2)
ncol(water_cdna_subset_nonzero)
names(water_cdna_subset_nonzero)
head(water_cdna_subset_nonzero)

#2 dimensionen
nmds_water_cdna<-metaMDS(water_cdna_subset_nonzero,try=100,autotransform=FALSE)
#extract data
names(nmds_water_cdna)
nmds_water_cdna_coord<-as.data.frame(nmds_water_cdna$points)
str(nmds_water_cdna_coord)
row.names(nmds_water_cdna_coord)<-row.names(water_cdna_subset_nonzero)
plot(nmds_water_cdna_coord)
text(nmds_water_cdna_coord,labels=row.names(nmds_water_cdna_coord),cex=0.7) #gives the sample ids, must be called when the plot is already open

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
meta_nmds_water_cdna<-merge(nmds_water_cdna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_nmds_water_cdna)
#make "row.names" to actual row.names
row.names(meta_nmds_water_cdna)<-meta_nmds_water_cdna$Row.names
meta_nmds_water_cdna<-meta_nmds_water_cdna[,-1]
head(meta_nmds_water_cdna)
#combining dna and cdna, marking samples for highlighting
write.csv(meta_nmds_water_dna,file="meta_nmds_water_dna.csv")
write.csv(meta_nmds_water_dna,file="meta_nmds_water_dna.csv")

combined_nmds<-read.csv(file.choose(),sep=";",row.names=1)

##plotting single in ggplot2
path<-meta_nmds_water_cdna
path$time<-as.numeric(as.character(path$time))
path<-path[order(path$time),]

ggplot() + 
  geom_point(data=meta_nmds_water_cdna,aes(x=MDS1,y=MDS2,shape=treatment,colour=time),size=4) + # add the point markers
  geom_path(data=path,aes(x=MDS1,y=MDS2,group=treatment,colour=as.numeric(as.character(time))),size=1.5)+
  geom_text(data=meta_nmds_water_cdna,aes(x=MDS1,y=MDS2,label=time,size=5,vjust=-0.7))+

#plotting combined

path_both<-combined_nmds
path_both$time<-as.numeric(as.character(path_both$time))
path_both<-path_both[order(path_both$time),]

##gefärbte punkte
ggplot() + 
  geom_point(data=subset(combined_nmds,group=="normal"),aes(x=MDS1,y=MDS2,shape=treatment,colour=time),size=2.5) +
  geom_point(data=subset(combined_nmds,group=="bloom"),aes(x=MDS1,y=MDS2,shape=treatment,colour=time),size=4,colour="red")+
  geom_point(data=subset(combined_nmds,group=="second"),aes(x=MDS1,y=MDS2,shape=treatment,colour=time),size=4,colour="green")+
  geom_path(data=path_both,aes(x=MDS1,y=MDS2,group=treatment,colour=as.numeric(as.character(time))),size=1)+
  geom_text(data=combined_nmds,aes(x=MDS1,y=MDS2,label=time,size=5,vjust=-0.7))+
  coord_equal() +
  facet_wrap(~nucleic_acid)

##gefärbte linien
ggplot() + 
  geom_path(data=path_both,aes(x=MDS1,y=MDS2,group=treatment,colour=group),size=2)+
  geom_point(data=combined_nmds,aes(x=MDS1,y=MDS2,shape=treatment),size=2.5)+
  geom_text(data=combined_nmds,aes(x=MDS1,y=MDS2,label=time,vjust=-0.7),size=3)+
  scale_colour_discrete(limits=levels(path_both$group))+
  coord_equal() +
  facet_wrap(~nucleic_acid)
  #theme_bw()
ggsave(file="nmds_water_paper_test.jpg", width=10, height=5)

	
##################nur DNA
#zero counts entfernen
water_dna_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start_dna, meta in sample name enthalten, nur wasser
separation<-colSums(water_dna_subset2)!=0
water_dna_subset_nonzero<-water_dna_subset2[,separation]
ncol(water_dna_subset2)
ncol(water_dna_subset_nonzero)
names(water_dna_subset_nonzero)
head(water_dna_subset_nonzero)

#2 dimensionen
nmds_water_dna<-metaMDS(water_dna_subset_nonzero,try=100,autotransform=FALSE)
#extract data
names(nmds_water_dna)
nmds_water_dna_coord<-as.data.frame(nmds_water_dna$points)
str(nmds_water_dna_coord)
row.names(nmds_water_dna_coord)<-row.names(water_dna_subset_nonzero)
plot(nmds_water_dna_coord)
text(nmds_water_dna_coord,labels=row.names(nmds_water_dna_coord),cex=0.7) #gives the sample ids, must be called when the plot is already open

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
meta_nmds_water_dna<-merge(nmds_water_dna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_nmds_water_dna)
#make "row.names" to actual row.names
row.names(meta_nmds_water_dna)<-meta_nmds_water_dna$Row.names
meta_nmds_water_dna<-meta_nmds_water_dna[,-1]
head(meta_nmds_water_dna)
##########################
#plotting in ggplot2




######################################################

#removing nmds coords for envfit
meta_nmds_water2<-meta_nmds_water[,3:10]
str(meta_nmds_water2)
#keine cellcounts, kein treatment, keine parallele, kein time, kein habitat
meta_nmds_water3<-meta_nmds_water2[,c(2,7,8)]
str(meta_nmds_water3)

#turning "time" into factor (after creating envfit-file)
meta_nmds_water$time<-as.factor(meta_nmds_water$time)

#subset(meta_data,nucleic_acid=="cdna"&habitat=="water"&treatment=="glyph")

#plotting (using color ramp for time)
colfunc <- colorRampPalette(c("yellow2","red","purple"))
colfunc(16)
plot(nmds_water,main="NMDS_water", type="p",display="sites")#,xlim=c(-1,1),ylim=c(-0.8,0.6))
with(meta_nmds_water, points(nmds_water,col=colfunc(16)[meta_nmds_water$time],pch=c(16,17)[meta_nmds_water$treatment]))
with(meta_nmds_water, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_nmds_water, legend("topright", legend = levels(time), pch=16, bty = "n", col=colfunc(16)))
text(nmds_water_coord,labels=row.names(nmds_water_coord),cex=0.7)
#ordihull(nmds_water,groups=interaction(meta_nmds_water$time,meta_nmds_water$treatment),label=TRUE,cex=0.4,draw="polygon",col="red") #falsche zuordnung
water_envfit<-envfit(nmds_water,meta_nmds_water3,permutations=100)
plot(water_envfit)


#plotting (using two colors for treatment, control)
plot(nmds_water,main="NMDS_water", type="p",display="sites")
with(meta_nmds_water, points(nmds_water,col=c("blue","red")[meta_nmds_water$treatment],pch=c(16,17)[meta_nmds_water$treatment]))
with(meta_nmds_water, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_nmds_water, legend("topright", legend = levels(time), pch=16, bty = "n", col=colfunc(16)))
ordihull(nmds_water,groups=interaction(meta_nmds_water$time,meta_nmds_water$treatment),label=FALSE,cex=0.4,draw="polygon",col="green")
#water_envfit<-envfit(nmds_water,meta_nmds_water3,permutations=100)
#plot(water_envfit)
text(nmds_water_coord,labels=row.names(nmds_water_coord),cex=0.5)