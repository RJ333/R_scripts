##control zeitpunkt 4 macht plot darstellung kaputt

#ohne umweltparameter
cca_water_cdna<-cca(water_cdna_subset_nonzero,try=100,autotransform=FALSE)
plot(cca_water_cdna,type="p", display="sites")
text(cca_water_cdna,labels=row.names(water_cdna_subset_nonzero),cex=0.7)

#mit umweltparameter 
meta_nmds_water_dna_ordered<-read.csv(file.choose(),sep="\t",row.names=1)
######################nur cDNA
#CCA
cca_water_cdna_meta<-cca(water_cdna_subset_nonzero~glyph,meta_nmds_water_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water_cdna_meta,type="p", display="sites")
text(cca_water_cdna_meta,labels=row.names(water_cdna_subset_nonzero),cex=0.7)

#extract data
cca_water_cdna_coord<-as.data.frame(scores(cca_water_cdna_meta)$sites)
str(cca_water_cdna_coord)


#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
#compare files
meta_file_vegan
cca_water_cdna_coord
meta_cca_water_cdna<-merge(cca_water_cdna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_cca_water_cdna)
#make "row.names" to actual row.names
row.names(meta_cca_water_cdna)<-meta_cca_water_cdna$Row.names
meta_cca_water_cdna<-meta_cca_water_cdna[,-1]
head(meta_cca_water_cdna)

##################nur DNA
#CCA
cca_water_dna_meta<-cca(water_dna_subset_nonzero~glyph,meta_nmds_water_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water_dna_meta,type="p", display="sites",ylim=c(-0.5,0.5))
text(cca_water_dna_meta,labels=row.names(water_dna_subset_nonzero),cex=0.7)

#extract data
cca_water_dna_coord<-as.data.frame(scores(cca_water_dna_meta)$sites)
str(cca_water_dna_coord)

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
#compare files
meta_file_vegan
cca_water_dna_coord
meta_cca_water_dna<-merge(cca_water_dna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_cca_water_dna)
#make "row.names" to actual row.names
row.names(meta_cca_water_dna)<-meta_cca_water_dna$Row.names
meta_cca_water_dna<-meta_cca_water_dna[,-1]
head(meta_cca_water_dna)

#combining dna and cdna, marking samples for highlighting
write.csv(meta_cca_water_dna,file="meta_cca_water_dna.csv")
write.csv(meta_cca_water_cdna,file="meta_cca_water_cdna.csv")

#change group (bloom) and color to orange
combined_cca5<-read.csv(file.choose(),sep="\t",row.names=1) #bei openoffice

#plotting combined

path_both5<-combined_cca5
path_both5$time<-as.numeric(as.character(path_both5$time))
path_both5$days<-as.numeric(as.character(path_both5$days))
path_both5<-path_both5[order(path_both5$days),]

##gefärbte punkte
ggplot() + 
  geom_point(data=subset(combined_cca5,group=="normal"),aes(x=CCA1,y=CA1,shape=treatment,colour="normal"),size=2.5)+
  geom_point(data=subset(combined_cca5,group=="bloom"),aes(x=CCA1,y=CA1,shape=treatment,colour="bloom"),size=4)+
  geom_path(data=path_both5,aes(x=CCA1,y=CA1,group=treatment,colour=as.numeric(as.character(days))),size=1)+
  geom_text(data=combined_cca5,aes(x=CCA1,y=CA1,label=days,size=5,vjust=-0.7))+
  coord_equal() +
  facet_wrap(~nucleic_acid)

##gefärbte linien
ggplot() + 
  geom_path(data=path_both5,aes(x=CCA1,y=CA1,group=treatment,colour=group),size=2)+
  geom_point(data=combined_cca5,aes(x=CCA1,y=CA1,shape=treatment),size=2.5)+
  geom_text(data=combined_cca5,aes(x=CCA1,y=CA1,label=days,vjust=-0.7),size=3)+
  scale_colour_discrete(limits=levels(path_both5$group))+
  coord_equal() +
  facet_wrap(~nucleic_acid)
  #theme_bw()
ggsave(file="cca_water_paper_test.jpg", width=10, height=5)
#welche mit metadaten beeinflusst? mod <- cca(varespec ~ Al + P + K, varechem)
#wie punkte extrahieren?
#"quali" beschreiben?
#wiederholen mit zeitpunkt 4 inklusive
#für cca wissen über glyphosatkonzentrationen notwendig...mit und ohne darstellen?2 x dna, 2 x cdna! ohne konz evtl in supplement?
#flip axis for comparable dna/cdna plots
m <- cca(dune)
plot(m, xlim = c(2, -3)
