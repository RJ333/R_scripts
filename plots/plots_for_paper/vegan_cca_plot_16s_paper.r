head(meta_nmds_water5_dna_ordered)
head(water5_cdna_subset_nonzero)

######################nur cDNA
#CCA
cca_water5_cdna_meta<-cca(water5_cdna_subset_nonzero~glyph,meta_nmds_water5_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water5_cdna_meta,type="p", display="sites")
text(cca_water5_cdna_meta,labels=row.names(water5_cdna_subset_nonzero),cex=0.7)

#extract data
cca_water5_cdna_coord<-as.data.frame(scores(cca_water5_cdna_meta)$sites)
str(cca_water5_cdna_coord)


#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
#compare files
meta_file_vegan
cca_water5_cdna_coord
meta_cca_water5_cdna<-merge(cca_water5_cdna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_cca_water5_cdna)
#make "row.names" to actual row.names
row.names(meta_cca_water5_cdna)<-meta_cca_water5_cdna$Row.names
meta_cca_water5_cdna<-meta_cca_water5_cdna[,-1]
head(meta_cca_water5_cdna)

##################nur DNA
#CCA
cca_water5_dna_meta<-cca(water5_dna_subset_nonzero~glyph,meta_nmds_water5_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water5_dna_meta,type="p", display="sites")
text(cca_water5_dna_meta,labels=row.names(water5_dna_subset_nonzero),cex=0.7)

#extract data
cca_water5_dna_coord<-as.data.frame(scores(cca_water5_dna_meta)$sites)
str(cca_water5_dna_coord)

#neue meta file notwendig, da neue row.names
meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
#compare files
meta_file_vegan
cca_water5_dna_coord
meta_cca_water5_dna<-merge(cca_water5_dna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_cca_water5_dna)
#make "row.names" to actual row.names
row.names(meta_cca_water5_dna)<-meta_cca_water5_dna$Row.names
meta_cca_water5_dna<-meta_cca_water5_dna[,-1]
head(meta_cca_water5_dna)

#combining dna and cdna, marking samples for highlighting
write.csv(meta_cca_water5_dna,file="meta_cca_water5_dna.csv")
write.csv(meta_cca_water5_cdna,file="meta_cca_water5_cdna.csv")

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
  geom_path(data=path_both5,aes(x=CCA1,y=CA1,group=treatment,colour=group),size=1)+
  geom_point(data=combined_cca5,aes(x=CCA1,y=CA1,shape=treatment),size=2.5)+
  geom_text(data=combined_cca5,aes(x=CCA1,y=CA1,label=days,vjust=-0.7),size=3)+
  scale_colour_discrete(limits=levels(path_both5$group))+
  scale_colour_manual(values=c("normal"="black","bloom"="orange"),
							name="phase",
							breaks=c("normal","bloom"),
							labels=c("normal","reaction"))+
  coord_equal() +
  facet_wrap(~nucleic_acid)
  #theme_bw()
ggsave(file="cca_water5_paper_test.jpg", width=10, height=5)
#welche mit metadaten beeinflusst? mod <- cca(varespec ~ Al + P + K, varechem)
#wie punkte extrahieren?
#"quali" beschreiben?
#wiederholen mit zeitpunkt 4 inklusive
#für cca wissen über glyphosatkonzentrationen notwendig...mit und ohne darstellen?2 x dna, 2 x cdna! ohne konz evtl in supplement?
#flip axis for comparable dna/cdna plots
m <- cca(dune)
plot(m, xlim = c(2, -3)

##anova einbauen?
anova(cca_water5_dna_meta) inertia cca1 dna 0.06427
anova(cca_water5_cdna_meta) intertia cca1 cdna 0.07661

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Number of permutations: 999

Model: cca(formula = water5_dna_subset_nonzero ~ glyph, data = meta_nmds_water5_dna_ordered, try = 100, autotransform = FALSE, test = "permutation", permutations = how(nperm = 2000))
         Df ChiSquare      F Pr(>F)   
Model     1   0.06427 4.6372  0.002 **
Residual 28   0.38809                 

Model: cca(formula = water5_cdna_subset_nonzero ~ glyph, data = meta_nmds_water5_dna_ordered, try = 100, autotransform = FALSE)
         Df ChiSquare      F Pr(>F)    
Model     1   0.07661 5.6647  0.001 ***
Residual 28   0.37867 