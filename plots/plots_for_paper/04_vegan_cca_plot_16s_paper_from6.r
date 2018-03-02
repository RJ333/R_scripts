#zero counts entfernen
water_cdna_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start_cdna, meta in sample name enthalten, nur wasser
separation<-colSums(water_cdna_subset2)!=0
water_cdna_subset_nonzero<-water_cdna_subset2[,separation]
ncol(water_cdna_subset2)
ncol(water_cdna_subset_nonzero)
names(water_cdna_subset_nonzero)
head(water_cdna_subset_nonzero)

meta_nmds_water6_dna_ordered<-subset(meta_nmds_water5_dna_ordered,time>5)
water6_cdna_subset_nonzero<-read.csv(file.choose(),sep=";",row.names=1)
water6_dna_subset_nonzero<-read.csv(file.choose(),sep=";",row.names=1)

head(meta_nmds_water6_dna_ordered)
head(meta_nmds_water5_dna_ordered)
head(water6_cdna_subset_nonzero)

######################nur cDNA
#CCA
cca_water6_cdna_meta<-cca(water6_cdna_subset_nonzero~glyph,meta_nmds_water6_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water6_cdna_meta,type="p", display="sites")
text(cca_water6_cdna_meta,labels=row.names(water6_cdna_subset_nonzero),cex=0.7)

#test mit autotransform (looks exactly the same)
#cca_water6_cdna_meta_true<-cca(water6_cdna_subset_nonzero~glyph,meta_nmds_water6_dna_ordered,try=100,autotransform=TRUE)
#plot(cca_water6_cdna_meta_true,type="p", display="sites")
#text(cca_water6_cdna_meta_true,labels=row.names(water6_cdna_subset_nonzero),cex=0.7)

#extract data
cca_water6_cdna_coord<-as.data.frame(scores(cca_water6_cdna_meta)$sites)
str(cca_water6_cdna_coord)


#neue meta file notwendig, da neue row.names
#meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
#compare files
meta_file_vegan
cca_water6_cdna_coord
meta_cca_water6_cdna<-merge(cca_water6_cdna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_cca_water6_cdna)
#make "row.names" to actual row.names
row.names(meta_cca_water6_cdna)<-meta_cca_water6_cdna$Row.names
meta_cca_water6_cdna<-meta_cca_water6_cdna[,-1]
head(meta_cca_water6_cdna)

##################nur DNA
#CCA
cca_water6_dna_meta<-cca(water6_dna_subset_nonzero~glyph,meta_nmds_water6_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water6_dna_meta,type="p", display="sites")
text(cca_water6_dna_meta,labels=row.names(water6_dna_subset_nonzero),cex=0.7)

#extract data
cca_water6_dna_coord<-as.data.frame(scores(cca_water6_dna_meta)$sites)
str(cca_water6_dna_coord)

#neue meta file notwendig, da neue row.names
#meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
#compare files
meta_file_vegan
cca_water6_dna_coord
meta_cca_water6_dna<-merge(cca_water6_dna_coord,meta_file_vegan,by="row.names",all.x=TRUE)
head(meta_cca_water6_dna)
#make "row.names" to actual row.names
row.names(meta_cca_water6_dna)<-meta_cca_water6_dna$Row.names
meta_cca_water6_dna<-meta_cca_water6_dna[,-1]
head(meta_cca_water6_dna)

#combining dna and cdna, marking samples for highlighting
write.csv(meta_cca_water6_dna,file="meta_cca_water6_dna.csv")
write.csv(meta_cca_water6_cdna,file="meta_cca_water6_cdna.csv")

#change group (bloom) and color to orange
combined_cca6<-read.csv(file.choose(),sep=";",row.names=1)

#plotting combined

path_both6<-combined_cca6
path_both6$time<-as.numeric(as.character(path_both6$time))
path_both6$days<-as.numeric(as.character(path_both6$days))
path_both6$new_day<-as.numeric(as.character(path_both6$new_day))
path_both6<-path_both6[order(path_both6$days),]

##gefärbte punkte
ggplot() + 
  geom_point(data=subset(combined_cca6,group=="normal"),aes(x=CCA1,y=CA1,shape=treatment,colour="normal"),size=2.5)+
  geom_point(data=subset(combined_cca6,group=="bloom"),aes(x=CCA1,y=CA1,shape=treatment,colour="bloom"),size=4)+
  geom_path(data=path_both6,aes(x=CCA1,y=CA1,group=treatment,colour=as.numeric(as.character(new_day))),size=1)+
  geom_text(data=combined_cca6,aes(x=CCA1,y=CA1,label=new_day,size=5,vjust=-0.7))+
  coord_equal() +
  facet_wrap(~nucleic_acid)

##gefärbte linien ab zeitpunkt6
path_both6<-read.csv(file.choose(),sep=";",row.names=1)
path_both6<-path_both6[with(path_both6, order(nucleic_acid)), ]
path_both6$nucleic_acid2<-factor(path_both6$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rDNA amplicons"))
combined_cca6$nucleic_acid2<-factor(combined_cca6$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rDNA amplicons"))
ggplot() + 
  geom_path(data=path_both6,aes(x=CCA1,y=CA1,group=treatment,colour=group),size=1.5)+
  geom_point(data=combined_cca6,aes(x=CCA1,y=CA1,shape=treatment),size=3,colour="black",fill="white")+
  geom_text(data=combined_cca6,aes(x=CCA1,y=CA1,label=new_day,vjust=-1),size=3)+
  scale_colour_discrete(limits=levels(path_both6$group))+
  scale_colour_manual(values=c("normal"="grey20","bloom"="grey60"),
							name="Community succession",
							breaks=c("normal","bloom"),
							labels=c("normal",'"changed"'))+
  scale_shape_manual(values=c("control"=21,"glyph"=22),
							name="Microcosm",
							breaks=c("control","glyph"),
							labels=c("Control","Treatment"))+
  coord_equal() +
  facet_wrap(~nucleic_acid2)+
  theme_bw()
ggsave(file="cca_water6_paper_test.png", width=12, height=7)
	facet_wrap(~material2,ncol=2,labeller=label_value)


#########gefärbte linien ab zeitpunkt 5
combined_cca5<-read.csv(file.choose(),sep=";",row.names=1)
path_both5<-read.csv(file.choose(),sep=";",row.names=1)
path_both5<-path_both5[with(path_both5, order(nucleic_acid)), ]
path_both5$nucleic_acid2<-factor(path_both5$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rDNA amplicons"))
combined_cca5$nucleic_acid2<-factor(combined_cca5$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rDNA amplicons"))
ggplot() + 
  geom_path(data=path_both5,aes(x=CCA1,y=CA1,group=treatment,colour=group),size=1.5)+
  geom_point(data=combined_cca5,aes(x=CCA1,y=CA1,shape=treatment),size=3,colour="black",fill="white")+
  geom_text(data=combined_cca5,aes(x=CCA1,y=CA1,label=new_day,vjust=-1),size=3)+
  scale_colour_discrete(limits=levels(path_both5$group))+
  scale_colour_manual(values=c("normal"="grey20","bloom"="grey50"),
							name="Community succession",
							breaks=c("normal","bloom"),
							labels=c("normal",'"changed"'))+
  scale_shape_manual(values=c("control"=21,"glyph"=22),
							name="Microcosm",
							breaks=c("control","glyph"),
							labels=c("Control","Treatment"))+
  coord_equal() +
  facet_wrap(~nucleic_acid2)+
  theme_bw()
ggsave(file="cca_water5_paper_test.png", width=12, height=7)
	facet_wrap(~material2,ncol=2,labeller=label_value)


##gefärbte linien + reverse y axis (gefärbte linien liegen falsch?)
ggplot() + 
  geom_path(data=path_both6,aes(x=CCA1,y=CA1,group=treatment,colour=group),size=1)+
  geom_point(data=combined_cca6,aes(x=CCA1,y=CA1,shape=treatment),size=2.5)+
  geom_text(data=combined_cca6,aes(x=CCA1,y=CA1,label=new_day,vjust=-0.7),size=3)+
  scale_colour_discrete(limits=levels(path_both6$group))+
  scale_colour_manual(values=c("normal"="black","bloom"="orange"),
							name="phase",
							breaks=c("normal","bloom"),
							labels=c("normal","reaction"))+
  coord_equal() +
  facet_wrap(~nucleic_acid)+
  scale_y_reverse(data=subset(path_both6,nucleic_acid=="cDNA")) # turns of course both plots :-\

##anova einbauen?
anova(cca_water6_dna_meta) inertia cca1 dna (ab t5) 0.06427 (ab t6) 0.08063
anova(cca_water6_cdna_meta) inertia cca1 cdna (ab t5) 0.07661 (ab t6) 0.08784

Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Number of permutations: 999

##ab t5
Model: cca(formula = water6_dna_subset_nonzero ~ glyph, data = meta_nmds_water6_dna_ordered, try = 100, autotransform = FALSE, test = "permutation", permutations = how(nperm = 2000))
         Df ChiSquare      F Pr(>F)   
Model     1   0.06427 4.6372  0.002 **
Residual 28   0.38809                 

Model: cca(formula = water6_cdna_subset_nonzero ~ glyph, data = meta_nmds_water6_dna_ordered, try = 100, autotransform = FALSE)
         Df ChiSquare      F Pr(>F)    
Model     1   0.07661 5.6647  0.001 ***
Residual 28   0.37867 

##ab t6
Model: cca(formula = water6_dna_subset_nonzero ~ glyph, data = meta_nmds_water6_dna_ordered, try = 100, autotransform = FALSE)
         Df ChiSquare      F Pr(>F)    
Model     1   0.08063 6.5928  0.001 ***
Residual 26   0.31800                  

Model: cca(formula = water6_cdna_subset_nonzero ~ glyph, data = meta_nmds_water6_dna_ordered, try = 100, autotransform = FALSE)
         Df ChiSquare      F Pr(>F)    
Model     1   0.08784 6.8397  0.001 ***
Residual 26   0.33391