library(vegan)
library(gridExtra)
library(ggplot2)
library(grid)

############# clean final version ################################################################
############################ für DNA only, Y Achse gespiegelt ########################################
dna_cca_5_plot<-ggplot() + 
  geom_path(data=subset(path_both5,nucleic_acid=="dna" & treatment == "glyph"),aes(x=CCA1,y=CA1,group=treatment),colour="black",size=2,alpha=1)+
  geom_path(data=subset(path_both5,nucleic_acid=="dna"& treatment == "control"),aes(x=CCA1,y=CA1,group=treatment),colour="grey70",size=2,alpha=1)+
  geom_point(data=subset(combined_cca5,nucleic_acid=="dna"),aes(x=CCA1,y=CA1,shape=treatment),size=4,colour="black",fill="white",alpha=1)+
  coord_cartesian(ylim = c(-2.3, 2.3),xlim=c(-1.5,3))+
  scale_y_reverse()+
  theme_bw()+
	theme(axis.text=element_text(size=30))+
	theme(legend.position="none")+
	theme(axis.title=element_blank())+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))
  
############################ für cDNA only, nicht gespiegelt ########################################  
 
 cdna_cca_5_plot<-ggplot() + 
  geom_path(data=subset(path_both5,nucleic_acid=="cdna" & treatment == "glyph"),aes(x=CCA1,y=CA1,group=treatment),colour="black",size = 2, alpha = 1)+
  geom_path(data=subset(path_both5,nucleic_acid=="cdna"& treatment == "control"),aes(x=CCA1,y=CA1,group=treatment),colour="grey70",size = 2, alpha = 1)+
  geom_point(data=subset(combined_cca5,nucleic_acid=="cdna"),aes(x=CCA1,y=CA1,shape=treatment),size=4,colour="black",fill="white",alpha = 1)+
  coord_cartesian(ylim = c(-2.3, 2.3),xlim=c(-1.5,3))+
  theme_bw()+
	theme(axis.text=element_text(size=30))+
	theme(legend.position="none")+
	theme(axis.title=element_blank())+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))
 
#########################################################################zusammenführen mit gridExtra
grid.arrange(dna_cca_5_plot, cdna_cca_5_plot, nrow = 1, bottom=textGrob("CCA1",gp=gpar(fontsize=35)),left=textGrob("CA1",rot=90,gp=gpar(fontsize=35)))


ggsave(file="Fig_04_cca_water5_paper_arrows.png", width = 20, height = 8.75, grid.arrange(dna_cca_5_plot, cdna_cca_5_plot, nrow = 1, bottom = textGrob("CCA1", gp = gpar(fontsize = 35)), left = textGrob("CA1", rot = 90, gp = gpar(fontsize = 35))))

################## end clean final version #########################################################################################################################

#zero counts entfernen
water_cdna_subset2<-read.csv(file.choose(),sep=";",row.names=1) #11_table_for_vegan_wo_start_cdna, meta in sample name enthalten, nur wasser
separation<-colSums(water_cdna_subset2)!=0
water_cdna_subset_nonzero<-water_cdna_subset2[,separation]
ncol(water_cdna_subset2)
ncol(water_cdna_subset_nonzero)
names(water_cdna_subset_nonzero)
head(water_cdna_subset_nonzero)

meta_nmds_water5_dna_ordered<-subset(meta_nmds_water5_dna_ordered,time>5)
water5_cdna_subset_nonzero<-read.csv(file.choose(),sep=";",row.names=1)
water5_dna_subset_nonzero<-read.csv(file.choose(),sep=";",row.names=1)

head(meta_nmds_water5_dna_ordered)
head(meta_nmds_water5_dna_ordered)
head(water5_cdna_subset_nonzero)

######################nur cDNA
#CCA
cca_water5_cdna_meta<-cca(water5_cdna_subset_nonzero~glyph,meta_nmds_water5_dna_ordered,try=100,autotransform=FALSE)
plot(cca_water5_cdna_meta,type="p", display="sites")
text(cca_water5_cdna_meta,labels=row.names(water5_cdna_subset_nonzero),cex=0.7)

#quality check
eigenvals(cca_water5_cdna_meta)
summary(cca_water5_cdna_meta)

#test mit autotransform (looks exactly the same)
#cca_water5_cdna_meta_true<-cca(water5_cdna_subset_nonzero~glyph,meta_nmds_water5_dna_ordered,try=100,autotransform=TRUE)
#plot(cca_water5_cdna_meta_true,type="p", display="sites")
#text(cca_water5_cdna_meta_true,labels=row.names(water5_cdna_subset_nonzero),cex=0.7)

#extract data
cca_water5_cdna_coord<-as.data.frame(scores(cca_water5_cdna_meta)$sites)
str(cca_water5_cdna_coord)


#neue meta file notwendig, da neue row.names
#meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
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
#quality check
eigenvals(cca_water5_dna_meta)
summary(cca_water5_dna_meta)
#extract data
cca_water5_dna_coord<-as.data.frame(scores(cca_water5_dna_meta)$sites)
str(cca_water5_dna_coord)

#neue meta file notwendig, da neue row.names
#meta_file_vegan<-read.csv(file.choose(),sep=";",row.names=1)
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
combined_cca5<-read.csv(file.choose(),sep=";",row.names=1)

#plotting combined

path_both5<-combined_cca5
path_both5$time<-as.numeric(as.character(path_both5$time))
path_both5$days<-as.numeric(as.character(path_both5$days))
path_both5$new_day<-as.numeric(as.character(path_both5$new_day))
path_both5<-path_both5[order(path_both5$days),]

##gefärbte punkte
ggplot() + 
  geom_point(data=subset(combined_cca5,group=="normal"),aes(x=CCA1,y=CA1,shape=treatment,colour="normal"),size=2.5)+
  geom_point(data=subset(combined_cca5,group=="bloom"),aes(x=CCA1,y=CA1,shape=treatment,colour="bloom"),size=4)+
  geom_path(data=path_both5,aes(x=CCA1,y=CA1,group=treatment,colour=as.numeric(as.character(new_day))),size=1)+
  geom_text(data=combined_cca5,aes(x=CCA1,y=CA1,label=new_day,size=5,vjust=-0.7))+
  coord_equal() +
  facet_wrap(~nucleic_acid)

##gefärbte linien ab zeitpunkt6
path_both5<-read.csv(file.choose(),sep=";",row.names=1)
path_both5<-path_both5[with(path_both5, order(nucleic_acid)), ]
path_both5$nucleic_acid2<-factor(path_both5$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rDNA amplicons"))
combined_cca5$nucleic_acid2<-factor(combined_cca5$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rDNA amplicons"))
ggplot() + 
  geom_path(data=path_both5,aes(x=CCA1,y=CA1,group=treatment),colour="grey",size=1.5)+
  geom_point(data=combined_cca5,aes(x=CCA1,y=CA1,shape=treatment),size=3,colour="black",fill="white")+
  geom_text(data=combined_cca5,aes(x=CCA1,y=CA1,label=days-69,vjust=-1),size=3)+
  #scale_colour_discrete(limits=levels(path_both5$group))+
  #scale_colour_manual(values=c("normal"="grey20","bloom"="grey60"),
						#	name="Community succession",
							#breaks=c("normal","bloom"),
							#labels=c("normal",'"changed"'))+
  scale_shape_manual(values=c("control"=21,"glyph"=22),
							name="Microcosm",
							breaks=c("control","glyph"),
							labels=c("Control","Treatment"))+
  coord_equal() +
  facet_wrap(~nucleic_acid2)+
  theme_bw()
ggsave(file="cca_water5_paper_test.png", width=12, height=7)
	facet_wrap(~material2,ncol=2,labeller=label_value)

############################################################################################################################################################################################
	
###############################################################für DNA only, y Achse gespiegelt
path_both5<-read.csv(file.choose(),sep=";",row.names=1)
path_both5<-path_both5[with(path_both5, order(nucleic_acid)), ]
path_both5$nucleic_acid2<-factor(path_both5$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rRNA gene amplicons"))
combined_cca5$nucleic_acid2<-factor(combined_cca5$nucleic_acid,labels=c("Communities by 16S rRNA amplicons","Communities by 16S rRNA gene amplicons"))


dna_cca_5_plot<-ggplot() + 
  geom_path(data=subset(path_both5,nucleic_acid=="dna"),aes(x=CCA1,y=CA1,group=treatment),colour="grey",size=1.5)+
  geom_point(data=subset(combined_cca5,nucleic_acid=="dna"),aes(x=CCA1,y=CA1,shape=treatment),size=3,colour="black",fill="white",alpha=0.8)+
  #geom_text(data=subset(combined_cca5,nucleic_acid=="dna"&treatment=="glyph"),fontface="bold",aes(x=CCA1,y=CA1,label=days-69,hjust=-0.5,vjust=0.7),size=4)+
  #geom_text(data=subset(combined_cca5,nucleic_acid=="dna"&treatment=="control"),fontface="bold",aes(x=CCA1,y=CA1,label=days-69,hjust=1.3,vjust=0.3),size=4)+
  #scale_colour_discrete(limits=levels(path_both5$group))+
  #scale_colour_manual(values=c("normal"="grey20","bloom"="grey60"),
						#	name="Community succession",
							#breaks=c("normal","bloom"),
							#labels=c("normal",'"changed"'))+
  scale_shape_manual(values=c("control"=15,"glyph"=16),
							name="Microcosm",
							breaks=c("control","glyph"),
							labels=c("Control","Treatment"))+
  coord_cartesian(ylim = c(-2, 2.5),xlim=c(-2,3.5))+
  #xlim(-1.5,3.5)+
  #ylim(-2.5,2)+
  scale_y_reverse()+
  #coord_equal() +
  theme_bw()+
  	#theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	#theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=13))+
	theme(legend.position="none")+
	theme(axis.title=element_blank())+
	#xlab("CCA1")+
	#ylab('CA1')+
	#theme(legend.title=element_text(size=13,face="bold"))+
	#theme(legend.text=element_text(size=11))+
  theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
  theme(panel.grid.minor=element_line(colour = NA, size = 0.5))
  #annotate("segment", x = 1.6, xend = 1.6, y = 0.1, yend = 0.30, colour = "black", size=1.2, alpha=1, arrow=arrow())+
  #annotate("segment", x = 1.5, xend = 0, y = 2, yend = 0.4, colour = "black", size=3, alpha=1, arrow=arrow())+
  #annotate("text", x = 0.95, y = 1.2, label = "average succession" , color="black", size=5 , angle=-53, fontface="bold")+
  #annotate("text", x = 1.6, y = 0, label = "glyphosate pulse" , color="black", size=4 , angle=0, fontface="bold")
  
####################################################################für cDNA only, nicht gespiegelt  
 combined_cca5$jit<-with(combined_cca5, ifelse(days == 79 | days == 133 | days == 91, 1,0))

 
 cdna_cca_5_plot<-ggplot() + 
  geom_path(data=subset(path_both5,nucleic_acid=="cdna"),aes(x=CCA1,y=CA1,group=treatment),colour="grey",size=1.5)+
  geom_point(data=subset(combined_cca5,nucleic_acid=="cdna"),aes(x=CCA1,y=CA1,shape=treatment),size=3,colour="black",fill="white",alpha=0.8)+
 # geom_text(data=subset(combined_cca5,nucleic_acid=="cdna"&treatment=="glyph"),fontface="bold",aes(x=CCA1,y=CA1,label=days-69,hjust=-0.2,vjust=-0.5),size=4)+
  #geom_text(data=subset(combined_cca5,nucleic_acid=="cdna"&treatment=="control"),fontface="bold",aes(x=CCA1,y=CA1,label=days-69,hjust=-0.5,vjust=0.5),size=4)+
  #scale_colour_discrete(limits=levels(path_both5$group))+
  #scale_colour_manual(values=c("normal"="grey20","bloom"="grey60"),
						#	name="Community succession",
							#breaks=c("normal","bloom"),
							#labels=c("normal",'"changed"'))+
  scale_shape_manual(values=c("control"=15,"glyph"=16),
							name="Microcosm",
							breaks=c("control","glyph"),
							labels=c("Control","Treatment"))+
  #scale_y_reverse()+
  coord_cartesian(ylim = c(-2, 2.5),xlim=c(-2,3.5))+
  #xlim(-1.5,3.5)+
  #ylim(-2,2)+
  #coord_equal() +
  theme_bw()+
  #theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
  #theme(axis.title = element_text(size=12,face="bold"))+
  #theme(axis.title.y = element_text(angle=90,vjust=0.5))+
  theme(axis.text=element_text(size=13))+
  theme(legend.position="none")+
  theme(axis.title=element_blank())+
  #xlab("CCA1")+
  #ylab('CA1')+
  #theme(legend.title=element_text(size=13,face="bold"))+
  #theme(legend.text=element_text(size=11))+
  theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
  theme(panel.grid.minor=element_line(colour = NA, size = 0.5))
  #annotate("segment", x = 1.4, xend = 1.8, y = 0.4, yend = 0.40, colour = "black", size=1.2, alpha=1, arrow=arrow())+
  #annotate("segment", x = 1.0, xend = 0, y = -1.5, yend = 0.3, colour = "black", size=3, alpha=1, arrow=arrow())+
  #annotate("text", x = 0.8, y = -0.5, label = "average succession" , color="black", size=5 , angle=-65, fontface="bold")+
  #annotate("text", x = 1.1, y = 0.5, label = "glyphosate pulse" , color="black", size=4 , angle=0, fontface="bold")
 
#########################################################################zusammenführen mit gridExtra
library(gridExtra)
grid.arrange(dna_cca_5_plot, cdna_cca_5_plot, nrow = 1, bottom=textGrob("CCA1"),left=textGrob("CA1",rot=90))

grid.arrange(dna_cca_5_plot, cdna_cca_5_plot, nrow = 1)
ggsave(file="cca_water5_paper_arrows.png", width=20, height=8.75,arrangeGrob(dna_cca_5_plot, cdna_cca_5_plot,nrow=1))
