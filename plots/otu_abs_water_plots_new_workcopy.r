#write.csv(abs_abu_molten_tax,file="abs_abu_molten_tax.csv")
abs_abu_molten_tax<-read.csv(file.choose(),sep=";")
###########otu plots with glyphosate addition
#Marinobacter
Pseudoalteromonas_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & genus == 'Pseudoalteromonas' & days > 40)
#Pseudoalteromonas_water_plot<-subset(Pseudoalteromonas_water_plot, treatment == "glyph")
#Pseudoalteromonas_water_plot<-subset(Pseudoalteromonas_water_plot, days >= 44)
species_title_Pseudoalteromonas<-expression(paste(,italic("Pseudoalteromonas")," sp."))
ggplot(data=Pseudoalteromonas_abs_water_plot, aes(x=days,y=value))+ 
	#coord_cartesian(ylim = c(0, 10000))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	#geom_point(data=more_cell_counts_44, aes(x= day, y=cells_ml,colour="cell_counts"),shape=21,fill="green",size=3,alpha=0.5)+
	geom_point(aes(group=nucleic_acid, colour=nucleic_acid))+
	stat_summary(aes(colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	facet_wrap(~treatment,nrow=2)+
	scale_colour_manual(values=c("blue","green","red"),
						name="Nucleic Acid",
						breaks=c("cdna","dna","cells_mL"),
						labels=c("16S-rRNA","16S-rDNA","cell_counts"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle(species_title_Pseudoalteromonas)+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	ylab("absolute\nabundance")
ggsave(file="absolute_genus_Pseudoalteromonas_water.jpg", width=14, height=8)

 [1] "_uncultured_Alteromonadales_Gammaproteobacteria_Proteobacteria"                                     "Alteromonas_Alteromonadaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                   
 [3] "Colwellia_Colwelliaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                         "Gallaecimonas_Unknown-Family_Gammaproteobacteria-Incertae-Sedis_Gammaproteobacteria_Proteobacteria"
 [5] "Idiomarina_Idiomarinaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                       "Marinobacter_Alteromonadaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                  
 [7] "Paraglaciecola_Alteromonadaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                 "Pseudoalteromonas_Pseudoalteromonadaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"       
 [9] "Salinimonas_Alteromonadaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                    "Shewanella_Shewanellaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                      
[11] "uncultured_Colwelliaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"                        "uncultured_Pseudoalteromonadaceae_Alteromonadales_Gammaproteobacteria_Proteobacteria"  

#uncultured
uncultured_Alteromonadales_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & genus == "uncultured" & order =="Alteromonadales" & days > 40)
#uncultured_Alteromonadales_water_plot<-subset(uncultured_Alteromonadales_water_plot, treatment == "glyph")
#uncultured_Alteromonadales_water_plot<-subset(uncultured_Alteromonadales_water_plot, days >= 44)
species_title_uncultured_Alteromonadales<-expression(paste(,italic("uncultured_Alteromonadales")," sp."))
ggplot(data=uncultured_Alteromonadales_abs_water_plot, aes(x=days,y=value))+ 
	#coord_cartesian(ylim = c(0, 1))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	#geom_point(data=more_cell_counts_44, aes(x= day, y=cells_ml,colour="cell_counts"),shape=21,fill="green",size=3,alpha=0.5)+
	geom_point(aes(group=nucleic_acid, colour=nucleic_acid))+
	stat_summary(aes(colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	facet_wrap(~treatment,nrow=2)+
	scale_colour_manual(values=c("blue","green","red"),
						name="Nucleic Acid",
						breaks=c("cdna","dna","cells_mL"),
						labels=c("16S-rRNA","16S-rDNA","cell_counts"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle(species_title_uncultured_Alteromonadales)+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	ylab("absolute\nabundance")
ggsave(file="absolute_genus_uncultured_Alteromonadales_water.jpg", width=14, height=8)

#
Hot-Creek_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & order =="Hot-Creek-32" & days > 40)
#Hot-Creek_water_plot<-subset(Hot-Creek_water_plot, treatment == "glyph")
#Hot-Creek_water_plot<-subset(Hot-Creek_water_plot, days >= 44)
species_title_Hot-Creek<-expression(paste(,italic("Hot-Creek")," sp."))
ggplot(data=Hot-Creek_abs_water_plot, aes(x=days,y=value))+ 
	#coord_cartesian(ylim = c(0, 1))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	#geom_point(data=more_cell_counts_44, aes(x= day, y=cells_ml,colour="cell_counts"),shape=21,fill="green",size=3,alpha=0.5)+
	geom_point(aes(group=nucleic_acid, colour=nucleic_acid))+
	stat_summary(aes(colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	facet_wrap(~treatment,nrow=2)+
	scale_colour_manual(values=c("blue","green","red"),
						name="Nucleic Acid",
						breaks=c("cdna","dna","cells_mL"),
						labels=c("16S-rRNA","16S-rDNA","cell_counts"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle(species_title_Hot-Creek)+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	ylab("absolute\nabundance")
ggsave(file="absolute_Hot-Creek_water.jpg", width=14, height=8)