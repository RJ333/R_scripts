#write.csv(abs_abu_substract_tax,file="abs_abu_substract_tax.csv")
###########otu plots with glyphosate addition
#rhizobium
more_cell_counts_44_glyph<-subset(more_cell_counts_44,treatment =="glyph")
Rhizobium_subtracted_plot<-subset(abs_abu_substract_tax, genus == "Rhizobium" & days > 40)
species_title_Rhizobium<-expression(paste(,italic("Rhizobium")," sp."))
ggplot(data=Rhizobium_subtracted_plot, aes(x=days,y=new_value))+ 
	#coord_cartesian(ylim = c(0, 1))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	#geom_point(data=more_cell_counts_44_glyph, aes(x= day, y=cells_ml,colour="cell_counts"),shape=21,fill="green",size=3,alpha=0.5)+
	geom_line(aes(group=nucleic_acid, colour=nucleic_acid))+
	scale_colour_manual(values=c("blue","green","red"),
						name="Nucleic Acid",
						breaks=c("cdna","dna","cells_mL"),
						labels=c("16S-rRNA","16S-rDNA","cell_counts"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle(species_title_Rhizobium)+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	ylab("absolute\ncorrected\nabundance")
ggsave(file="absolute_genus_Rhizobium_subtract.jpg", width=14, height=8)