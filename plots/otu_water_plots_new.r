###########otu plots with glyphosate addition
#Gallaecimonas	
Gallaecimonas_water_plot<-subset(final_table_tax, habitat == "water" & genus == "Gallaecimonas")
species_title_Gallaecimonas<-expression(paste(,italic("Gallaecimonas")," sp."))
ggplot(data=Gallaecimonas_water_plot, aes(x=days,y=value))+ 
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	geom_point(data=more_cell_counts_44, aes(x= day, y=cells_ml*max(Gallaecimonas_water_plot$value)/(450000000),colour="cell_counts"),shape=21,fill="green",size=3,alpha=0.5)+
	geom_point(aes(group=nucleic_acid, colour=nucleic_acid))+
	stat_summary(aes(colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	facet_wrap(~treatment,nrow=2)+
	scale_colour_manual(values=c("blue","green","red"),
						name="Nucleic Acid",
						breaks=c("cdna","dna","cells_mL"),
						labels=c("16S-rRNA","16S-rDNA","cell_counts"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle(species_title_Gallaecimonas)+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	ylab("relative\nabundance")
#ggsave(file="genus_Gallaecimonas_water.jpg", width=14, height=8)
	
	
	
