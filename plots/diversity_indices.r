####plot diversity indices
diversity_16S<-read.csv(file.choose(),sep=";",row.names=1)

test_div_plot<-subset(diversity_16S, habitat == "water")
##shannon in water column
ggplot(data=test_div_plot, aes(x=days,colour=nucleic_acid))+ 
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	geom_point(aes(y=shannon))+
	geom_point(aes(y=max_shannon))+
	geom_point(aes(y=evenness*3))+
	stat_summary(aes(y=shannon,colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	stat_summary(aes(y=max_shannon,colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	stat_summary(aes(y=evenness*3,colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	geom_point(data=more_cell_counts_44, aes(x= day, y=cells_ml/(15000000),colour="cell_counts"),shape=21,fill="green",size=3,alpha=1)+
	facet_wrap(~treatment,nrow=2)
ggsave(file="shannon_pielou_water.jpg", width=14, height=8)
	
	+
	scale_colour_manual(values=c("blue","green","red"),
						name="Nucleic Acid",
						breaks=c("cdna","dna","cells_mL"),
						labels=c("16S-rRNA","16S-rDNA","cell_counts"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle(species_title_test)+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	ylab("relative\nabundance")
ggsave(file="genus_test_div.jpg", width=14, height=8)

##richness in water column
ggplot(data=test_div_plot, aes(x=days,colour=nucleic_acid))+ 
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	geom_point(aes(y=richness))+
	stat_summary(aes(y=richness,colour=nucleic_acid),fun.y="mean", geom="line", size=1.5)+
	geom_point(data=more_cell_counts_44, aes(x= day, y=cells_ml/(300000),colour="cell_counts"),shape=21,fill="green",size=3,alpha=1)+
	facet_wrap(~treatment,nrow=2)
ggsave(file="richness_water.jpg", width=14, height=8)