meta_omics<-read.csv(file.choose(),sep=";")
meta_omics_glyph<-subset(meta_omics,treatment=="glyph")
diversity_16S<-read.csv(file.choose(),sep=";")
div_dna_water<-subset(diversity_16S,habitat=="water" & nucleic_acid == "dna" & days >= 60)
div_dna_water_glyph<-subset(diversity_16S,habitat=="water" & nucleic_acid == "dna" & days >= 60 & treatment == "glyph")
more_cell_counts<-read.csv(file.choose(),sep=";")
more_cell_counts_60_glyph<-subset(more_cell_counts,treatment=="glyph" & day >= 60)

all_unique_contigs_barplot<-ggplot(meta_omics_glyph, aes(x = day))+						
	geom_bar(data=subset(meta_omics_glyph,group=1),aes(y=unique_contigs, colour="unique\ncontigs\n"), stat="identity", width=2)+
	geom_bar(data=subset(meta_omics_glyph,group=2),aes(y=unique_contigs, colour="unique\ncontigs\n",fill=line), position="stack",stat="identity", width=2)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#ylim(0,1000000)+
	ggtitle("contig, gene and species richness")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	geom_vline(aes(xintercept=70),linetype="dashed", size=1.2)+
	xlab("day")+
	ylab("relative units")+
	geom_point(data=subset(more_cell_counts_60_glyph,group!="bloom"),aes(y=cells_ml/900,colour="cell_counts\n"),size=5,alpha=1)+
	geom_point(data=subset(more_cell_counts_60_glyph,group=="bloom"),aes(y=cells_ml/900,colour="cell_counts\nbloom"),size=5,alpha=1)+
	geom_line(aes(y=gene_amount/2.5,colour="gene_amount\n"),size=1.2)+
	#geom_line(aes(y=unique_gene,colour="unique_gene\n"))+
	#geom_line(aes(y=total_gene_richness*10,colour="total_gene_\nrichness\n"))+
	stat_summary(data=subset(div_dna_water_glyph,group=="standard"),aes(x=days,y=richness*1000,colour="species_richness\n"),fun.y="mean", geom="line", size=1.2)+
	stat_summary(data=subset(div_dna_water_glyph,group!="standard"),aes(x=days,y=richness*1000,colour="species_richness\nbloom"),fun.y="mean", geom="line", size=1.2)
	#geom_point(data=Gallaecimonas_water_dna_plot,aes(x=days, y=value*30000,colour="Gallaecimonas sp.\nDNA"))+
	#stat_summary(data=Gallaecimonas_water_dna_plot,aes(x=days, y=value*30000,colour="Gallaecimonas sp.\nDNA"),fun.y="mean", geom="line", size=1.5)+
	#theme(legend.position="none")+
all_unique_contigs_barplot
 ggsave(file="contig_and_gene_and_species_richness_paper.jpg", width=14, height=8)
 
#total gene richness: number of different gene classes per sample --> almost no change
#unique genes: very similar to unique contigs what is unique genes????
#unique contigs: number of unique contigs per sample, 1000.1 and 1000.2 parts count as 1 contig
#gene amount: number of gene classes * number of instances

#glyphosat noch miteinabuen? oder aus ersten plot erwähnen? zuviel für letzten plot

#colour for cell counts differ during bloom!
#colour for species richness should do the same!
	
	
	