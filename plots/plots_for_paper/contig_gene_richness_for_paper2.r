meta_omics<-read.csv(file.choose(),sep=";")
meta_omics_glyph<-subset(meta_omics,treatment=="glyph"&day>60)
diversity_16S<-read.csv(file.choose(),sep=";")
div_dna_water<-subset(diversity_16S,habitat=="water" & nucleic_acid == "dna" & days >= 60)
div_dna_water_glyph<-subset(diversity_16S,habitat=="water" & nucleic_acid == "dna" & days >= 60 & treatment == "glyph")
more_cell_counts<-read.csv(file.choose(),sep=";")
more_cell_counts_60_glyph<-subset(more_cell_counts,treatment=="glyph" & day >= 60)

ggplot(meta_omics_glyph, aes(x = new_day))+						
	geom_line(data=subset(meta_omics_glyph,sequence=="metagenome"),aes(y=unique_contigs, colour="unique\ncontigs\n"))+
	#geom_line(data=subset(meta_omics_glyph,sequence=="metagenome"),aes(y=unique_gene,colour="unique_gene\n"))+ #das gleich wie unique contigs
	geom_line(data=subset(meta_omics_glyph,sequence=="metagenome"),aes(y=total_gene_richness*30,colour="total_gene_\nrichness\n"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	geom_line(data=subset(meta_omics_glyph,sequence=="metagenome"),aes(y=gene_amount/1.5,colour="gene_amount\n"))+
	geom_line(data=subset(meta_omics_glyph,sequence=="amplicon"&nucleic_acid=="cdna"),aes(y=evenness*100000,colour="evenness\ncDNA\n"))+
	geom_line(data=subset(meta_omics_glyph,sequence=="amplicon"&nucleic_acid=="cdna"),aes(y=species.richness*400,colour="species_richness\ncDNA\n"))+
	geom_line(data=subset(meta_omics_glyph,sequence=="amplicon"&nucleic_acid=="dna"),aes(y=evenness*100000,colour="evenness\nDNA\n"))+
	geom_line(data=subset(meta_omics_glyph,sequence=="amplicon"&nucleic_acid=="cdna"),aes(y=shannon*20000,colour="Shannon\ncDNA\n"))+
	geom_vline(aes(xintercept=0),linetype="dashed", size=1.2,color="blue")+
	ggtitle("contig, gene and species richness")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("day")+
	ylab("relative units")+
	facet_wrap(~sequence,nrow=2)
 ggsave(file="contig_and_gene_and_species_richness_paper.jpg", width=14, height=8)
 
#total gene richness: number of different gene classes per sample --> almost no change
#unique genes: very similar to unique contigs what is unique genes????
#unique contigs: number of unique contigs per sample, 1000.1 and 1000.2 parts count as 1 contig
#gene amount: number of gene classes * number of instances

#glyphosat noch miteinabuen? oder aus ersten plot erwähnen? zuviel für letzten plot

#colour for cell counts differ during bloom!
#colour for species richness should do the same!
	
	
	