############combining 16S plots and omics plot
#types:
single otu
otu diversity
otu richness

gene abundance
gene diversity

genes belonging to one organism (via concoct or metaxa) with single otu

#additional data:
cell counts
glyphosate concentration


#you can further subset table or use everything for plotting. the plotting itself also provides subsetting options
sample_subset<-subset(final_table_tax_mean,habitat=="water"&nucleic_acid=="dna")
more_cell_counts<-read.csv(file.choose(),sep=";")



#specific otu
#gallaeci
gtest<-subset(final_table_tax,habitat == "water" & nucleic_acid == "dna" & grepl("Pseudomonas",variable))
meta_data_test<-subset(meta_data,habitat =="water" & nucleic_acid == "dna")
gg<-ggplot(gtest,aes(x=days))+
	geom_point(aes(y=value,colour="relative abundance"))+
	geom_point(data=meta_data_test, aes(y=shannon*5,colour="Shannon Index"))+
	geom_point(data=meta_data_test, aes(y=richness/10,colour="Species Richness"))+
	geom_line(data=more_cell_counts, aes(x=day, y=cells_ml/5000000,colour="cell counts per mL"),size=1.3)+
	geom_line(data=more_cell_counts, aes(x=day,y=glyph_mg_L,colour="glyphosate concentration"),alpha=0.8,linetype="solid", size=1.3)+
	geom_line(data=more_cell_counts, aes(x=day,y=glyph_theor,colour="glyphosate dilution"),alpha=0.5,linetype="F1", size=1.3)+
	facet_wrap(~treatment,nrow=2)+
	ggtitle("gtest")+
	stat_summary(data=gtest,fun.y="mean",geom="line",aes(y=value,colour="relative abundance"),size=1.3)+
	stat_summary(data=meta_data_test, fun.y="mean",geom="line",aes(y=shannon*5,colour="Shannon Index"),size=1.3)+
	stat_summary(data=meta_data_test, fun.y="mean",geom="line",aes(y=richness/10,colour="Species Richness"),size=1.3)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 30))+
	theme(plot.title = element_text(size=18,face="bold"))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)

##only plot diversity
divtest<-ggplot(meta_data,aes(x=days))+
	geom_point(data=meta_data,aes(y=meta_data$shannon*100,colour=nucleic_acid))+
	geom_point(data=meta_data,aes(y=meta_data$richness,shape=nucleic_acid))+
	facet_wrap(~treatment*habitat,nrow=2,ncol=2)
divtest	


#####################gene subsets
#check for glyph genes
#create gene subset
#check for bins and genes
subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("Pseudomonas",wholetax))
pseudomo_bins<-subset(kallisto_prokka_concoct2_metaxa_checkm, genus =="Pseudomonas")
pseudomo_bins<-droplevels(pseudomo_bins)
pseudomo_bins2<-pseudomo_bins[,c(1,5,6,7,9,11,13,15,17,19,21,23,25)]
library(reshape2)
molten_pseudomo<-melt(pseudomo_bins2,id=c("contig_id","ec_number","genes"))
#meta data for plotting: treatment/control and time
meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
merged_pseudomo_bins<-merge(molten_pseudomo,meta_omics,by.x="variable",by.y="row.names",all.x=T)
write.csv(merged2_pseudomo_bins,file="merged_pseudomo_bins.csv")

#pseudomo_plot<-merged_pseudomo_bins[order(merged_pseudomo_bins[,2]),]
#nrow(pseudomo_plot)
gtest<-subset(final_table_tax,habitat == "water" & nucleic_acid == "dna" & grepl("Pseudomonas",variable))

test_bins<-ggplot(merged_pseudomo_bins, aes(x = day))+
	coord_cartesian(ylim = c(0, 100))+
	geom_bar(width = 3, aes(y=value, fill=contig_id),stat = "identity",position="dodge")+
	geom_point(data=gtest,aes(x=days,y=value*5,colour="relative abundance"))+
	stat_summary(data=gtest,fun.y="mean",geom="line",aes(x=days,y=value*5,colour="relative abundance"),size=1.0)+
	geom_vline(aes(xintercept=69),linetype="dashed", size=1.2)+
	facet_wrap( ~treatment,nrow=2)
test_bins
	
	

	
	sox_plot<-subset(gallaeci_bins, grepl("sox", genes))
sox_plot<-sox_plot[order(sox_plot[,4]),]
nrow(sox_plot)

sox_genes<-ggplot(sox_plot, aes(x = day, y = value, fill=genes))+						
	geom_bar(width = 1.7, stat = "identity")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	geom_text(data=subset(sox_plot, value > 3),aes(label=contig_id),size=rel(2.7), position = position_stack(vjust = 0.5))+
	#geom_label_repel(data=subset(sox_plot, value > 50),position="stack",min.segment.length = unit(0.1, "lines"),aes(label=genes,size=1.4))+
	ggtitle("sox_genes")+
	theme(plot.title = element_text(size=18,face="bold"))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	geom_vline(aes(xintercept=70),linetype="dashed", size=1.2)+
	xlab("day")+
	ylab("rpm")+
	geom_point(aes(y=cells_ml/1000000),size=3,alpha=0.4)+
	facet_wrap( ~treatment,nrow=2)
sox_genes


#################check abundance behaviour of contigs in one bin

subset(kallisto_prokka_concoct2_metaxa_checkm, bin_gt1000 == "99")
x99_bins<-subset(kallisto_prokka_concoct2_metaxa_checkm, bin_gt1000 == "99")
nrow(x99_bins)
x99_bins<-droplevels(x99_bins)
x99_bins2<-x99_bins[,c(1,5,6,7,9,11,13,15,17,19,21,23,25)]
library(reshape2)
molten_x99<-melt(x99_bins2,id=c("contig_id","ec_number","genes"))
#meta data for plotting: treatment/control and time
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
merged_x99_bins<-merge(molten_x99,meta_omics,by.x="variable",by.y="row.names",all.x=T)
#write.csv(merged2_x99_bins,file="merged_x99_bins.csv")

#x99_plot<-merged_x99_bins[order(merged_x99_bins[,2]),]
#nrow(x99_plot)
#gtest<-subset(final_table_tax,habitat == "water" & nucleic_acid == "dna" & treatment == "glyph" & grepl("Pseudomonas",variable))

test_x99_bins<-ggplot(merged_x99_bins, aes(x = day))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=0.8)+
	#coord_cartesian(ylim = c(0, 30))+
	geom_bar(width = 3, aes(y=value, fill=contig_id),stat = "identity",position="dodge")+
	#geom_point(data=gtest,aes(x=days,y=value*5,colour="relative abundance"))+
	#stat_summary(data=gtest,fun.y="mean",geom="line",aes(x=days,y=value*5,colour="relative abundance"),size=1.0)+
	facet_wrap( ~treatment,nrow=2)+
	theme(legend.position="none")
test_x99_bins


########################################comparing two bins
subset(kallisto_prokka_concoct2_metaxa_checkm, bin_gt1000 == "84" | bin_gt1000 == "49")
x84_49_bins<-subset(kallisto_prokka_concoct2_metaxa_checkm, bin_gt1000 == "84" | bin_gt1000 == "49")
nrow(x84_49_bins)
x84_49_bins<-droplevels(x84_49_bins)
x84_49_bins2<-x84_49_bins[,c(1,5,6,7,9,11,13,15,17,19,21,23,25,26)]
library(reshape2)
molten_x84_49<-melt(x84_49_bins2,id=c("contig_id","ec_number","genes","bin_gt1000"))
#meta data for plotting: treatment/control and time
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
merged_x84_49_bins<-merge(molten_x84_49,meta_omics,by.x="variable",by.y="row.names",all.x=T)
#write.csv(merged2_x84_49_bins,file="merged_x84_49_bins.csv")

#x84_49_plot<-merged_x84_49_bins[order(merged_x84_49_bins[,2]),]
#nrow(x84_49_plot)
#gtest<-subset(final_table_tax,habitat == "water" & nucleic_acid == "dna" & treatment == "glyph" & grepl("Pseudomonas",variable))

test_x84_49_bins<-ggplot(merged_x84_49_bins, aes(x = day))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=0.8)+
	coord_cartesian(ylim = c(0, 150))+
	geom_bar(width = 3, aes(y=value, fill=contig_id),stat = "identity",position="dodge")+
	#geom_point(data=gtest,aes(x=days,y=value*5,colour="relative abundance"))+
	#stat_summary(data=gtest,fun.y="mean",geom="line",aes(x=days,y=value*5,colour="relative abundance"),size=1.0)+
	facet_wrap( ~bin_gt1000*treatment,nrow=2)+
	theme(legend.position="none")
test_x84_49_bins

#theme(legend.position="none")

#################check abundance behaviour of cut up contigs in one bin
k141_100246.0
subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("k141_100246",contig_id))
cut_up_bins<-subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("k141_100246",contig_id))
nrow(cut_up_bins)
cut_up_bins<-droplevels(cut_up_bins)
cut_up_bins2<-cut_up_bins[,c(1,5,6,7,9,11,13,15,17,19,21,23,25)]
library(reshape2)
molten_cut_up<-melt(cut_up_bins2,id=c("contig_id","ec_number","genes"))
#meta data for plotting: treatment/control and time
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
merged_cut_up_bins<-merge(molten_cut_up,meta_omics,by.x="variable",by.y="row.names",all.x=T)
#write.csv(merged2_cut_up_bins,file="merged_cut_up_bins.csv")

#cut_up_plot<-merged_cut_up_bins[order(merged_cut_up_bins[,2]),]
#nrow(cut_up_plot)
#gtest<-subset(final_table_tax,habitat == "water" & nucleic_acid == "dna" & treatment == "glyph" & grepl("Pseudomonas",variable))

test_cut_up_bins<-ggplot(merged_cut_up_bins, aes(x = day))+
	geom_vline(aes(xintercept=69),linetype="dashed", size=0.8)+
	#coord_cartesian(ylim = c(0, 30))+
	geom_bar(width = 3, aes(y=value, fill=contig_id),stat = "identity",position="dodge")+
	#geom_point(data=gtest,aes(x=days,y=value*5,colour="relative abundance"))+
	#stat_summary(data=gtest,fun.y="mean",geom="line",aes(x=days,y=value*5,colour="relative abundance"),size=1.0)+
	facet_wrap( ~treatment,nrow=2)+
	theme(legend.position="none")
test_cut_up_bins