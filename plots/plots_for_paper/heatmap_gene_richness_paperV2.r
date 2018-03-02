##heatmap 
library(ggplot2)
library(reshape2)
richness_change<-read.csv(file.choose(),sep=";") #selected_richness_genes_ggplotxx
head(richness_change)
richness_change2<-melt(richness_change,id=c("gene","instances_0"))
head(richness_change2)
str(richness_change2)

ggplot(richness_change2, aes(x=variable,y=gene))+
	geom_tile(aes(fill=value))+
	geom_point(aes(x="", size=instances_0),shape=16)	
	
#one colorscale for negative and one for positive changes would be cool
#colour for treatment/control separation
#include "phnM-similars"
#radius for instances_0?
	
	
