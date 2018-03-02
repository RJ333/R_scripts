###test
#richness_change$radius<-sqrt(richness_change$change/pi) #creates a new variable (radius) from the abundance data which is better suited for visualisation in bubble plots

richness_change<-read.csv(file.choose(),sep=";")
head(richness_change)
str(richness_change)

library(RColorBrewer)
heatmap.palette1 <- colorRampPalette(brewer.pal(9, 'Oranges'), space='Lab')
heatmap.palette2 <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')

cols_orange<-heatmap.palette1(177)
cols_blue<-heatmap.palette2(44)

ggplot(richness_change, aes(y=gene))+
	geom_tile(data=subset(richness_change,pos>2),aes(x=pos,fill=percent))+
	geom_tile(data=subset(richness_change,pos==2),aes(x=pos,fill=percent*0.6))+
	geom_point(aes(x=pos, size= max),shape=16)

str(richness_change)	
richness_change$diffmax5<-as.factor(richness_change$diffmax5)
richness_change$diff56<-as.factor(richness_change$diff56)		
richness_change$diff67<-as.factor(richness_change$diff67)		
richness_change$diffmax5<-as.factor(richness_change$diffmax5)		
str(richness_change)
				
ggplot(richness_change, aes(y=gene))+
	geom_tile(richness_change,aes(fill=percent))+
	geom_tile(richness_change,aes(x=pos,fill=diffmax5))+
	geom_tile(richness_change,aes(x=pos,fill=diff56))+
	geom_tile(richness_change,aes(x=pos,fill=diff67))+
	geom_point(aes(x=pos, size= max),shape=16)	
	
	
	
	
library(RColorBrewer)
heatmap.palette <- colorRampPalette(brewer.pal(9, 'Oranges'), space='Lab')
cols<-heatmap.palette(16)
#g<-guide_legend("title")
+
	
	
	scale_fill_manual(values=cols,
	labels=c("90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105"),name="range of\nchange")+
	scale_size_discrete(labels=c("0.01-0.49","0.5-9.9","10-30"),name="maximum relative\nabundance",range = c(4, 9))+
	theme_classic()+
	xlab("")+
	ylab("OTU")+
	ggtitle("Changes in OTU abundance")+
	coord_fixed(ratio = 0.4)+
	scale_x_discrete(breaks=c("","cdna","dna"),labels=c("max_rel\n_abu","cdna\nchange","dna\nchange"))+
	theme(plot.title = element_text(size=18,face="bold"))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.text.y = element_text(face = "italic"))+
	#theme(axis.text.x = element_blank())+
	theme(axis.ticks = element_blank())
	
	
#orange gradient for heatmap
(colours = c("#FEECDE", "#FDD0A1", "#FDAD6B", "#FD8C3B", "#E5540C", "#A53502"))
	
ggplot(df,aes(x = Var1,y = Var2, fill = z)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = c("white", "green", "red"), values = c(0,0.1,1))

  #creating data intervals  
   You can create cuts using tmp <- cut(richness_change2$change, seq(0, .5, .1)). 
 scale_fill_gradient(breaks = seq(0, .5, .1), 
   labels = paste(levels(tmp), table(tmp))). – 
	
	http://ggplot2.tidyverse.org/reference/guides.html
	
	For the following bar graph:

x <- ggplot(foo, aes(x=variety, y=percent)) + geom_bar()

#italics
The following italicizes all of my x-axis text:

opts(axis.text.y=theme_text(face='italic'))


ggplot(CO2, aes(y=uptake,x=Type, group=Type))+
  geom_point()+
  scale_x_discrete("Location", labels=expression(Quebec, italic(Mississippi)))
  breaks=c("Quebec", "Mississippi")
	
	

	