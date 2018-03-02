###test
#change_tab$radius<-sqrt(change_tab$change/pi) #creates a new variable (radius) from the abundance data which is better suited for visualisation in bubble plots

change_tab<-read.csv(file.choose(),sep=";")
change_tab2<-change_tab
change_tab2$otu <- factor(change_tab2$otu, levels=unique(change_tab2[order(change_tab2$abu,decreasing=T),]$otu))
library(RColorBrewer)
heatmap.palette <- colorRampPalette(brewer.pal(6, 'Oranges'), space='Lab')
ggplot(change_tab2, aes(y=otu))+
	geom_tile(aes(x=nucleic_acid,fill = change))+
	geom_point(aes(x="",size=abu*5),shape=21)+
	geom_vline(aes(xintercept=6.5),linetype="dashed", size=1.2)+
	geom_vline(aes(xintercept=12.5),linetype="dashed", size=1.2)+
	scale_fill_gradientn(colours = heatmap.palette(200))+
	theme_classic()
ggsave(file="otu_changes_heatmap.jpg", width=7, height=9)	


#now make the plot nice
str(change_tab2)
change_tab2$abu<-as.factor(change_tab2$abu)
change_tab2$change<-as.factor(change_tab2$change)
str(change_tab2)

library(RColorBrewer)
heatmap.palette <- colorRampPalette(brewer.pal(5, 'Oranges'), space='Lab')
cols<-heatmap.palette(5)
#g<-guide_legend("title")
ggplot(change_tab2, aes(y=otu))+
	geom_tile(aes(x=nucleic_acid,fill = change))+
	geom_point(aes(x="",size=abu),shape=16,colour="black")+
	scale_fill_manual(values=cols,labels=c("4 - 10","10 - 20","20 - 50 ","50 - 100","100 - 150"),name="range of\nchange")+
	scale_size_discrete(labels=c("0.5-9.9","10-30"),name="maximum relative\nabundance",range = c(4, 9))+
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

# für separate DNA/cDNA plots mit zwei OTU-Listen
##only one title ("top35 otu reactions)
# somewhere DNA and cDNA

#DNA
heatmap.palette <- colorRampPalette(brewer.pal(4, 'Oranges'), space='Lab')
cols<-heatmap.palette(4)
cols = c("#FEECDE","#FDD0A1", "#FDAD6B", "#FD8C3B")
change_tab2_dna<-subset(change_tab2, nucleic_acid=="dna")	
change_tab2_dna$otu <- factor(change_tab2_dna$otu, levels=unique(change_tab2_dna[order(change_tab2_dna$change,decreasing=F),]$otu))
ggplot(change_tab2_dna, aes(y=otu))+
	geom_tile(aes(x=nucleic_acid,fill = change))+
	geom_point(aes(x="",size=abu),shape=16,colour="black")+
	scale_fill_manual(values=cols,labels=c("4 - 10","10 - 20","20 - 50 ","50 - 100","100 - 150"),name="range of\nchange")+
	scale_size_discrete(labels=c("0.005-0.01","0.01-0.1","0.1-0.5","0.5-3","3-15"),name="maximum relative\nabundance",range = c(2, 8))+
	theme_classic()+
	xlab("")+
	ylab("OTU")+
	ggtitle("Changes in OTU abundance")+
	coord_fixed(ratio = 0.4)+
	scale_x_discrete(breaks=c("","dna"),labels=c("max_rel\n_abu","dna\nchange"))+
	theme(plot.title = element_text(size=18,face="bold"))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.text.y = element_text(face = "italic"))+
	#theme(axis.text.x = element_blank())+
	theme(axis.ticks = element_blank())

	
# für separate DNA/cDNA plots mit zwei OTU-Listen
#cDNA
cols = c("#FEECDE", "#FDD0A1", "#FDAD6B", "#FD8C3B", "#E5540C")
change_tab2_cdna<-subset(change_tab2, nucleic_acid=="cdna")	
change_tab2_cdna$otu <- factor(change_tab2_cdna$otu, levels=unique(change_tab2_cdna[order(change_tab2_cdna$change,decreasing=F),]$otu))
ggplot(change_tab2_cdna, aes(y=otu))+
	geom_tile(aes(x=nucleic_acid,fill = change))+
	geom_point(aes(x="",size=abu),shape=16,colour="black")+
	scale_fill_manual(values=cols,labels=c("7 - 10","10 - 20","20 - 50 ","50 - 100","100 - 150"),name="range of\nchange")+
	scale_size_discrete(labels=c("0.005-0.01","0.01-0.1","0.1-0.5","0.5-3","3-15"),name="maximum relative\nabundance",range = c(2, 8))+
	theme_classic()+
	xlab("")+
	ylab("OTU")+
	ggtitle("Changes in OTU abundance")+
	coord_fixed(ratio = 0.4)+
	scale_x_discrete(breaks=c("cdna",""),labels=c("cdna\nchange","max_rel\n_abu"))+
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
cols = c("#FEECDE","#FDD0A1", "#FDAD6B", "#FD8C3B")
cols = c("#FEECDE", "#FDD0A1", "#FDAD6B", "#FD8C3B", "#E5540C")
cols = c("#FEECDE", "#FDD0A1", "#FDAD6B", "#FD8C3B", "#E5540C", "#A53502")	
ggplot(df,aes(x = Var1,y = Var2, fill = z)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = c("white", "green", "red"), values = c(0,0.1,1))
  
   You can create cuts using tmp <- cut(change_tab2$change, seq(0, .5, .1)). 
 scale_fill_gradient(breaks = seq(0, .5, .1), 
   labels = paste(levels(tmp), table(tmp))). – 
	
	http://ggplot2.tidyverse.org/reference/guides.html
	
	For the following bar graph:

x <- ggplot(foo, aes(x=variety, y=percent)) + geom_bar()

The following italicizes all of my x-axis text:

opts(axis.text.y=theme_text(face='italic'))


ggplot(CO2, aes(y=uptake,x=Type, group=Type))+
  geom_point()+
  scale_x_discrete("Location", labels=expression(Quebec, italic(Mississippi)))
  breaks=c("Quebec", "Mississippi")
	
	

	