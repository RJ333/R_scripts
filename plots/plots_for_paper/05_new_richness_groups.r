richness3<-read.csv(file.choose(),row.names=1,sep=";")
trichness3<-t(richness3)
#meta_omics<-read.csv(file.choose(),sep=";")
trichness_meta<-merge(trichness3,meta_omics,by.x="row.names",by.y="X",all.x=TRUE)
row.names(trichness_meta)<-trichness_meta$Row.names
trichness_meta<-trichness_meta[,-c(34:43)]
trichness_meta<-trichness_meta[-8,] #row.name "group" bringt hier noch nix
library(reshape2)
richm<-melt(trichness_meta,id=c("Row.names","time","day"))
write.csv(richm,file="richm.csv")
#in excel 4 gruppen hinzugefügt
richm2<-read.csv(file.choose(),row.names=1,sep=";")

#excluding rpoB
richm3<-subset(richm2,variable!="rpoB")

richm_background<-richm3[,-6]

#excluding unsorted
richm3<-subset(richm2,variable!="rpoB"&group!="unsorted")

phn_richm<-ggplot(richm2, aes(x=day-69,y=value,group=variable,shape=variable))+
	scale_shape_manual(values=c(0:25,50:70),name="genes")+
	#scale_linetype_manual(values=c(0,1,"13","14","15","24","25","26", "38", "48", "58", "3142","2211"),labels=c("gyrA","phnC"),name="legend title")+
	#geom_line(data=richm_background,colour="grey",alpha=.5,aes(x=day-69,y=value,group=variable),linetype=1)+
	geom_line(data=subset(richm2,group=="one"),size=1.2,alpha=0.3)+
	geom_point(data=subset(richm2,group=="one"),size=3.5)+
	geom_line(data=subset(richm2,group=="two"),size=1.2,alpha=0.3)+
	geom_point(data=subset(richm2,group=="two"),size=3.5)+
	geom_line(data=subset(richm2,group=="no_similarity"),size=1.2,alpha=0.3)+
	geom_point(data=subset(richm2,group=="no_similarity"),size=3.5)+
	geom_line(data=subset(richm2,group=="four"),size=1.2,alpha=0.2)+
	geom_point(data=subset(richm2,group=="four"),size=3.5)+
	facet_wrap(~group,nrow=2,ncol=2)+
	theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	#ggtitle(goi)+
	#scale_fill_discrete(name = "genes")+
	#scale_colour_discrete(name = "data series")+
	#scale_x_discrete(drop=TRUE)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.text=element_text(size=11))+
	xlab("Days after glyphosate addition")+
	#theme(legend.position="none")+
	#theme(strip.text = element_text(face = "italic"))+
	#guides(size=FALSE)+
	ylab('Gene richness in relation to day 0 [%]')
phn_richm
ggsave(file="phn_richness_groups.png", width=14, height=8)


count_function<-function(n){
ggtitle(paste("group",n))
}
pl <- lapply(seq_len(4), count_function)
library(gridExtra)
out<- by(data = richm2, INDICES = richm2$group, FUN = function(m) {
      m <- droplevels(m)
      m <- ggplot(m, aes(x=day-69,y=value,group=variable,shape=variable))+
	scale_shape_manual(values=c(15:20,21:25,11:13),name="Genes")+
	#scale_linetype_manual(values=c(0,1,"13","14","15","24","25","26", "38", "48", "58", "3142","2211"),labels=c("gyrA","phnC"),name="legend title")+
	#geom_line(data=richm_background,colour="grey",alpha=.5,aes(x=day-69,y=value,group=variable),linetype=1)+
	geom_line(size=1.2,alpha=0.3)+
	geom_point(size=3.5)+
	#facet_wrap(~group,nrow=2,ncol=2)+
	theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	coord_cartesian(ylim = c(70, 170),xlim=c(0,70))+
	#scale_fill_discrete(name = "genes")+
	#scale_colour_discrete(name = "data series")+
	#scale_x_discrete(drop=TRUE)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title=element_blank())+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	#theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.text=element_text(size=11))+
	guides(colour=FALSE, size=FALSE, width=FALSE,
		   shape=guide_legend(ncol=1,
						keyheight=1.5,
						label.theme=element_text(size=10,
												face="italic",
												angle=0),
						(title = NULL)))+
	#xlab("Days after glyphosate addition")+
	lapply(seq_len(4), count_function)
	#theme(legend.position="none")+
	#theme(strip.text = element_text(face = "italic"))+
	#guides(size=FALSE)+
	#ylab('Gene richness in relation to day 0 [%]')
   })
do.call(grid.arrange, out)

g<-arrangeGrob(do.call(grid.arrange, out))

ggsave(file="phn_richness_supplement.png",g, width=20, height=8.75)


# The inner arrangeGrob() function arranges the four plots, the main title, 
#   and the global y-axis title.
# The outer grid.arrange() function arranges and draws the arrangeGrob object and the legend.
g2<-grid.arrange(arrangeGrob(g, 
                         #top = textGrob("Main Title", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
                         bottom = textGrob("Days after glyphosate addition", rot = 0, vjust = 1),
                         left = textGrob('Gene richness in relation to day 0 [%]', rot = 90, vjust = 1)))
						 
ggsave(file="phn_richness_supplement2.png",g2, width=20, height=8.75)
################################################################
richm4<-read.csv(file.choose(),row.names=1,sep=";")


phn_richm_mean<-ggplot(richm4, aes(x=day-69,y=value,linetype=variable))+
	geom_line(size=1.2,alpha=1)+
	geom_point(size=3.5)+
	theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	#ggtitle(goi)+
	#scale_discrete(name = "genes")+
	#scale_colour_discrete(name = "data series")+
	#scale_x_discrete(drop=TRUE)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.text=element_text(size=11))+
	xlab("Days after glyphosate addition")+
	theme(legend.position="none")+
	#theme(strip.text = element_text(face = "italic"))+
	#guides(size=FALSE)+
	ylab('Gene richness in relation to day 0 [%]')
phn_richm_mean
ggsave(file="phn_richness_mean.png", width=10, height=9)