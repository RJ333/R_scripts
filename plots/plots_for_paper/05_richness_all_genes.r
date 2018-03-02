##for all genes passing the threshold, selection in Excel
######################################################## data preparation ###################################################
richness3<-read.csv(file.choose(),row.names=1,sep=";") #reading table with all selected genes per sample in wide format
meta_omics<-read.csv(file.choose(),sep=";") #reading meta data for samples
trichness3<-t(richness3) #transpose for merging
trichness_meta<-merge(trichness3,meta_omics,by.x="row.names",by.y="X",all.x=TRUE) #merge meta data with richness counts
names(trichness_meta) #check number of columns
trichness_meta<-trichness_meta[,-c(1600:1608)] #removing superfluous meta data columns
names(trichness_meta)
library(reshape2)
richness_melt<-melt(trichness_meta,id=c("Row.names","time","day","new_day")) #melt into long format for plotting, using meta data as id columns
write.csv(richness_melt,file="richness_melt_all_genes.csv") #save melt file and 
#add two factors the different groups (group: group1, group2_sensu_lato and group2_sensu_stricto) (group2:group1, group2)
richness_melt2<-read.csv(file.choose(),row.names=1,sep=";") 
head(richness_melt2)
richness_melt_summary<-summarySE(richness_melt4,measurevar="value",groupvars=c("new_day","group","group2")) #creating statistical data based on script summarySE
#richness_melt_background<-richness_melt2[,-c(7,8)]
#richness_melt_background<-subset(richness_melt_background,variable!="pknB")
#richness_melt_background<-subset(richness_melt_background,variable!="lrsA")
richness_melt4<-subset(richness_melt2,variable!="pknB") #removing outsiders which stretch the plot area (they were still included in summarySE step)
richness_melt4<-subset(richness_melt4,variable!="lrsA") #removing outsiders which stretch the plot area (they were still included in summarySE step)


############################################################ clean plots ###################################################

##plot with all genes in detail for supplement
phn_richness_melt_all<-ggplot(richness_melt4, aes(x=day-69,y=value,colour=group))+	
	geom_line(data=subset(richness_melt4,group=="group1"),aes(group=variable),size=1.5,alpha=0.1)+
	geom_line(data=subset(richness_melt2,group=="group2_sensu_lato"),aes(group=variable),size=1.5,alpha=0.4)+
	geom_line(data=subset(richness_melt2,group=="group2_sensu_stricto"),aes(group=variable),size=1.5,alpha=0.4)+
	geom_ribbon(data=subset(richness_melt_summary,group2=="group1"), 
		aes(x=new_day,ymax=value+ci,ymin=value-ci,colour=group),fill="black",alpha=0.7,size=1)+
	geom_ribbon(data=subset(richness_melt_summary,group2!="group1"), 
		aes(x=new_day,ymax=value+ci,ymin=value-ci,group=group),fill="grey30",colour=NA,alpha=0.6,size=1)+
	facet_wrap(~group2,ncol=2)+
	theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size=20,face="bold"))+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=17))+
	theme(legend.position="none")+
	theme(strip.background = element_blank(),
		strip.text.x = element_blank())+
	xlab("Days after glyphosate addition")+
	ylab('Gene richness in relation to day 0 [%]')
phn_richness_melt_all
ggsave(file="SI_fig_02_detailed_richness_group_mean_with_ci.png", width=20, height=8.75)



#plot with averaged groups and confidence interval
plot_richness_melt_mean_ci<-ggplot(richness_melt_summary, aes(x=new_day,y=value,group=group,linetype=group))+
geom_line(size=2)+
geom_ribbon(aes(ymax=value+ci,ymin=value-ci,group=group),alpha=0.1)+
theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size=20,face="bold"))+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=17))+
	theme(legend.position="none")+
	theme(strip.background = element_blank(),
		strip.text.x = element_blank())+
	xlab("Days after glyphosate addition")+
	ylab('Gene richness in relation to day 0 [%]')
plot_richness_melt_mean_ci
ggsave(file="fig_05_richness_confidence_intervall_95.png", width=10, height=8)


################################################################# devel ########################################################
##plot with all genes in detail for supplement
phn_richness_melt_all<-ggplot(richness_melt4, aes(x=day-69,y=value,colour=group))+	
	#scale_linetype_manual(values=c(0,1,"13","14","15","24","25","26", "38", "48", "58", "3142","2211"),labels=c("gyrA","phnC"),name="legend title")+
	#geom_line(data=richness_melt_background,colour="grey",alpha=.2,aes(x=day-69,y=value,group=variable),linetype=1,size=1.5)+
	geom_line(data=subset(richness_melt4,group=="group1"),aes(group=variable),size=1.5,alpha=0.1)+
	#geom_point(data=subset(richness_melt4,group=="group1"),aes(group=variable),size=2.0)+
	geom_line(data=subset(richness_melt2,group=="group2_sensu_lato"),aes(group=variable),size=1.5,alpha=0.4)+
	#geom_point(data=subset(richness_melt2,group=="group2_sensu_lato"),aes(group=variable),size=2.0)+
	geom_line(data=subset(richness_melt2,group=="group2_sensu_stricto"),aes(group=variable),size=1.5,alpha=0.4)+
	#geom_point(data=subset(richness_melt2,group=="group2_sensu_stricto"),aes(group=variable),size=2.0)+
	#stat_summary(data=subset(richness_melt2,group=="group2_sensu_lato"),aes(group=group),colour="black",fun.y="mean", geom="line",size=2,alpha=1)+
	#stat_summary(data=subset(richness_melt2,group=="group2_sensu_stricto"),aes(group=group),colour="black",fun.y="mean", geom="line", size=2,alpha=1)+
	#stat_summary(data=subset(richness_melt2,group=="group1"),aes(group=group),colour="black",fun.y="mean", geom="line", size=2,alpha=1)+
	geom_ribbon(data=subset(richness_melt_summary,group2=="group1"), aes(x=new_day,ymax=value+ci,ymin=value-ci,colour=group),fill="grey20",alpha=0.7,size=1)+
	geom_ribbon(data=subset(richness_melt_summary,group2!="group1"), aes(x=new_day,ymax=value+ci,ymin=value-ci,group=group,colour=group),fill="grey20",alpha=0.7,size=1)+
	#geom_ribbon(data=subset(richness_melt_summary,group=="group2_sensu_lato"), aes(x=new_day,ymax=value+ci,ymin=value-ci,fill=group),alpha=0.3)+
	#stat_summary(data=subset(richness_melt2,group=="three"),aes(group=group),colour="black",fun.y="mean", geom="line", size=1.5,alpha=1)+
	#scale_linetype_manual("mean",values=c("mean1"=4,"mean2"=3,"mean3"=1),labels=c("group2_sensu_lato","group2_sensu_stricto","group1"))+
	facet_wrap(~group2,nrow=2,ncol=2)+
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
	guides(linetype=guide_legend(override.aes=list(size=1.2)))+
	ylab('Gene richness in relation to day 0 [%]')
phn_richness_melt_all
ggsave(file="detailed_richness_group_mean_with_background_and_ci.png", width=14, height=8)


######gene_richness plots only mean with statistics
#confidence interval
plot_richness_melt_mean_ci<-ggplot(richness_melt_summary, aes(x=new_day,y=value,group=group,linetype=group))+
geom_line()+
geom_ribbon(aes(ymax=value+ci,ymin=value-ci,group=group),alpha=0.1)
plot_richness_melt_mean_ci
ggsave(file="richness_confidence_intervall_95.png", width=14, height=8)
#standard error
plot_richness_melt_mean_se<-ggplot(richness_melt_summary, aes(x=new_day,y=value,group=group,linetype=group))+
geom_line()+
geom_ribbon(aes(ymax=value+se,ymin=value-se,group=group),alpha=0.1)
plot_richness_melt_mean_se
ggsave(file="richness_confidence_standard_error.png", width=14, height=8)
#standard deviation
plot_richness_melt_mean_sd<-ggplot(richness_melt_summary, aes(x=new_day,y=value,group=group,linetype=group))+
geom_line()+
geom_ribbon(aes(ymax=value+sd,ymin=value-sd,group=group),alpha=0.3)
plot_richness_melt_mean_sd
ggsave(file="richness_confidence_standard_deviation.png", width=14, height=8)