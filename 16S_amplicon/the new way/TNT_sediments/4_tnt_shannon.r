##diversity plots

jo<-subset(meta_data_tnt10,area!="lab_control"&nucleic_acid=="cdna")
divtest<-ggplot(jo,aes(x=as.factor(parallel)))+
	#geom_point(aes(y=shannon,colour="Shannon Diversity"),colour="black",size=3)+
	#geom_point(aes(y=species_richness,colour="Species Richness"))+
	geom_point(data=subset(jo,cruise=="L17_07"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="blue") +
	geom_point(data=subset(jo,cruise=="L17_08"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="red")+
	geom_point(data=subset(jo,cruise=="L16_14"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="green")+
	geom_point(data=subset(jo,station=="Mo7"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="black")+
	theme_bw()+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Parallels", y="Shannon Diversity")+	
	facet_wrap(~station_ord*cruise*cruise_date,nrow=1)+
	theme(legend.position='none')
divtest	
ggsave(file="shannon_diversity_otu10_cdna.png", width=20, height=6)


 
  geom_point(data=subset(jo,cruise=="L17_07"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=2.5,colour="blue") +
  geom_point(data=subset(jo,cruise=="L17_08"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="red")+
  geom_point(data=subset(jo,cruise=="L16_14"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="green")+
  geom_point(data=subset(jo,station=="Mo7"),aes(y=shannon,shape=as.factor(parallel),colour=cruise),size=4,colour="yellow")+
  geom_text(data=meta_nmds_tnt10_cdna,aes(x=MDS1,y=MDS2,label=station,size=3,vjust=-1))+
  coord_equal()