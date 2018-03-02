##contig_abundance plots


#check for c-p lyase genes
#create gene subset

test1<-table(gsub("\\..*$", "", phnM_select$contig_id))
unique(gsub("\\..*$", "", phnM_select$contig_id))
table(phnM_select$contig_id)

#test<-subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("k141_186291|k141_225608|k141_37311.1|k141_81456|k141_158295", contig_id))
test<-subset(kallisto_prokka_concoct2_metaxa_checkm, contig_id == "k141_186291.1" | grepl("k141_37311.1|k141_186291.0", contig_id))
test2<-test[,c(1,5,6,8,10,12,14,16,18,20,22,24,26)]
library(reshape2)
molten_test<-melt(test2,id=c("contig_id","ec_number","genes"))
#meta data for plotting: treatment/control and time
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
merged_test<-merge(molten_test,meta_omics,by.x="variable",by.y="row.names",all.x=T)
#write.csv(merged_test,file="merged_test.csv")

##all phn genes
#test_plot<-subset(merged_test, grepl("test", genes))
#test_plot<-test_plot[order(test_plot[,4]),]
#nrow(test_plot)

test_bar<-ggplot(merged_test, aes(x = day, y = value, fill=contig_id))+				
	geom_bar(width = 1.7, stat = "identity",position="dodge")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#geom_text(data=subset(merged_test, value > 5),aes(label=contig_id),size=rel(2.7), position = position_stack(vjust = .5))+
	#geom_label_repel(data=subset(merged_test, value > 1),position="stack",min.segment.length = unit(0.1, "lines"),aes(label=genes,size=1.4))+
	ggtitle("phn_genes")+
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
	#theme(legend.position="none")
test_bar
 ggsave(file="phn_genes.jpg", width=14, height=8)
 