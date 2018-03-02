###bubble_plot 16S changes	
change_tab<-read.csv(file.choose(),sep=";")

change_tab$radius<-sqrt(change_tab$change/pi) #creates a new variable (radius) from the abundance data which is better suited for visualisation in bubble plots

bubble_change<-ggplot(change_tab,aes(x=method,y=otu,fill=otu))+
	geom_point(aes(size=radius_log*7.5),shape=21)+
	scale_size_identity()+
	theme(legend.position="none")+
	theme(axis.text.x=element_text(size=16))+
	theme(axis.text.y=element_text(size=16))+
	facet_wrap(~nucleic_acid)
bubble_change


test<-change_tab
test[]<-apply(test[,2],2,function(x) ifelse(x<200&x>20,"e",x))


flevels<-levels(bubblemyt$y)
bubblemyt+scale_y_discrete(limits=rev(flevels(bubblemyt$y)))	#seems only to work with alpha characters

ggsave("bubble_myt.tiff", width=35, height=30, dpi=150)

myt2<-read.csv(file.choose(), sep=";")
myt2$radius<-sqrt(myt2$z/pi)

mytall<-ggplot(myt2,aes(x,y))+
	geom_point(aes(size=radius*7.5),shape=21)+
	scale_size_identity()+
	theme_classic()+
	theme(axis.text.x=element_text(size=16, angle=45))+
	theme(axis.text.y=element_text(size=16))

ggsave("mytall3.svg",width=38, height=34)

ggplot(myt,aes(x,y))+
	geom_point(aes(size=radius*7.5),shape=21)+
	scale_size_identity()+
	theme_grey()+
	theme(axis.text.x=element_text(size=16, angle=45))+
	theme(axis.text.y=element_text(size=16))


myt.legend<-read.csv(file.choose(), sep=";")
myt.legend$radius<-sqrt(myt.legend$z/pi)

myt_legend<-ggplot(myt.legend,aes(x,y))+
	geom_point(aes(size=radius*7.5),shape=21)+
	scale_size_identity()+
	theme_classic()

ggsave("myt_legend3.svg",width=38, height=34)