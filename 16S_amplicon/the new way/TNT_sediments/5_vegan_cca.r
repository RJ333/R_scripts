muncomp_mean<-read.csv(file.choose(),sep=";",row.names=1)
#merge tnt concentrations with meta data
head(meta_data_tnt10)
head(tnt10_vegan_cdna)
head(muncomp_mean)
muncomp_mean<-muncomp_mean[,c(1:8)]
head(muncomp_mean)

cdna10_vegan_meta<-merge(meta_data_tnt10,tnt10_vegan_cdna,by="row.names",all.y=TRUE)
head(cdna10_vegan_meta)
names(cdna10_vegan_meta)


cdna10_mun<-merge(muncomp_mean,cdna10_vegan_meta,by="station_ord",all=TRUE)
names(cdna10_mun)
cdna10_all_meta<-cdna10_mun[,c(1:28)]
names(cdna10_all_meta)
row.names(cdna10_all_meta)<-cdna10_all_meta$Row.names													#turns "sample_name" into row names
cdna10_all_meta<-cdna10_all_meta[,-9]
cdna10_all_meta<-cdna10_all_meta[,-2]	

head(tnt10_vegan_cdna,40)
head(cdna10_all_meta,40)
#sort both tables the same way
cdna10_all_meta<-cdna10_all_meta[order(row.names(cdna10_all_meta)),]
#remove Mo7 due to missing values from both tables
cdna10_all_meta2<-subset(cdna10_all_meta,station.y!="Mo7")
head(cdna10_all_meta2,40)
tnt10_vegan_cdna2<-tnt10_vegan_cdna[c(1:36),]
head(tnt10_vegan_cdna2,40)
#table and meta data ready for cca
str(cdna10_all_meta2)

#plot munition concentrations along the area
ggplot(cdna10_all_meta,aes(x=parallel))+
geom_bar(aes(y=sum_all,colour="sum expl"),stat="identity", width=1)+
geom_point(aes(y=TNT,colour="TNT"))+
geom_point(aes(y=x4ADNT,colour="2-ADNT"))+
geom_point(aes(y=x2ADNT,colour="4-ADNT"))+
facet_wrap(~station_ord,nrow=1)

#TNT constrained
tnt_cca_cdna10<-cca(tnt10_vegan_cdna2~TNT,data=cdna10_all_meta2,autotransform=FALSE)

colfunc <- colorRampPalette(c("yellow2","red","green","blue"))
colfunc(4)
plot(tnt_cca_cdna10,main="tnt_cca_cdna10", type="p",display="sites")
with(cdna10_all_meta2, points(tnt_cca_cdna10,col=colfunc(4)[cdna10_all_meta2$cruise],pch=c(16,17,18)[cdna10_all_meta2$parallel]))
with(cdna10_all_meta2, legend("bottomright", legend = levels(as.factor(parallel)), bty = "n", pch = c(16,17,18)))
with(cdna10_all_meta2, legend("topright", legend = levels(cruise), pch=16, bty = "n", col=colfunc(4)))
text(tnt_cca_cdna10,labels=cdna10_all_meta2$station_ord,cex=0.7)

plot(tnt_cca_cdna10, type="p",display="sites")
with(cdna10_all_meta2, points(tnt_cca_cdna10,col=colfunc(4)[cdna10_all_meta2$cruise],pch=c(16,17,18)[cdna10_all_meta2$parallel]))
with(cdna10_all_meta2, legend("bottomright", legend = levels(as.factor(parallel)), bty = "n", pch = c(16,17,18)))
with(cdna10_all_meta2, legend("topright", legend = levels(cruise), pch=16, bty = "n", col=colfunc(4)))
text(tnt_cca_cdna10,labels=cdna10_all_meta2$station_ord,cex=0.5)

#no constrained
tnt_cca_cdna10<-cca(tnt10_vegan_cdna2,autotransform=FALSE)
#no constrained including Mo7
tnt_cca_cdna10<-cca(tnt10_vegan_cdna,autotransform=FALSE)

#sum of compounds
tnt_cca_cdna10<-cca(tnt10_vegan_cdna2~sum_all,data=cdna10_all_meta2,autotransform=FALSE)

#sum of compounds + shannon
tnt_cca_cdna10<-cca(tnt10_vegan_cdna2~sum_all+shannon,data=cdna10_all_meta2,autotransform=FALSE)

#sum of degraded compounds
tnt_cca_cdna10<-cca(tnt10_vegan_cdna2~sum_deg,data=cdna10_all_meta2,autotransform=FALSE)

#ratio TNT to degraded compounds
tnt_cca_cdna10<-cca(tnt10_vegan_cdna2~ratio_deg,data=cdna10_all_meta2,autotransform=FALSE)
