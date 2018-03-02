#for newcomer plot
#merge relativ and absolut abundances of newcomer otus
newcomer_abs<-read.csv(file.choose(),sep=";",row.names=1)
water_otus_max_rel_abu<-read.csv(file.choose(),sep=";",row.names=1)
max_abus_newcomer<-merge(newcomer_abs,water_otus_max_rel_abu,by="row.names",all.x=TRUE)
row.names(max_abus_newcomer)<-max_abus_newcomer$Row.names
max_abus_newcomer<-max_abus_newcomer[,-1]
head(max_abus_newcomer)
write.csv(max_abus_newcomer,file="max_abus_newcomer.csv")

#base plot
hist(max_abus_newcomer$maxdna_rel)
hist_dna<-subset(max_abus_newcomer,max_abus_newcomer$maxdna_rel > 0)
nrow(hist_dna)
#90
hist(hist_dna$maxdna_rel)


hist(max_abus_newcomer$maxcdna_rel)
hist_cdna<-subset(max_abus_newcomer,max_abus_newcomer$maxcdna_rel > 0)
nrow(hist_cdna)
#117
hist(hist_cdna$maxcdna_rel)

#ggplot2 ??
ggplot(max_abus_newcomer)+
	geom_bar(aes(x="",y=maxdna_rel))