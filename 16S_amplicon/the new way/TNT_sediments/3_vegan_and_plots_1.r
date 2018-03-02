#vegan analysis and plots
##for vegan
##otus in spalten/proben als rownames 
head(cast_tnt_otu1_table3)
cast_tnt_otu1_table3_dna<-subset(cast_tnt_otu1_table3,grepl("dna",row.names))

#diversity analysis (pielou in excel)
h_cast_tnt_otu1_table3<-diversity(cast_tnt_otu1_table3)
s_cast_tnt_otu1_table3<-specnumber(cast_tnt_otu1_table3)
write.csv(h_cast_tnt_otu1_table3,file="h_cast_tnt_otu1_table3.csv")
write.csv(s_cast_tnt_otu1_table3,file="s_cast_tnt_otu1_table3.csv")

#2 dimensionen
nmds_tnt1_notrans<-metaMDS(cast_tnt_otu1_table3,try=100,autotransform=FALSE)
nmds_tnt1_trans<-metaMDS(cast_tnt_otu1_table3,try=100,autotransform=TRUE) #sehr merkwürdige ergebnisse
plot(nmds_tnt1_notrans,main="NMDS TNT sediments untransformed", type="t",display="sites")
plot(nmds_tnt1_trans,main="NMDS TNT sediments untransformed", type="t",display="sites",xlim=c(-297,-299),ylim=c(-0.5,0.5)) #sehr merkwürdige ergebnisse
#extract data
names(nmds_water_dna)
nmds_water_dna_coord<-as.data.frame(nmds_water_dna$points)
str(nmds_water_dna_coord)
row.names(nmds_water_dna_coord)<-row.names(cast_tnt_otu1_table3)
plot(nmds_water_dna_coord)
text(nmds_water_dna_coord,labels=row.names(nmds_water_dna_coord),cex=0.5) #gives the sample ids, must be called when the plot is already open
#meta_file<-read.csv(file.choose(),sep=";",row.names=1)
meta_nmds_water_dna<-merge(nmds_water_dna_coord,meta_file,by="row.names",all.x=TRUE)
head(meta_nmds_water_dna)

#make "row.names" to actual row.names
row.names(meta_nmds_water_dna)<-meta_nmds_water_dna$Row.names
meta_nmds_water_dna<-meta_nmds_water_dna[,-1]
head(meta_nmds_water_dna)

#removing nmds coords for envfit
meta_nmds_water_dna2<-meta_nmds_water_dna[,3:10]
str(meta_nmds_water_dna2)
#keine cellcounts, kein treatment, keine parallele, kein time, kein habitat
meta_nmds_water_dna3<-meta_nmds_water_dna2[,c(2,7,8)]
str(meta_nmds_water_dna3)

#turning "days" and "time" into factors (after creating envfit-file)
meta_nmds_water_dna$days<-as.factor(meta_nmds_water_dna$days)
meta_nmds_water_dna$time<-as.factor(meta_nmds_water_dna$time)

#plotting
colfunc <- colorRampPalette(c("yellow2","red","purple"))
colfunc(16)
plot(nmds_water_dna,main="NMDS_water_dna", type="p",display="sites",xlim=c(-1,1),ylim=c(-1,1))
with(meta_nmds_water_dna, points(nmds_water_dna,col=colfunc(16)[meta_nmds_water_dna$days],pch=c(16,17)[meta_nmds_water_dna2$treatment]))
with(meta_nmds_water_dna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_nmds_water_dna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
water_dna_envfit<-envfit(nmds_water_dna,meta_nmds_water_dna3,permutations=100)
plot(water_dna_envfit)
text(nmds_water_dna_coord,labels=row.names(nmds_water_dna_coord),cex=0.5)

#########oder

#plotting
colfunc <- colorRampPalette(c("lightgrey","yellow","blue"))
colfunc(16)
plot(nmds_water_dna,main="NMDS_water_dna", type="t",display="species",cex=0.4)
with(meta_nmds_water_dna, points(nmds_water_dna,col=colfunc(16)[meta_nmds_water_dna$days],pch=c(16,17)[meta_nmds_water_dna2$treatment]))
with(meta_nmds_water_dna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_nmds_water_dna, legend("bottomleft", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
water_dna_envfit<-envfit(nmds_water_dna,meta_nmds_water_dna3,permutations=100)
plot(water_dna_envfit)
text(nmds_water_dna_coord,labels=row.names(nmds_water_dna_coord),cex=0.7)

#anosim 
anosim_treatment_water_dna<-anosim(cast_tnt_otu1_table3,meta_nmds_water_dna$treatment,permutations=999,distance="bray")

plot(anosim_treatment)

#cca
str(meta_water_cdna)
meta_water_cdna$days<-as.factor(meta_water_cdna$days)
str(meta_water_cdna)
cca_water_cdna<-cca(water_cdna_subset_nonzero~disturbance + days,data=meta_water_cdna)

colfunc <- colorRampPalette(c("yellow2","red","purple"))
colfunc(16)
plot(cca_water_cdna,main="CCA_water_cdna", type="p",display="sites")
with(meta_water_cdna, points(cca_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(cca_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_cdna,labels=row.names(meta_water_cdna),cex=0.7)

names(cca_water_cdna)
summary(cca_water_cdna)

#dca
dca_water_cdna<-decorana(water_cdna_subset_nonzero)
plot(dca_water_cdna,type="p",display="sites",main="dca_water_cdna")
with(meta_water_cdna, points(dca_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(dca_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_cdna,labels=row.names(meta_water_cdna),cex=0.7)

names(dca_water_cdna)
summary(dca_water_cdna)

#rda
rda_water_cdna<-rda(water_cdna_subset_nonzero~days+disturbance,meta_water_cdna)
plot(rda_water_cdna,type="p",display="sites",main="rda_water_cdna")
with(meta_water_cdna, points(rda_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("bottomleft", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(rda_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
names(rda_water_cdna)
summary(rda_water_cdna)

#pca
pca_water_cdna<-rda(water_cdna_subset_nonzero)
plot(pca_water_cdna,main="PCA_water_cdna", type="p",display="sites")
with(meta_water_cdna, points(pca_water_cdna,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(pca_water_cdna,groups=interaction(meta_water_cdna$days,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
#text(meta_water_cdna,labels=row.names(meta_water_cdna),cex=0.7)

#pca
pca_water_cdna_scaled<-rda(water_cdna_subset_nonzero,scale = TRUE)
plot(pca_water_cdna_scaled,main="PCA_water_cdna_scaled", type="p",display="sites")
with(meta_water_cdna, points(pca_water_cdna_scaled,col=colfunc(16)[meta_water_cdna$days],pch=c(16,17)[meta_water_cdna$treatment]))
with(meta_water_cdna, legend("bottomright", legend = levels(treatment), bty = "n", pch = c(16,17)))
with(meta_water_cdna, legend("topright", legend = levels(days), pch=16, bty = "n", col=colfunc(16)))
ordihull(pca_water_cdna_scaled,groups=interaction(meta_water_cdna$disturbance,meta_water_cdna$treatment),label=TRUE,cex=0.4,draw="polygon",col="red")
text(pca_water_cdna_scaled,labels=row.names(meta_water_cdna),cex=0.7)