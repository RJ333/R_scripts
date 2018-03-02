#meta-subset erstellen
meta_water_cdna<-subset(meta_data,habitat=="water"&nucleic_acid=="cdna"&days>=45)
head(meta_water_cdna)
head(water_cdna_subset_nonzero)

#turning "days" and "time" into factors (after creating envfit-file)
meta_water_cdna$days<-as.factor(meta_water_cdna$days)
meta_water_cdna$time<-as.factor(meta_water_cdna$time)

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