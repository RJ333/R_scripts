#shannon indeces
#otus in spalten/proben als rownames

Examples

data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")
r.2 <- rarefy(BCI, 2)
alpha <- fisher.alpha(BCI)
pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
J <- H/log(S)


#transformed normalized count tables
head(tnorm_water_dna)
h_water_dna<-diversity(tnorm_water_dna)
s_water_dna<-specnumber(tnorm_water_dna)
j_water_dna<-h_water_dna/log(s_water_dna)
write.csv(h_water_dna,file="h_water_dna.csv")
write.csv(s_water_dna,file="s_water_dna.csv")
write.csv(j_water_dna,file="j_water_dna.csv")


#gemittelte Proben
vegan_dna_mean_V2<-read.csv(file.choose(),sep=";",row.names=1)
h_water_dna_mean_V2_paper<-diversity(vegan_dna_mean_V2)
s_water_dna_mean_V2_paper<-specnumber(vegan_dna_mean_V2)
j_water_dna_mean_V2_paper<-h_water_dna_mean_V2_paper/log(s_water_dna_mean_V2_paper)
write.csv(h_water_dna_mean_V2_paper,file="h_water_dna_mean_V2_paper.csv")
write.csv(s_water_dna_mean_V2_paper,file="s_water_dna_mean_V2_paper.csv")
write.csv(j_water_dna_mean_V2_paper,file="j_water_dna_mean_V2_paper.csv")

vegan_cdna_mean_V2<-read.csv(file.choose(),sep=";",row.names=1)
h_water_cdna_mean_V2_paper<-diversity(vegan_cdna_mean_V2)
s_water_cdna_mean_V2_paper<-specnumber(vegan_cdna_mean_V2)
j_water_cdna_mean_V2_paper<-h_water_cdna_mean_V2_paper/log(s_water_cdna_mean_V2_paper)
write.csv(h_water_cdna_mean_V2_paper,file="h_water_cdna_mean_V2_paper.csv")
write.csv(s_water_cdna_mean_V2_paper,file="s_water_cdna_mean_V2_paper.csv")
write.csv(j_water_cdna_mean_V2_paper,file="j_water_cdna_mean_V2_paper.csv")